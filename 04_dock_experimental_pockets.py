# -*- coding: utf-8 -*-
"""Site-directed docking of ART/DHA on experimental-group pockets; redock RMSD < 2 A gate."""
from __future__ import annotations

import json
import math
import re
import shutil
import ssl
import subprocess
from pathlib import Path

import numpy as np

from constants import DOCK, HA, P0_BCL2, ARTHCC_NEWDATA, RESULTS, SEEDS, VINA

SSL_CTX = ssl._create_unverified_context()
SKIP_HET = {
    "HOH", "WAT", "H2O", "SO4", "PO4", "GOL", "EDO", "PEG", "NAG", "MAN", "BMA",
    "FUC", "NDG", "ACT", "CL", "NA", "ZN", "CA", "MG", "DMS", "EPE", "TRS",
}

TARGETS = [
    # gene, pdb, preferred ligand residue, literature note
    ("GBA", "2V3E", "IFG", "Human GBA catalytic domain with iminosugar; align to Chen et al. GBA active-site residues"),
    ("MMP9", "1GKC", None, "Human MMP9 catalytic domain with hydroxamate inhibitor"),
    ("MGEA5", "5VVV", None, "Human OGA catalytic domain with active-site ligand; OGA=MGEA5"),
    ("TRYPSIN_NEG", "3PTB", "BEN", "Unrelated protein negative control (trypsin/benzamidine)"),
]


def http_get(url: str) -> bytes:
    import urllib.request
    req = urllib.request.Request(url, headers={"User-Agent": "art-hcc-twogroup/1.0"})
    with urllib.request.urlopen(req, timeout=120, context=SSL_CTX) as resp:
        return resp.read()


def dominant_het(pdb_text: str, prefer: str | None):
    counts = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        rn = line[17:20].strip()
        if rn in SKIP_HET or rn == "UNX":
            continue
        counts[rn] = counts.get(rn, 0) + 1
    if prefer and prefer in counts:
        return prefer, counts
    if not counts:
        return None, counts
    return max(counts, key=counts.get), counts


def write_receptor_ligand(pdb_text: str, ligcode: str, rec_path: Path, lig_path: Path):
    lig_all = []
    rec_all = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            rec_all.append(line)
        elif line.startswith("HETATM") and line[17:20].strip() == ligcode:
            lig_all.append(line)
    if not lig_all:
        rec_path.write_text("\n".join(rec_all) + "\nEND\n", encoding="ascii")
        lig_path.write_text("END\n", encoding="ascii")
        return np.zeros((0, 3))
    chain = lig_all[0][21]
    resi = lig_all[0][22:26]
    lig = [ln for ln in lig_all if ln[21] == chain and ln[22:26] == resi]
    rec = [ln for ln in rec_all if ln[21] == chain] or rec_all
    rec_path.write_text("\n".join(rec) + "\nEND\n", encoding="ascii")
    lig_path.write_text("\n".join(lig) + "\nEND\n", encoding="ascii")
    xyz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in lig], float)
    return xyz


def box_from_xyz(xyz: np.ndarray, pad=8.0):
    mn, mx = xyz.min(0), xyz.max(0)
    center = (mn + mx) / 2
    size = np.maximum(mx - mn + pad, 16.0)
    size = np.minimum(size, 30.0)
    return tuple(np.round(center, 2)), tuple(np.round(size, 1))


def prepare_pdbqt(src_pdb: Path, out_pdbqt: Path, is_lig: bool):
    py = Path(r"C:\Users\20917\AppData\Local\Programs\Python\Python310\python.exe")
    scripts = py.parent / "Scripts"
    if is_lig:
        exe = scripts / "mk_prepare_ligand.exe"
        cmd = [str(exe), "-i", str(src_pdb), "-o", str(out_pdbqt)]
    else:
        exe = scripts / "mk_prepare_receptor.exe"
        cmd = [str(exe), "-p", str(src_pdb), "-o", str(out_pdbqt)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not out_pdbqt.exists():
        # Open Babel fallback
        ob = scripts / "obabel.exe"
        extra = ["-xh"] if is_lig else ["-xr"]
        proc2 = subprocess.run(
            [str(ob), str(src_pdb), "-O", str(out_pdbqt), *extra],
            capture_output=True, text=True,
        )
        if not out_pdbqt.exists():
            raise RuntimeError(f"pdbqt failed {src_pdb}: {proc.stderr}\n{proc2.stderr}")
    if is_lig:
        text = out_pdbqt.read_text(encoding="utf-8", errors="replace")
        if text.count("ROOT") > 1:
            # keep first molecule only (Vina 1.1.2 cannot parse concatenated pdbqt)
            keep = []
            seen_root = 0
            for line in text.splitlines():
                if line.startswith("ROOT"):
                    seen_root += 1
                    if seen_root > 1:
                        break
                keep.append(line)
            out_pdbqt.write_text("\n".join(keep) + "\n", encoding="ascii")
    return proc.returncode


def parse_pdbqt_models(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    models, cur, score = [], [], None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            cur, score = [], None
        elif line.startswith("REMARK VINA RESULT:"):
            score = float(line.split()[3])
        elif line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elem = line[77:79].strip() if len(line) >= 79 else ""
            if not elem:
                elem = re.sub(r"[0-9]", "", name)[:1] or "C"
            if elem.upper().startswith("H"):
                continue
            cur.append((elem.upper()[0], xyz, name))
        elif line.startswith("ENDMDL"):
            models.append({"score": score, "atoms": cur})
    if not models and cur:
        models.append({"score": score, "atoms": cur})
    return models


def parse_pdb_het(path: Path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        name = line[12:16].strip()
        elem = line[76:78].strip() if len(line) >= 78 else re.sub(r"[0-9]", "", name)[:1]
        if elem.upper().startswith("H"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        atoms.append((elem.upper()[0], xyz, name))
    return atoms


def hungarian_inplace_rmsd(a, b):
    from collections import Counter
    ea, eb = [x[0] for x in a], [x[0] for x in b]
    Pa, Qa = [], []
    try:
        from scipy.optimize import linear_sum_assignment
    except Exception:
        linear_sum_assignment = None
    for e in sorted(set(ea) & set(eb)):
        A = [x[1] for x in a if x[0] == e]
        B = [x[1] for x in b if x[0] == e]
        n = min(len(A), len(B))
        if n == 0:
            continue
        if linear_sum_assignment is None:
            used = set()
            for p in A[:n]:
                best = None
                for i, q in enumerate(B):
                    if i in used:
                        continue
                    d = np.linalg.norm(p - q)
                    if best is None or d < best[0]:
                        best = (d, i, q)
                used.add(best[1])
                Pa.append(p)
                Qa.append(best[2])
        else:
            cost = np.zeros((len(A), len(B)))
            for i, p in enumerate(A):
                for j, q in enumerate(B):
                    cost[i, j] = np.linalg.norm(p - q)
            r, c = linear_sum_assignment(cost)
            for i, j in zip(r, c):
                Pa.append(A[i])
                Qa.append(B[j])
    if len(Pa) < 3:
        return None
    P = np.asarray(Pa)
    Q = np.asarray(Qa)
    return float(np.sqrt(((P - Q) ** 2).sum() / len(P)))


def getbestrms(a, b):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign
    except Exception:
        return None

    def to_mol(atoms):
        xyz = "\n".join(
            f"{e} {p[0]:.3f} {p[1]:.3f} {p[2]:.3f}" for e, p, _ in atoms
        )
        mol = Chem.MolFromXYZBlock(f"{len(atoms)}\n\n{xyz}")
        return mol

    ma, mb = to_mol(a), to_mol(b)
    if ma is None or mb is None:
        return None
    try:
        return float(rdMolAlign.GetBestRMS(ma, mb))
    except Exception:
        return None


def write_conf(path, receptor, ligand, center, size, seed, out_name):
    path.write_text(
        f"receptor = {receptor}\nligand = {ligand}\n"
        f"center_x = {center[0]}\ncenter_y = {center[1]}\ncenter_z = {center[2]}\n"
        f"size_x = {size[0]}\nsize_y = {size[1]}\nsize_z = {size[2]}\n"
        f"exhaustiveness = 32\nnum_modes = 9\nseed = {seed}\n"
        f"out = {out_name}\n",
        encoding="utf-8",
    )


def run_vina(work: Path, conf: Path):
    log = work / (conf.stem + "_run.log")
    with open(log, "w", encoding="utf-8") as fh:
        proc = subprocess.run(
            [str(VINA), "--config", conf.name],
            cwd=str(work), stdout=fh, stderr=subprocess.STDOUT, text=True,
        )
    return proc.returncode


def copy_ligands(work: Path):
    art = P0_BCL2 / "Artesunate.pdbqt"
    dha = P0_BCL2 / "DHA.pdbqt"
    shutil.copy2(art, work / "Artesunate.pdbqt")
    if dha.exists():
        shutil.copy2(dha, work / "DHA.pdbqt")
    else:
        shutil.copy2(ARTHCC_NEWDATA / "DHA_pubchem.sdf", work / "DHA.sdf")


def second_engine_recheck(work: Path, center, size):
    """GNINA/smina if present; else Vina Python API; else record unavailable."""
    for exe_name in ("gnina", "smina"):
        exe = shutil.which(exe_name)
        if not exe:
            continue
        out = work / f"{exe_name}_art_recheck.pdbqt"
        cmd = [
            exe, "--receptor", str(work / "receptor.pdbqt"),
            "--ligand", str(work / "Artesunate.pdbqt"),
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--exhaustiveness", "16", "--out", str(out),
        ]
        log = work / f"{exe_name}_recheck.log"
        with open(log, "w", encoding="utf-8") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
        return {"engine": exe_name, "returncode": proc.returncode, "out": str(out)}
    try:
        from vina import Vina
        v = Vina(sf_name="vina")
        v.set_receptor(str(work / "receptor.pdbqt"))
        v.set_ligand_from_file(str(work / "Artesunate.pdbqt"))
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=8, n_poses=5)
        score = float(v.energies()[0][0])
        return {"engine": "python_vina_1_2", "score": score}
    except Exception as exc:
        return {"engine": None, "note": f"GNINA/smina/python-vina unavailable: {exc}"}


def process_target(gene, pdb_id, prefer, note):
    work = DOCK / f"{gene}_{pdb_id}"
    work.mkdir(parents=True, exist_ok=True)
    pdb_path = work / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path.write_bytes(http_get(f"https://files.rcsb.org/download/{pdb_id}.pdb"))
    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    ligcode, counts = dominant_het(text, prefer)
    rec = {"gene": gene, "pdb": pdb_id, "note": note, "het_counts": counts, "crystal_ligand": ligcode}
    if not ligcode:
        rec["status"] = "no_organic_ligand"
        rec["main_figure"] = False
        return rec
    xyz = write_receptor_ligand(text, ligcode, work / "receptor.pdb", work / "crystal_ligand.pdb")
    if len(xyz) < 5:
        rec["status"] = "ligand_too_small"
        rec["main_figure"] = False
        return rec
    center, size = box_from_xyz(xyz, pad=8.0)
    rec["box_center"] = center
    rec["box_size"] = size
    prepare_pdbqt(work / "receptor.pdb", work / "receptor.pdbqt", is_lig=False)
    prepare_pdbqt(work / "crystal_ligand.pdb", work / "crystal_ligand.pdbqt", is_lig=True)
    copy_ligands(work)

    # fail-fast redock seed 0
    seed0 = SEEDS[0]
    conf = work / f"vina_redock_ex32_s{seed0}.conf"
    write_conf(conf, "receptor.pdbqt", "crystal_ligand.pdbqt", center, size, seed0, f"redock_ex32_s{seed0}_out.pdbqt")
    run_vina(work, conf)
    models = parse_pdbqt_models(work / f"redock_ex32_s{seed0}_out.pdbqt")
    crystal = parse_pdb_het(work / "crystal_ligand.pdb")
    rmsd_ip = hungarian_inplace_rmsd(models[0]["atoms"], crystal) if models else None
    rmsd_gb = None
    try:
        rmsd_gb = getbestrms(models[0]["atoms"], crystal) if models else None
    except Exception:
        rmsd_gb = None
    rec["redock_seed0_score"] = models[0]["score"] if models else None
    rec["redock_rmsd_inplace"] = rmsd_ip
    rec["redock_rmsd_getbestrms"] = rmsd_gb
    passed = (rmsd_ip is not None and rmsd_ip < 2.0) or (rmsd_gb is not None and rmsd_gb < 2.0)
    rec["pass_2A"] = bool(passed)
    rec["main_figure"] = bool(passed) and gene != "TRYPSIN_NEG"
    if not passed:
        rec["status"] = "redock_fail"
        rec["art_dha"] = "not_run_pocket_failed_gate"
        return rec

    rec["status"] = "redock_pass"
    rec["scores"] = {}
    for lig in ("crystal_ligand", "Artesunate", "DHA"):
        if not (work / f"{lig}.pdbqt").exists() and lig != "crystal_ligand":
            continue
        ligfile = "crystal_ligand.pdbqt" if lig == "crystal_ligand" else f"{lig}.pdbqt"
        scores = []
        for seed in SEEDS:
            out = f"{lig}_ex32_s{seed}_out.pdbqt"
            conf = work / f"vina_{lig}_ex32_s{seed}.conf"
            write_conf(conf, "receptor.pdbqt", ligfile, center, size, seed, out)
            run_vina(work, conf)
            ms = parse_pdbqt_models(work / out)
            if ms:
                scores.append(ms[0]["score"])
        rec["scores"][lig] = scores
    rec["second_engine"] = second_engine_recheck(work, center, size)
    return rec


def write_md_stub(passed):
    md = DOCK / "md"
    md.mkdir(exist_ok=True)
    if not passed:
        (md / "MD_SECTION_OMITTED.txt").write_text(
            "All new experimental-group pockets failed in-place redocking RMSD < 2.0 A.\n"
            "Per protocol the MD section is omitted from the manuscript.\n"
            "BCL2/4LVT remains the only protocol-valid pocket already reported.\n",
            encoding="utf-8",
        )
        return
    (md / "README_50ns.md").write_text(
        "\n".join(
            [
                "# Conditional 50 ns MD",
                "Run only on the first redock-pass experimental pocket.",
                "Ligands: artesunate and DHA; 3 replicates x 50 ns.",
                "Report ligand RMSD mean +/- SD and pocket occupancy.",
                "Use CHARMM-GUI or the existing HA_docking GROMACS workflow.",
                "This session does not start 300 ns of production MD on the workstation.",
            ]
        ),
        encoding="utf-8",
    )


def main():
    rows = []
    for gene, pdb_id, prefer, note in TARGETS:
        print("TARGET", gene, pdb_id)
        try:
            rec = process_target(gene, pdb_id, prefer, note)
        except Exception as exc:
            rec = {"gene": gene, "pdb": pdb_id, "status": "error", "error": str(exc), "main_figure": False, "pass_2A": False}
        rows.append(rec)
        print(json.dumps({k: rec[k] for k in rec if k != "het_counts"}, default=str))
    (RESULTS / "docking_experimental_pockets.json").write_text(json.dumps(rows, indent=2, default=str), encoding="utf-8")
    passed = [r for r in rows if r.get("pass_2A") and r.get("gene") != "TRYPSIN_NEG"]
    write_md_stub(passed)
    (RESULTS / "docking_gate.json").write_text(
        json.dumps(
            {
                "n_passed_experimental": len(passed),
                "md_section": "include" if passed else "omit",
                "bcl2_4lvt": {
                    "reuse": True,
                    "GetBestRMS": 1.54,
                    "ART": -7.5,
                    "DHA": -7.9,
                    "source": "HA_docking/_p0_tmp/BCL2_4LVT",
                },
                "CASP3_PARP1": "supplement only",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
