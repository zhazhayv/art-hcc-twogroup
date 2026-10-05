# -*- coding: utf-8 -*-
"""Copy Swiss/CTD/GeneCards/STRING sources into locked_inputs and freeze MMP9 probability."""
from __future__ import annotations

import hashlib
import json
import shutil
from datetime import date
from pathlib import Path

import pandas as pd

from constants import ARTHCC_NEWDATA, EXPERIMENTAL3, FULL67, LOCKED, NETWORK6, STRICT14


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_art_root() -> Path:
    from constants import SOURCES, ROOT
    if SOURCES.exists():
        for d in SOURCES.iterdir():
            if d.is_dir() and "Art&HCC" in d.name:
                return d
        for d in SOURCES.iterdir():
            if d.is_dir() and "\u809d\u764c" in d.name:
                return d
    desk = Path(r"C:\Users\20917\Desktop")
    for d in desk.iterdir():
        if d.is_dir() and "Art&HCC" in d.name:
            return d
    raise FileNotFoundError("Art&HCC folder not found under art&hcc/01_network_pharmacology_source")


def main():
    root = find_art_root()
    target_xlsx = next(root.rglob("14*/2.*/*.xlsx"))
    disease_xlsx = next(root.rglob("14*/3.*/*.xlsx"))
    string14_src = next(p for p in root.rglob("14*/6.*/string_interactions_short_14_nodes.tsv"))
    string67_src = ARTHCC_NEWDATA / "Full67_STRING_raw.tsv"
    if not string67_src.exists():
        raise FileNotFoundError(string67_src)

    swiss = pd.read_excel(target_xlsx, "SWISS raw")
    swiss_conv = pd.read_excel(target_xlsx, "SWISS conversion")
    ctd = pd.read_excel(target_xlsx, "CTD filtered")
    gc_raw = pd.read_excel(disease_xlsx, header=None)
    header_row = gc_raw.index[gc_raw[0] == "Symbol"][0]
    gc = gc_raw.iloc[header_row + 1 :].copy()
    gc.columns = ["Symbol", "Name", "Type", "Relevance_Score", "Knowledge"]
    gc = gc.dropna(subset=["Symbol"])
    gc = gc[gc["Symbol"].astype(str).str.lower() != "nan"]
    gc["Relevance_Score"] = pd.to_numeric(gc["Relevance_Score"], errors="coerce")
    gc = gc.dropna(subset=["Relevance_Score"])

    mmp9 = swiss.loc[swiss["Common name"].astype(str).str.upper() == "MMP9"]
    if mmp9.empty:
        raise RuntimeError("MMP9 missing from SWISS raw")
    mmp9_p = float(mmp9.iloc[0]["Probability*"])

    swiss_out = LOCKED / "swiss_raw.csv"
    conv_out = LOCKED / "swiss_conversion.csv"
    ctd_out = LOCKED / "ctd_filtered.csv"
    gc_out = LOCKED / "genecards_hcc.csv"
    s14_out = LOCKED / "Strict14_STRING_raw.tsv"
    s67_out = LOCKED / "Full67_STRING_raw.tsv"
    swiss.to_csv(swiss_out, index=False)
    swiss_conv.to_csv(conv_out, index=False)
    ctd.to_csv(ctd_out, index=False)
    gc.to_csv(gc_out, index=False)
    shutil.copy2(string14_src, s14_out)
    shutil.copy2(string67_src, s67_out)

    n67 = pd.read_csv(s67_out, sep="\t")
    n_edges = len(n67)
    if n_edges != 660:
        raise RuntimeError(f"Full_67 STRING must be the 660-edge export, got {n_edges}")

    gc_lookup = dict(zip(gc["Symbol"].astype(str), gc["Relevance_Score"]))
    swiss_lookup = {
        str(r["Common name"]).upper(): float(r["Probability*"])
        for _, r in swiss.iterrows()
        if pd.notna(r["Common name"])
    }
    ctd_genes = sorted({str(x).upper() for x in ctd["Gene Symbol"].dropna()})

    miss = []
    aliases = {"GBA": ["GBA", "GBA1"], "MMP9": ["MMP9"], "MGEA5": ["MGEA5", "OGA"]}
    for gene, names in aliases.items():
        gc_hits = {n: gc_lookup.get(n) for n in names if n in gc_lookup}
        miss.append(
            {
                "gene": gene,
                "in_swiss": any(n in swiss_lookup for n in names),
                "swiss_probability": next((swiss_lookup[n] for n in names if n in swiss_lookup), None),
                "in_ctd": any(n in ctd_genes for n in names),
                "genecards_symbols_present": gc_hits,
                "genecards_relevance": max(gc_hits.values()) if gc_hits else None,
            }
        )

    manifest = {
        "lock_date": str(date.today()),
        "source_target_xlsx": str(target_xlsx),
        "source_disease_xlsx": str(disease_xlsx),
        "source_string14": str(string14_src),
        "source_string67": str(string67_src),
        "Full67_STRING_edges": n_edges,
        "STRING_536_edge_snapshot": "never used",
        "MMP9_swiss_probability": mmp9_p,
        "MMP9_note": "Single number from SWISS raw Common name==MMP9 Probability*; methods and results must reuse this field.",
        "swiss_n_rows": int(len(swiss)),
        "ctd_n_rows": int(len(ctd)),
        "ctd_genes": ctd_genes,
        "genecards_n": int(len(gc)),
        "genecards_score_gt2": int((gc["Relevance_Score"] > 2).sum()),
        "network6": NETWORK6,
        "experimental3": EXPERIMENTAL3,
        "strict14": STRICT14,
        "full67": FULL67,
        "experimental_miss_preview": miss,
        "ctd_interaction_typing": "Harmonizome CTD gene-chemical association export; no binding/expression/inferred split. Network-group CTD labels are untyped.",
        "sha256": {p.name: sha256(p) for p in [swiss_out, conv_out, ctd_out, gc_out, s14_out, s67_out]},
    }
    (LOCKED / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("MMP9_swiss_probability", mmp9_p)
    print("STRING edges", n_edges)
    print("locked", LOCKED)
    for row in miss:
        print("miss", row)


if __name__ == "__main__":
    main()
