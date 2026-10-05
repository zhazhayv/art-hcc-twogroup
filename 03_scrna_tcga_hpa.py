# -*- coding: utf-8 -*-
"""Nine-gene cell-type and bulk validation: GSE149614, GSE125449, TCGA-LIHC, HPA, Cox (supplement)."""
from __future__ import annotations

import gzip
import io
import json
import ssl
import tarfile
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata

from constants import EXPERIMENTAL3, FIGURES, NINE, NINE_EXPR, NETWORK6, RESULTS

SSL_CTX = ssl._create_unverified_context()
UA = {"User-Agent": "Mozilla/5.0 art-hcc-twogroup/1.0"}
CACHE = RESULTS / "cohort_cache"
CACHE.mkdir(parents=True, exist_ok=True)

CLASS_RULES = [
    ("Hepatocyte_cholangiocyte", ("hepato", "cholangio", "epithelial", "malignant", "tumor", "hcc", "bi-potent", "bipotent", "liver cell")),
    ("T_NK", ("t cell", "tcell", "nk", "nkt", "cd8", "cd4", "treg", "ilc")),
    ("Myeloid", ("macro", "mono", "dc", "dendritic", "neutrophil", "kupffer", "myeloid", "tam", "mast")),
    ("B", ("b cell", "bcell", "plasma", "b_cell")),
    ("Endothelial", ("endoth", "lsec", "vec")),
    ("Fibroblast", ("fibro", "stellate", "caf", "hsc", "mesench", "pericyte", "smooth muscle", "myofibro")),
]


def http_get(url, dest: Path | None = None, timeout=300):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
        data = resp.read()
    if dest:
        dest.write_bytes(data)
    return data


def canonical_gene(sym: str) -> str:
    s = str(sym).upper()
    if s in ("GBA", "GBA1"):
        return "GBA"
    if s in ("MGEA5", "OGA"):
        return "MGEA5"
    return s


def map_class(label: str) -> str:
    s = str(label).lower()
    for cls, keys in CLASS_RULES:
        if any(k in s for k in keys):
            return cls
    return "Other"


def try_urls(urls):
    last = None
    for url in urls:
        try:
            print("GET", url[:120])
            return http_get(url)
        except Exception as exc:
            last = exc
            print(" fail", exc)
    raise RuntimeError(last)


def load_tisch_like():
    """TISCH2 static files if present; otherwise GEO supplementary / processed tables."""
    out = {}
    tisch_meta = [
        "https://tisch.comp-genomics.org/static/data/LIHC_GSE149614/LIHC_GSE149614_CellMetainfo_table.tsv",
        "https://tisch.comp-genomics.org/static/data/GSE149614/GSE149614_CellMetainfo_table.tsv",
        "https://tisch.comp-genomics.org/static/data/LIHC_GSE149614/LIHC_GSE149614_cellinfo.csv",
    ]
    tisch_expr = [
        "https://tisch.comp-genomics.org/static/data/LIHC_GSE149614/LIHC_GSE149614_Expression.txt",
        "https://tisch.comp-genomics.org/static/data/GSE149614/GSE149614_Expression.txt",
    ]
    meta_path = CACHE / "GSE149614_meta.tsv"
    expr_path = CACHE / "GSE149614_expr_9genes.tsv"
    if not meta_path.exists():
        try:
            meta_path.write_bytes(try_urls(tisch_meta))
        except Exception:
            meta_path.write_text("", encoding="utf-8")
    return meta_path, expr_path, tisch_expr


def download_geo_soft(acc: str) -> Path:
    dest = CACHE / f"{acc}_family.soft.gz"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz"
    try:
        http_get(url, dest)
    except Exception as exc:
        dest.write_text(str(exc), encoding="utf-8")
    return dest


def parse_soft_sample_titles(path: Path):
    rows = []
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            sample, title, source = None, None, None
            for line in fh:
                if line.startswith("^SAMPLE"):
                    if sample:
                        rows.append({"sample": sample, "title": title, "source": source})
                    sample = line.split("=")[-1].strip()
                    title = source = None
                elif line.startswith("!Sample_title"):
                    title = line.split("=")[-1].strip()
                elif line.startswith("!Sample_source_name"):
                    source = line.split("=")[-1].strip()
            if sample:
                rows.append({"sample": sample, "title": title, "source": source})
    except Exception:
        pass
    return pd.DataFrame(rows)


def hpa_json(gene: str):
    cache = CACHE / f"hpa_{gene}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    ens = {
        "PARP1": "ENSG00000143799",
        "BCL2": "ENSG00000171791",
        "CASP3": "ENSG00000164305",
        "BAX": "ENSG00000087088",
        "FAS": "ENSG00000026103",
        "CASP9": "ENSG00000132906",
        "GBA": "ENSG00000177628",
        "GBA1": "ENSG00000177628",
        "MMP9": "ENSG00000100985",
        "MGEA5": "ENSG00000198408",
        "OGA": "ENSG00000198408",
    }
    url = f"https://www.proteinatlas.org/{ens[gene]}.json"
    try:
        rec = json.loads(http_get(url).decode("utf-8"))
    except Exception as exc:
        rec = {"error": str(exc)}
    cache.write_text(json.dumps(rec), encoding="utf-8")
    time.sleep(0.2)
    return rec


def hpa_ihc_liver(gene: str):
    rec = hpa_json(gene)
    path = rec.get("ihc") or rec.get("antibody") or []
    # pathology / tissue
    tissues = rec.get("tissueExpression") or rec.get("tissue") or []
    ihc_entries = []
    def walk(obj, acc):
        if isinstance(obj, dict):
            name = str(obj.get("name") or obj.get("tissue") or obj.get("organ") or "")
            if "liver" in name.lower() or "hepat" in name.lower():
                acc.append(
                    {
                        "gene": gene,
                        "tissue": name,
                        "level": obj.get("level") or obj.get("expression") or obj.get("ihc"),
                        "keys": list(obj.keys())[:12],
                    }
                )
            for v in obj.values():
                walk(v, acc)
        elif isinstance(obj, list):
            for v in obj:
                walk(v, acc)
    walk(rec, ihc_entries)
    # RNA tissue liver
    rna_liver = None
    for item in rec.get("rnaTissueSpecificity") or rec.get("rnaTissue") or []:
        if isinstance(item, dict) and "liver" in str(item).lower():
            rna_liver = item
            break
    return ihc_entries[:8], rec.get("error")


def download_xena_lihc():
    """UCSC Xena: gene expression (log2+1 fpkm) and phenotype for TCGA-LIHC."""
    expr_path = CACHE / "TCGA-LIHC.htseq_fpkm.tsv.gz"
    pheno_path = CACHE / "TCGA-LIHC.GDC_phenotype.tsv.gz"
    surv_path = CACHE / "TCGA-LIHC.survival.tsv.gz"
    urls = {
        expr_path: "https://gdc.xenahubs.net/download/TCGA-LIHC.htseq_fpkm.tsv.gz",
        pheno_path: "https://gdc.xenahubs.net/download/TCGA-LIHC.GDC_phenotype.tsv.gz",
        surv_path: "https://gdc.xenahubs.net/download/TCGA-LIHC.survival.tsv.gz",
    }
    # Smaller: use TOIL RSEM gene expected_count? htseq_fpkm is large (~50MB). Prefer HiSeq gene.
    alt = {
        expr_path: "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FHiSeqV2.gz",
        pheno_path: "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix.gz",
        surv_path: "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLIHC_survival.txt.gz",
    }
    for dest, url in alt.items():
        if dest.exists() and dest.stat().st_size > 10000:
            continue
        try:
            print("Xena", url)
            http_get(url, dest, timeout=600)
        except Exception as e:
            print("xena fail", dest.name, e)
    return expr_path, pheno_path, surv_path


def read_xena_expr_9genes(path: Path) -> pd.DataFrame:
    # first column gene, rest samples; file may be gz
    open_fn = gzip.open if str(path).endswith(".gz") else open
    wanted = set(NINE_EXPR)
    rows = []
    header = None
    with open_fn(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            gene = line.split("\t", 1)[0]
            if gene in wanted:
                parts = line.rstrip("\n").split("\t")
                rows.append(parts)
                if len(rows) == len(wanted):
                    break
    if not rows:
        raise RuntimeError(f"No locked genes in {path}")
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(header[0]).astype(float)
    df.index = [canonical_gene(i) for i in df.index]
    df = df.groupby(level=0).mean()
    return df


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n)
    prev = 1.0
    for i, idx in enumerate(order[::-1], start=0):
        rank = n - i
        val = min(prev, p[idx] * n / rank)
        ranked[idx] = val
        prev = val
    return ranked


def tcga_tumor_normal(expr: pd.DataFrame):
    samples = list(expr.columns)
    rec = []
    for gene in expr.index:
        tumor, normal = [], []
        for s in samples:
            code = s[13:15] if len(s) >= 15 and s.startswith("TCGA") else ""
            val = expr.loc[gene, s]
            if code == "01" or (len(s) >= 15 and s[13:15] == "01"):
                tumor.append(val)
            elif code == "11" or (len(s) >= 15 and s[13:15] == "11"):
                normal.append(val)
            else:
                # Xena HiSeqV2 uses TCGA barcodes
                if "-11" in s[12:16]:
                    normal.append(val)
                elif "-01" in s[12:16]:
                    tumor.append(val)
        # fallback by 14th-15th
        if not tumor and not normal:
            for s in samples:
                bits = s.split("-")
                if len(bits) >= 4 and bits[3].startswith("11"):
                    normal.append(expr.loc[gene, s])
                elif len(bits) >= 4 and bits[3].startswith("01"):
                    tumor.append(expr.loc[gene, s])
        if len(tumor) < 3 or len(normal) < 3:
            rec.append({"gene": gene, "n_tumor": len(tumor), "n_normal": len(normal), "p": np.nan, "delta_mean": np.nan})
            continue
        stat, p = mannwhitneyu(tumor, normal, alternative="two-sided")
        rec.append(
            {
                "gene": gene,
                "group": "network" if gene in NETWORK6 else "experimental",
                "n_tumor": len(tumor),
                "n_normal": len(normal),
                "mean_tumor": float(np.mean(tumor)),
                "mean_normal": float(np.mean(normal)),
                "delta_mean": float(np.mean(tumor) - np.mean(normal)),
                "p": float(p),
            }
        )
    df = pd.DataFrame(rec)
    df["fdr"] = bh_fdr(df["p"].fillna(1).values)
    return df


def cox_os(expr: pd.DataFrame, surv_path: Path):
    """Univariable Cox via partial likelihood with Newton; same model for both groups."""
    if not surv_path.exists() or surv_path.stat().st_size < 100:
        return pd.DataFrame()
    open_fn = gzip.open if str(surv_path).endswith(".gz") else open
    try:
        surv = pd.read_csv(open_fn(surv_path, "rt", encoding="utf-8", errors="replace"), sep="\t")
    except Exception:
        return pd.DataFrame()
    # standardize columns
    cols = {c.lower(): c for c in surv.columns}
    idcol = surv.columns[0]
    time_col = None
    event_col = None
    for k, c in cols.items():
        if k in ("os.time", "os_time", "_os_time", "os.t"):
            time_col = c
        if k in ("os", "os.event", "_os_event", "os_event"):
            event_col = c
    if time_col is None:
        for c in surv.columns:
            if "time" in c.lower() and "os" in c.lower():
                time_col = c
    if event_col is None:
        for c in surv.columns:
            if c.lower() in ("os", "event", "status"):
                event_col = c
    if time_col is None or event_col is None:
        return pd.DataFrame({"note": ["survival file parsed but OS columns not found"], "columns": [",".join(surv.columns[:12])]})

    def cox_one(x, t, e):
        x = np.asarray(x, float)
        t = np.asarray(t, float)
        e = np.asarray(e, float)
        mask = np.isfinite(x) & np.isfinite(t) & np.isfinite(e) & (t > 0)
        x, t, e = x[mask], t[mask], e[mask]
        if e.sum() < 10:
            return np.nan, np.nan, np.nan
        x = (x - x.mean()) / (x.std() + 1e-8)
        beta = 0.0
        for _ in range(25):
            ll1 = 0.0
            ll2 = 0.0
            order = np.argsort(-t)
            x_o, e_o = x[order], e[order]
            risk_sum = 0.0
            risk_x = 0.0
            risk_x2 = 0.0
            # process from shortest remaining? Use Efron-unspecified Breslow from longest time
            # Standard: sort time ascending, at each death sum over R(t)
            order = np.argsort(t)
            x, t, e = x[order], t[order], e[order]
            n = len(x)
            expx = np.exp(beta * x)
            for i in range(n):
                if e[i] != 1:
                    continue
                R = expx[i:]
                S0 = R.sum()
                S1 = (x[i:] * R).sum()
                S2 = ((x[i:] ** 2) * R).sum()
                ll1 += x[i] - S1 / S0
                ll2 += -(S2 / S0 - (S1 / S0) ** 2)
            if abs(ll2) < 1e-12:
                break
            step = ll1 / (-ll2) if ll2 < 0 else ll1
            beta = beta + 0.5 * step
        se = np.sqrt(1 / max(-ll2, 1e-8))
        z = beta / se
        # two-sided normal
        from math import erfc, sqrt
        p = erfc(abs(z) / sqrt(2))
        return float(beta), float(se), float(p)

    # align samples: expr columns vs survival ids
    surv[idcol] = surv[idcol].astype(str)
    rec = []
    for gene in expr.index:
        xs, ts, es = [], [], []
        for _, row in surv.iterrows():
            sid = str(row[idcol])
            # match 15-char barcode to expr
            hits = [c for c in expr.columns if c.startswith(sid[:12]) and (len(c) < 14 or c[13:15] == "01")]
            if not hits:
                hits = [c for c in expr.columns if c[:12] == sid[:12]]
            if not hits:
                continue
            xs.append(float(expr.loc[gene, hits[0]]))
            ts.append(float(row[time_col]))
            es.append(float(row[event_col]))
        beta, se, p = cox_one(xs, ts, es) if xs else (np.nan, np.nan, np.nan)
        rec.append(
            {
                "gene": gene,
                "group": "network" if gene in NETWORK6 else "experimental",
                "cox_beta_per_sd": beta,
                "cox_se": se,
                "cox_p": p,
                "n": len(xs),
                "model": "univariable Cox OS, same specification for both groups",
            }
        )
    df = pd.DataFrame(rec)
    if "cox_p" in df.columns:
        df["cox_fdr"] = bh_fdr(df["cox_p"].fillna(1).values)
    return df


def fetch_cellxgene_or_hpa_scrna():
    """
    Prefer TISCH/CELLxGENE; fall back to HPA RNA single-cell type resource
    (not a substitute for GSE accessions if those download).
    """
    hpa_sc = CACHE / "rna_single_cell_type.tsv.zip"
    if not hpa_sc.exists():
        try:
            http_get("https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip", hpa_sc, timeout=300)
        except Exception as e:
            print("HPA sc fail", e)
    return hpa_sc


def parse_hpa_sc(path: Path):
    if not path.exists() or path.stat().st_size < 1000:
        return pd.DataFrame()
    with zipfile.ZipFile(path) as z:
        name = z.namelist()[0]
        with z.open(name) as fh:
            df = pd.read_csv(fh, sep="\t")
    gene_col = [c for c in df.columns if "gene" in c.lower() and "name" in c.lower()]
    gene_col = gene_col[0] if gene_col else df.columns[1]
    tissue_col = [c for c in df.columns if "tissue" in c.lower() or "cluster" in c.lower() or "cell" in c.lower()]
    ntp = [c for c in df.columns if "ntpm" in c.lower() or "ptpm" in c.lower() or "nTPM" in c]
    if not ntp:
        ntp = [df.columns[-1]]
    sub = df[df[gene_col].isin(NINE_EXPR)].copy()
    sub["major_class"] = sub[tissue_col[0]].map(map_class) if tissue_col else "Other"
    return sub


def try_tisch_gene_avg(acc_tag: str):
    """Download TISCH 'gene.in.celltype' style tables for the 9 genes if the static API exists."""
    rows = []
    bases = [
        f"https://tisch.comp-genomics.org/static/data/{acc_tag}/",
        f"https://tisch.comp-genomics.org/gallery/?cancer=LIHC&dataset={acc_tag}",
    ]
    meta_urls = [
        f"https://tisch.comp-genomics.org/static/data/{acc_tag}/{acc_tag}_CellMetainfo_table.tsv",
        f"https://tisch.comp-genomics.org/static/data/{acc_tag}/{acc_tag}_metainfo.tsv",
        f"https://tisch.comp-genomics.org/static/data/{acc_tag}/{acc_tag}_cellinfo.csv",
    ]
    meta_file = CACHE / f"{acc_tag}_meta.tsv"
    if not meta_file.exists():
        for url in meta_urls:
            try:
                http_get(url, meta_file)
                break
            except Exception:
                continue
    return meta_file


def peak_from_patient_celltype(df):
    """df columns: patient, major_class, gene, mean_expr, detection_rate"""
    peaks = []
    for gene, g in df.groupby("gene"):
        # majority of patients' peak class
        patient_peak = []
        for pat, pg in g.groupby("patient"):
            idx = pg["mean_expr"].idxmax()
            patient_peak.append(pg.loc[idx, "major_class"])
        mode = Counter(patient_peak).most_common(1)[0][0] if patient_peak else "NA"
        agree = Counter(patient_peak).most_common(1)[0][1] / len(patient_peak) if patient_peak else np.nan
        peaks.append(
            {
                "gene": gene,
                "group": "network" if gene in NETWORK6 else "experimental",
                "peak_class_mode": mode,
                "patients_with_that_peak": int(Counter(patient_peak).most_common(1)[0][1] if patient_peak else 0),
                "n_patients": len(patient_peak),
                "fraction_patients_agree": agree,
                "patient_peaks": ",".join(patient_peak),
            }
        )
    return pd.DataFrame(peaks)


def build_pseudobulk_from_long(long_df):
    return peak_from_patient_celltype(long_df)


def stop_rule(peaks_a: pd.DataFrame, peaks_b: pd.DataFrame | None):
    def mode_group(peaks, genes):
        sub = peaks[peaks["gene"].isin(genes)]
        if sub.empty:
            return None
        return Counter(sub["peak_class_mode"]).most_common(1)[0][0]

    net_a = mode_group(peaks_a, NETWORK6)
    exp_a = mode_group(peaks_a, EXPERIMENTAL3)
    separated_a = bool(net_a and exp_a and net_a != exp_a)
    separated_b = None
    if peaks_b is not None and not peaks_b.empty:
        net_b = mode_group(peaks_b, NETWORK6)
        exp_b = mode_group(peaks_b, EXPERIMENTAL3)
        separated_b = bool(net_b and exp_b and net_b != exp_b)
        direction = separated_a and separated_b and net_a == net_b and exp_a == exp_b
        # plan: majority peak class differs AND GSE125449 direction consistent
        ok = separated_a and separated_b
    else:
        ok = False
        direction = False
    triggered_stop = not ok
    return {
        "GSE149614_network_mode_class": net_a,
        "GSE149614_experimental_mode_class": exp_a,
        "GSE149614_separated": separated_a,
        "GSE125449_separated": separated_b,
        "stop_rule_triggered": triggered_stop,
        "continue_celltype_story": (not triggered_stop),
        "rule": "Main-cohort mode peak class of network6 != experimental3 AND GSE125449 also separated. Else stop at threshold grid + miss path.",
    }


def plot_dotplot(long_df, title, path):
    if long_df.empty:
        return
    genes = [g for g in NINE if g in set(long_df["gene"])]
    classes = [c for c in ["Hepatocyte_cholangiocyte", "T_NK", "Myeloid", "B", "Endothelial", "Fibroblast", "Other"] if c in set(long_df["major_class"])]
    # patient-averaged
    avg = long_df.groupby(["gene", "major_class"], as_index=False).agg(mean_expr=("mean_expr", "mean"), detection_rate=("detection_rate", "mean"))
    fig, ax = plt.subplots(figsize=(7.4, 5.2), dpi=300)
    for i, gene in enumerate(genes):
        for j, cls in enumerate(classes):
            hit = avg[(avg.gene == gene) & (avg.major_class == cls)]
            if hit.empty:
                continue
            size = 40 + 180 * float(hit.detection_rate.iloc[0] if pd.notna(hit.detection_rate.iloc[0]) else 0)
            val = float(hit.mean_expr.iloc[0])
            ax.scatter(j, i, s=size, c=[val], cmap="viridis", vmin=0, vmax=max(avg.mean_expr.max(), 1e-6), edgecolors="k", linewidths=0.3)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def synthesize_from_hpa_if_needed(hpa_sc_df: pd.DataFrame, tag: str):
    """Map HPA liver-related clusters into the six major classes as a documented fallback."""
    if hpa_sc_df.empty:
        return pd.DataFrame()
    gene_col = [c for c in hpa_sc_df.columns if str(c).lower() in ("gene name", "genename", "gene")]
    # already filtered
    rec = []
    # identify columns
    gcol = "Gene name" if "Gene name" in hpa_sc_df.columns else hpa_sc_df.columns[1]
    ccol = [c for c in hpa_sc_df.columns if "Cell type" in c or "cell type" in c.lower() or "Cluster" in c]
    ccol = ccol[0] if ccol else hpa_sc_df.columns[2]
    vcol = [c for c in hpa_sc_df.columns if "nTPM" in c or "pTPM" in c or "NX" in c]
    vcol = vcol[0] if vcol else hpa_sc_df.columns[-1]
    # restrict to liver if a tissue column exists
    tissue_cols = [c for c in hpa_sc_df.columns if "tissue" in c.lower()]
    df = hpa_sc_df
    if tissue_cols:
        df = df[df[tissue_cols[0]].astype(str).str.contains("liver|Liver", case=False, na=False)]
        if df.empty:
            df = hpa_sc_df
    for _, r in df.iterrows():
        gene = str(r[gcol])
        if gene not in NINE:
            continue
        rec.append(
            {
                "dataset": tag,
                "patient": "HPA_reference",
                "major_class": map_class(str(r[ccol])),
                "gene": gene,
                "mean_expr": float(r[vcol]) if pd.notna(r[vcol]) else 0.0,
                "detection_rate": np.nan,
                "source_note": "HPA RNA single-cell type fallback; not GSE patient pseudobulk",
            }
        )
    return pd.DataFrame(rec)


def try_cellxgene_census():
    try:
        import cellxgene_census  # noqa: F401
    except Exception:
        return pd.DataFrame()
    return pd.DataFrame()


def plot_tcga_box(expr, stats, path):
    fig, axes = plt.subplots(3, 3, figsize=(10, 8), dpi=300)
    axes = axes.ravel()
    for ax, gene in zip(axes, NINE):
        if gene not in expr.index:
            ax.axis("off")
            continue
        tumor, normal = [], []
        for s in expr.columns:
            bits = s.split("-")
            val = expr.loc[gene, s]
            if len(bits) >= 4 and bits[3].startswith("11"):
                normal.append(val)
            elif len(bits) >= 4 and bits[3].startswith("01"):
                tumor.append(val)
        ax.boxplot([normal, tumor], tick_labels=["NT", "TP"], showfliers=False)
        row = stats[stats.gene == gene]
        ttl = gene
        if len(row):
            ttl += f" FDR={row.iloc[0].fdr:.2e}"
        ax.set_title(ttl, fontsize=9)
    fig.suptitle("TCGA-LIHC tumor vs adjacent (Wilcoxon, BH)")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def stream_nine_genes_gz(url: str, dest_tsv: Path, gene_col_index=0):
    """Download a gzip gene x cell matrix and keep only locked genes without writing the full matrix."""
    if dest_tsv.exists() and dest_tsv.stat().st_size > 100:
        return dest_tsv
    req = urllib.request.Request(url, headers=UA)
    wanted = set(NINE_EXPR)
    found = {}
    tmp = dest_tsv.with_suffix(dest_tsv.suffix + ".partial")
    try:
        with urllib.request.urlopen(req, timeout=600, context=SSL_CTX) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                header = gz.readline().decode("utf-8", errors="replace").rstrip("\n")
                tmp.write_text(header + "\n", encoding="utf-8")
                for raw in gz:
                    line = raw.decode("utf-8", errors="replace")
                    gene = line.split("\t", 1)[0].strip().strip('"')
                    gene = gene.split(".")[0]
                    if gene in wanted and gene not in found:
                        found[gene] = True
                        with tmp.open("a", encoding="utf-8") as fh:
                            fh.write(line if line.endswith("\n") else line + "\n")
                    if len(found) == len(wanted):
                        break
        tmp.replace(dest_tsv)
    except Exception as exc:
        dest_tsv.write_text(f"ERROR\t{exc}\n", encoding="utf-8")
    return dest_tsv


def metadata_to_long(meta_path: Path, expr_path: Path, dataset: str):
    if not meta_path.exists() or not expr_path.exists():
        return pd.DataFrame()
    if expr_path.read_text(encoding="utf-8", errors="replace")[:5] == "ERROR":
        return pd.DataFrame()
    try:
        meta = pd.read_csv(meta_path, sep="\t")
    except Exception:
        try:
            meta = pd.read_csv(meta_path, sep=",")
        except Exception:
            return pd.DataFrame()
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    # cell id column
    cell_col = None
    for c in meta.columns:
        if str(c).lower() in {"cell", "barcode", "cellid", "cell_id", "cell.name"}:
            cell_col = c
            break
    if cell_col is None:
        cell_col = meta.columns[0]
    patient_col = next((c for c in meta.columns if str(c).lower() in {"patient", "sample", "donor", "orig.ident", "orig_ident", "patient_id"}), None)
    type_col = next(
        (
            c
            for c in meta.columns
            if any(k in str(c).lower() for k in ("celltype", "cell_type", "annotation", "cluster", "lineage", "major"))
        ),
        None,
    )
    if type_col is None:
        type_col = meta.columns[min(1, len(meta.columns) - 1)]
    meta = meta.copy()
    meta["_cell"] = meta[cell_col].astype(str)
    meta["_patient"] = meta[patient_col].astype(str) if patient_col else "unknown"
    meta["_class"] = meta[type_col].map(map_class)
    expr.columns = [str(c) for c in expr.columns]
    rec = []
    cells = list(expr.columns)
    meta_idx = meta.set_index("_cell")
    for gene in expr.index:
        gname = str(gene)
        if gname not in NINE_EXPR:
            continue
        vals = expr.loc[gname]
        grouped = defaultdict(list)
        det = defaultdict(lambda: [0, 0])
        for cell in cells:
            if cell not in meta_idx.index:
                continue
            row = meta_idx.loc[cell]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            key = (row["_patient"], row["_class"])
            v = float(vals[cell])
            grouped[key].append(v)
            det[key][1] += 1
            if v > 0:
                det[key][0] += 1
        for (pat, cls), arr in grouped.items():
            rec.append(
                {
                    "dataset": dataset,
                    "patient": pat,
                    "major_class": cls,
                    "gene": canonical_gene(gname),
                    "mean_expr": float(np.mean(arr)),
                    "detection_rate": det[(pat, cls)][0] / max(det[(pat, cls)][1], 1),
                    "n_cells": len(arr),
                    "source_note": dataset,
                }
            )
    return pd.DataFrame(rec)


def download_gse149614():
    meta = CACHE / "GSE149614_HCC.metadata.updated.txt.gz"
    expr9 = CACHE / "GSE149614_nine_genes.tsv"
    urls_meta = [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614_HCC.metadata.updated.txt.gz",
    ]
    urls_count = [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
    ]
    if not meta.exists() or meta.stat().st_size < 100:
        for url in urls_meta:
            try:
                http_get(url, meta, timeout=300)
                break
            except Exception as e:
                print("meta fail", e)
    if meta.exists() and str(meta).endswith(".gz"):
        out_meta = CACHE / "GSE149614_metadata.tsv"
        if not out_meta.exists():
            try:
                with gzip.open(meta, "rt", encoding="utf-8", errors="replace") as fh:
                    out_meta.write_text(fh.read(), encoding="utf-8")
            except Exception:
                out_meta = meta
    else:
        out_meta = meta
    for url in urls_count:
        if expr9.exists() and expr9.stat().st_size > 100 and not expr9.read_text(encoding="utf-8", errors="replace").startswith("ERROR"):
            break
        stream_nine_genes_gz(url, expr9)
        if expr9.exists() and not expr9.read_text(encoding="utf-8", errors="replace").startswith("ERROR"):
            break
    return out_meta, expr9


def download_gse125449_tisch():
    tag = "LIHC_GSE125449_aPDL1aCTLA4"
    meta = CACHE / f"{tag}_CellMetainfo_table.tsv"
    avg = CACHE / f"{tag}_expression_majorlineage.txt"
    urls_meta = [
        f"https://tisch.comp-genomics.org/static/data/{tag}/{tag}_CellMetainfo_table.tsv",
    ]
    urls_avg = [
        f"https://tisch.comp-genomics.org/static/data/{tag}/{tag}_expression_Celltype_majorlineage.txt",
        f"https://tisch.comp-genomics.org/static/data/{tag}/{tag}_expression.txt",
    ]
    if not meta.exists():
        for url in urls_meta:
            try:
                http_get(url, meta)
                break
            except Exception as e:
                print("tisch meta", e)
    if not avg.exists():
        for url in urls_avg:
            try:
                http_get(url, avg)
                break
            except Exception as e:
                print("tisch expr", e)
    return meta, avg


def tisch_avg_to_long(avg_path: Path, dataset: str):
    if not avg_path.exists() or avg_path.stat().st_size < 50:
        return pd.DataFrame()
    try:
        df = pd.read_csv(avg_path, sep="\t", index_col=0)
    except Exception:
        return pd.DataFrame()
    rec = []
    for gene in NINE_EXPR:
        if gene not in df.index:
            continue
        for cls in df.columns:
            rec.append(
                {
                    "dataset": dataset,
                    "patient": "cohort_average",
                    "major_class": map_class(str(cls)),
                    "gene": canonical_gene(gene),
                    "mean_expr": float(df.loc[gene, cls]),
                    "detection_rate": np.nan,
                    "source_note": "TISCH2 major-lineage average; not patient-level",
                }
            )
    return pd.DataFrame(rec)


def main():
    notes = []
    # GEO SOFT for sample inventory (always)
    soft149 = download_geo_soft("GSE149614")
    soft125 = download_geo_soft("GSE125449")
    inv149 = parse_soft_sample_titles(soft149)
    inv125 = parse_soft_sample_titles(soft125)
    inv149.to_csv(RESULTS / "GSE149614_GEO_sample_inventory.csv", index=False)
    inv125.to_csv(RESULTS / "GSE125449_GEO_sample_inventory.csv", index=False)
    notes.append(f"GSE149614 GEO samples parsed: {len(inv149)}")
    notes.append(f"GSE125449 GEO samples parsed: {len(inv125)}")

    long149 = pd.DataFrame()
    long125 = pd.DataFrame()
    try:
        meta149, expr149 = download_gse149614()
        notes.append(f"GSE149614 meta={meta149} expr={expr149} bytes={expr149.stat().st_size if expr149.exists() else 0}")
        long149 = metadata_to_long(meta149, expr149, "GSE149614")
    except Exception as exc:
        notes.append(f"GSE149614 download/parse failed: {exc}")

    try:
        meta125, avg125 = download_gse125449_tisch()
        long125 = tisch_avg_to_long(avg125, "GSE125449_TISCH")
        if long125.empty:
            notes.append("GSE125449 TISCH average empty")
        else:
            notes.append(f"GSE125449 TISCH rows={len(long125)}")
    except Exception as exc:
        notes.append(f"GSE125449 TISCH failed: {exc}")

    hpa_sc_zip = fetch_cellxgene_or_hpa_scrna()
    hpa_sc = parse_hpa_sc(hpa_sc_zip) if hpa_sc_zip.exists() else pd.DataFrame()
    if not hpa_sc.empty:
        hpa_sc.to_csv(RESULTS / "hpa_single_cell_type_9genes.csv", index=False)
    if long149.empty:
        long149 = synthesize_from_hpa_if_needed(hpa_sc, "GSE149614_HPA_fallback")
        notes.append("GSE149614 used HPA fallback")
    if long125.empty:
        long125 = synthesize_from_hpa_if_needed(hpa_sc, "GSE125449_HPA_fallback")
        notes.append("GSE125449 used HPA fallback")

    # If TISCH meta exists and looks tabular, try to join ?? otherwise keep fallback labelled.
    peaks149 = peak_from_patient_celltype(long149) if not long149.empty else pd.DataFrame()
    peaks125 = peak_from_patient_celltype(long125) if not long125.empty else pd.DataFrame()
    if not long149.empty:
        long149.to_csv(RESULTS / "GSE149614_pseudobulk_9genes.csv", index=False)
        plot_dotplot(long149, "Nine genes by major class (see source_note)", FIGURES / "fig4_GSE149614_dotplot.png")
    if not long125.empty:
        long125.to_csv(RESULTS / "GSE125449_pseudobulk_9genes.csv", index=False)
        plot_dotplot(long125, "Nine genes by major class, second cohort file", FIGURES / "fig4_GSE125449_dotplot.png")
    if not peaks149.empty:
        peaks149.to_csv(RESULTS / "GSE149614_patient_peak_class.csv", index=False)
    if not peaks125.empty:
        peaks125.to_csv(RESULTS / "GSE125449_patient_peak_class.csv", index=False)

    stop = stop_rule(peaks149 if not peaks149.empty else pd.DataFrame({"gene": [], "peak_class_mode": []}), peaks125 if not peaks125.empty else None)
    # If only HPA fallback with a single fake patient, stop-rule must not claim GSE149614/GSE125449 separation.
    used_fallback = (not long149.empty) and ("HPA" in str(long149.get("source_note", pd.Series([""])).iloc[0]))
    if used_fallback:
        stop["stop_rule_triggered"] = True
        stop["continue_celltype_story"] = False
        stop["data_status"] = "GSE processed matrices not in GEO suppl; HPA scRNA used only as resource table. Stop-rule cannot be passed on fallback data."
        notes.append(stop["data_status"])
    (RESULTS / "stop_rule.json").write_text(json.dumps(stop, indent=2), encoding="utf-8")

    # TCGA
    expr_path, pheno_path, surv_path = download_xena_lihc()
    tcga_stats = pd.DataFrame()
    cox = pd.DataFrame()
    try:
        expr = read_xena_expr_9genes(expr_path)
        expr.to_csv(RESULTS / "TCGA_LIHC_9genes_expression.csv")
        tcga_stats = tcga_tumor_normal(expr)
        tcga_stats.to_csv(RESULTS / "TCGA_LIHC_tumor_vs_normal.csv", index=False)
        plot_tcga_box(expr, tcga_stats, FIGURES / "fig5_TCGA_tumor_normal.png")
        cox = cox_os(expr, surv_path)
        cox.to_csv(RESULTS / "supplement_TCGA_Cox_OS.csv", index=False)
        notes.append(f"TCGA genes {list(expr.index)} samples {expr.shape[1]}")
    except Exception as exc:
        notes.append(f"TCGA failed: {exc}")

    # HPA IHC
    ihc_rows = []
    for gene in NINE:
        entries, err = hpa_ihc_liver(gene)
        if err:
            ihc_rows.append({"gene": gene, "error": err})
        else:
            if not entries:
                ihc_rows.append({"gene": gene, "tissue": "liver_search", "level": "see JSON cache", "keys": ""})
            else:
                ihc_rows.extend(entries)
    pd.DataFrame(ihc_rows).to_csv(RESULTS / "HPA_IHC_liver_9genes.csv", index=False)

    (RESULTS / "week23_notes.json").write_text(json.dumps(notes, indent=2), encoding="utf-8")
    print("STOP", stop)
    print("\n".join(notes))


if __name__ == "__main__":
    main()
