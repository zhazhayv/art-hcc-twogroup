# -*- coding: utf-8 -*-
"""Pre-registered Swiss x GeneCards x STRING grid, miss-path, CTD labels, g:Profiler dumbbell."""
from __future__ import annotations

import hashlib
import json
import ssl
import time
import urllib.parse
import urllib.request
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

from constants import (
    EXPERIMENTAL3,
    FIGURES,
    GC_SCORE,
    LOCKED,
    NETWORK6,
    RESULTS,
    STRING_SCORE,
    SWISS_P,
)

SSL_CTX = ssl._create_unverified_context()
UA = {"User-Agent": "art-hcc-twogroup/1.0"}


def http_get(url, timeout=180):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
        return resp.read()


def http_post_json(url, payload, timeout=180):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url, data=data, headers={**UA, "Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
        return json.loads(resp.read().decode())


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def load_locked():
    swiss = pd.read_csv(LOCKED / "swiss_raw.csv")
    ctd = pd.read_csv(LOCKED / "ctd_filtered.csv")
    gc = pd.read_csv(LOCKED / "genecards_hcc.csv")
    manifest = json.loads((LOCKED / "MANIFEST.json").read_text(encoding="utf-8"))
    return swiss, ctd, gc, manifest


def compound_targets(swiss, ctd, pmin: float):
    ctd_genes = {str(x).upper() for x in ctd["Gene Symbol"].dropna()}
    swiss_hi = {
        str(r["Common name"]).upper()
        for _, r in swiss.iterrows()
        if pd.notna(r["Common name"]) and float(r["Probability*"]) >= pmin
    }
    return sorted(ctd_genes | swiss_hi)


def disease_genes(gc, score_min: float):
    sub = gc[gc["Relevance_Score"] > score_min]
    return set(sub["Symbol"].astype(str).str.upper())


def string_network(genes, min_score: float, cache_dir: Path):
    genes = sorted(set(genes))
    key = hashlib.sha1((",".join(genes) + f"|{min_score:.3f}").encode()).hexdigest()[:16]
    cache = cache_dir / f"string_{key}.tsv"
    if cache.exists():
        txt = cache.read_text(encoding="utf-8")
    else:
        ids = "%0d".join(genes)
        # required_score is 0-1000
        req = int(round(min_score * 1000))
        url = (
            "https://string-db.org/api/tsv/network?"
            + urllib.parse.urlencode(
                {"identifiers": "\r".join(genes), "species": 9606, "required_score": req, "caller_identity": "art_hcc"}
            )
        )
        try:
            raw = http_get(url.replace("+", "%20") if False else (
                "https://string-db.org/api/tsv/network?identifiers="
                + "%0d".join(genes)
                + f"&species=9606&required_score={req}&caller_identity=art_hcc"
            ))
            txt = raw.decode("utf-8", errors="replace")
        except Exception as exc:
            txt = f"ERROR\t{exc}\n"
        cache.write_text(txt, encoding="utf-8")
        time.sleep(0.4)
    if txt.startswith("ERROR") or "preferredName_A" not in txt.splitlines()[0] if txt.strip() else True:
        # fallback: local Full67 STRING filtered to genes present
        local = pd.read_csv(LOCKED / "Full67_STRING_raw.tsv", sep="\t")
        a = "preferredName_A" if "preferredName_A" in local.columns else None
        if a:
            sub = local[
                (local["preferredName_A"].isin(genes))
                & (local["preferredName_B"].isin(genes))
                & (local["score"] >= min_score)
            ]
            return sub
        return pd.DataFrame()
    df = pd.read_csv(StringIO(txt), sep="\t")
    if "score" in df.columns:
        df = df[df["score"] >= min_score]
    return df


def degree_top5(edge_df, genes):
    g = nx.Graph()
    g.add_nodes_from(genes)
    if edge_df is None or edge_df.empty:
        return [], {}
    col_a = "preferredName_A" if "preferredName_A" in edge_df.columns else edge_df.columns[2]
    col_b = "preferredName_B" if "preferredName_B" in edge_df.columns else edge_df.columns[3]
    for _, r in edge_df.iterrows():
        a, b = str(r[col_a]), str(r[col_b])
        if a != b:
            g.add_edge(a, b)
    deg = dict(g.degree())
    ranked = sorted(deg.items(), key=lambda x: (-x[1], x[0]))
    return [n for n, _ in ranked[:5]], deg


def gprofiler_kegg(genes, cache_dir: Path):
    genes = sorted(set(genes))
    key = hashlib.sha1(",".join(genes).encode()).hexdigest()[:16]
    cache = cache_dir / f"gprofiler_{key}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": ["KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": True,
    }
    try:
        out = http_post_json("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", payload)
        cache.write_text(json.dumps(out), encoding="utf-8")
        time.sleep(0.4)
        return out
    except Exception as exc:
        rec = {"error": str(exc), "result": []}
        cache.write_text(json.dumps(rec), encoding="utf-8")
        return rec


def first_kegg(gp):
    rows = gp.get("result") or []
    kegg = [r for r in rows if str(r.get("source", "")).startswith("KEGG") or r.get("source") == "KEGG"]
    if not kegg:
        return "", None, None
    kegg.sort(key=lambda r: float(r.get("p_value", 1)))
    top = kegg[0]
    return top.get("name", ""), top.get("native", ""), top.get("p_value")


def kegg_rank_map(gp):
    rows = gp.get("result") or []
    kegg = [r for r in rows if r.get("source") == "KEGG" or str(r.get("native", "")).startswith("KEGG")]
    kegg.sort(key=lambda r: float(r.get("p_value", 1)))
    return {r.get("name"): i + 1 for i, r in enumerate(kegg)}


def miss_path(swiss, ctd, gc, manifest):
    swiss_p = {
        str(r["Common name"]).upper(): float(r["Probability*"])
        for _, r in swiss.iterrows()
        if pd.notna(r["Common name"])
    }
    ctd_genes = {str(x).upper() for x in ctd["Gene Symbol"].dropna()}
    gc_map = dict(zip(gc["Symbol"].astype(str).str.upper(), gc["Relevance_Score"]))
    from constants import EXPERIMENTAL_ALIASES

    rows = []
    for gene in EXPERIMENTAL3:
        names = EXPERIMENTAL_ALIASES[gene]
        in_swiss = any(n in swiss_p for n in names)
        p = next((swiss_p[n] for n in names if n in swiss_p), None)
        in_ctd = any(n in ctd_genes for n in names)
        gc_hits = {n: gc_map.get(n) for n in names if n in gc_map}
        gc_score = max(gc_hits.values()) if gc_hits else None
        gc_used = max(gc_hits, key=gc_hits.get) if gc_hits else None
        steps = []
        if not in_swiss:
            steps.append("absent_from_swiss")
        elif p < 0.1:
            steps.append(f"swiss_probability_below_0.1 (p={p})")
        if not in_ctd:
            steps.append("absent_from_ctd")
        if gc_score is None:
            steps.append("absent_from_genecards")
        elif gc_score <= 2:
            steps.append(f"genecards_score_le_2 ({gc_used}={gc_score})")
        else:
            steps.append(f"survives_genecards_gt2 as {gc_used}={gc_score}")
        if gene == "MMP9":
            if p != manifest["MMP9_swiss_probability"]:
                raise RuntimeError("MMP9 probability drifted from locked MANIFEST")
        rows.append(
            {
                "gene": gene,
                "aliases_checked": ",".join(names),
                "in_swiss": in_swiss,
                "swiss_probability": manifest["MMP9_swiss_probability"] if gene == "MMP9" else p,
                "in_ctd": in_ctd,
                "genecards_symbol_used": gc_used,
                "genecards_relevance": None if gc_score is None or pd.isna(gc_score) else float(gc_score),
                "miss_steps": "; ".join(steps) if steps else "recovered_at_primary_cutoffs",
            }
        )
    return rows


def ctd_labels(ctd):
    rows = []
    for gene in NETWORK6:
        hit = ctd[ctd["Gene Symbol"].astype(str).str.upper() == gene]
        if hit.empty:
            rows.append(
                {
                    "gene": gene,
                    "in_ctd": False,
                    "interaction_class": "not_in_export",
                    "raw_interaction": "",
                }
            )
            continue
        raw = str(hit.iloc[0].get("Interaction", ""))
        rows.append(
            {
                "gene": gene,
                "in_ctd": True,
                "interaction_class": "untyped_gene_chemical_record",
                "raw_interaction": raw,
                "note": "Export is Harmonizome CTD association; binding vs expression vs inferred not labelled.",
            }
        )
    return rows


def plot_grid(df, path: Path):
    fig, ax = plt.subplots(figsize=(10.2, 4.8), dpi=300)
    labels = [
        f"P>={r.swiss_p:g}\nGC>{r.gc_score:g}\nSTR {r.string_score:.3f}"
        for r in df.itertuples()
    ]
    ax.bar(range(len(df)), df["n_intersection"], color="#4C72B0")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Intersection gene count")
    ax.set_title("Pre-registered threshold grid")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_dumbbell(rank_a, rank_b, path: Path):
    names = sorted(set(rank_a) | set(rank_b), key=lambda n: (rank_a.get(n, 99) + rank_b.get(n, 99), n))
    names = names[:15]
    fig, ax = plt.subplots(figsize=(7.2, 5.5), dpi=300)
    y = list(range(len(names)))[::-1]
    xa = [rank_a.get(n, None) for n in names]
    xb = [rank_b.get(n, None) for n in names]
    for yi, a, b in zip(y, xa, xb):
        if a and b:
            ax.plot([a, b], [yi, yi], color="#888888", lw=1.2, zorder=1)
        if a:
            ax.scatter([a], [yi], color="#4C72B0", s=36, zorder=2, label="Strict_14" if yi == y[0] else "")
        if b:
            ax.scatter([b], [yi], color="#C44E52", s=36, zorder=2, label="Full_67" if yi == y[0] else "")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("KEGG rank (g:Profiler g:SCS, same run engine)")
    ax.invert_xaxis()
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("KEGG rank shift after evidence filtering")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    swiss, ctd, gc, manifest = load_locked()
    cache = RESULTS / "api_cache"
    cache.mkdir(exist_ok=True)
    grid_rows = []
    gene_sets = {}
    for pmin in SWISS_P:
        art = compound_targets(swiss, ctd, pmin)
        for gmin in GC_SCORE:
            disease = disease_genes(gc, gmin)
            inter = sorted(set(art) & disease)
            key_genes = tuple(inter)
            for smin in STRING_SCORE:
                edges = string_network(inter, smin, cache)
                top5, deg = degree_top5(edges, inter)
                gp = gprofiler_kegg(inter, cache)
                kegg_name, kegg_id, kegg_p = first_kegg(gp)
                alias_flat = {"GBA": ["GBA", "GBA1"], "MMP9": ["MMP9"], "MGEA5": ["MGEA5", "OGA"]}
                present_exp = [g for g, ns in alias_flat.items() if any(n in inter for n in ns)]
                present_net = [g for g in NETWORK6 if g in inter]
                grid_rows.append(
                    {
                        "swiss_p": pmin,
                        "gc_score": gmin,
                        "string_score": smin,
                        "n_compound": len(art),
                        "n_disease": len(disease),
                        "n_intersection": len(inter),
                        "degree_top5": ",".join(top5),
                        "jaccard_intersection_vs_network6": jaccard(inter, NETWORK6),
                        "jaccard_degree_top5_vs_network6": jaccard(top5, NETWORK6),
                        "jaccard_intersection_vs_experimental3": jaccard(inter, EXPERIMENTAL3),
                        "network6_recovered": ",".join(present_net),
                        "experimental3_in_set": ",".join(present_exp) if present_exp else "",
                        "GBA_in_set": any(n in inter for n in ["GBA", "GBA1"]),
                        "MMP9_in_set": "MMP9" in inter,
                        "MGEA5_in_set": any(n in inter for n in ["MGEA5", "OGA"]),
                        "gprofiler_first_kegg": kegg_name,
                        "gprofiler_first_kegg_id": kegg_id,
                        "gprofiler_first_kegg_p": kegg_p,
                    }
                )
                gene_sets[f"P{pmin}_GC{gmin}_STR{smin}"] = inter

    grid = pd.DataFrame(grid_rows)
    grid.to_csv(RESULTS / "threshold_grid.csv", index=False)
    (RESULTS / "grid_gene_sets.json").write_text(json.dumps(gene_sets, indent=2), encoding="utf-8")

    miss = miss_path(swiss, ctd, gc, manifest)
    pd.DataFrame(miss).to_csv(RESULTS / "experimental_miss_path.csv", index=False)
    labels = ctd_labels(ctd)
    pd.DataFrame(labels).to_csv(RESULTS / "network6_ctd_labels.csv", index=False)

    from constants import FULL67, STRICT14

    gp_s = gprofiler_kegg(STRICT14, cache)
    gp_f = gprofiler_kegg(FULL67, cache)
    (RESULTS / "gprofiler_strict14.json").write_text(json.dumps(gp_s), encoding="utf-8")
    (RESULTS / "gprofiler_full67.json").write_text(json.dumps(gp_f), encoding="utf-8")
    plot_dumbbell(kegg_rank_map(gp_s), kegg_rank_map(gp_f), FIGURES / "fig3_kegg_dumbbell.png")
    plot_grid(grid, FIGURES / "fig2_threshold_grid.png")

    # Jaccard heatmap-style table figure for experimental miss
    fig, ax = plt.subplots(figsize=(8.4, 2.6), dpi=300)
    ax.axis("off")
    cell = [[r["gene"], r["swiss_probability"], r["in_ctd"], r["genecards_relevance"], r["miss_steps"]] for r in miss]
    table = ax.table(
        cellText=[[str(c) if c is not None else "" for c in row] for row in cell],
        colLabels=["Gene", "Swiss P (locked)", "In CTD", "GeneCards score", "Miss path"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.auto_set_column_width(list(range(5)))
    ax.set_title("Experimental-group miss path (GBA / MMP9 / OGA=MGEA5)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_miss_path.png", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "n_grid_cells": len(grid),
        "swiss_thresholds_collapse": bool(
            grid.groupby(["gc_score", "string_score"])["n_intersection"].nunique().max() == 1
        ),
        "MMP9_swiss_probability": manifest["MMP9_swiss_probability"],
        "primary_strict14_first_kegg": first_kegg(gp_s)[0],
        "primary_full67_first_kegg": first_kegg(gp_f)[0],
        "ctd_typing": "untyped_gene_chemical_record",
    }
    (RESULTS / "week1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(grid.to_string(index=False))
    print("summary", summary)


if __name__ == "__main__":
    main()
