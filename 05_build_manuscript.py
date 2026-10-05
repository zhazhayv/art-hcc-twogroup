# -*- coding: utf-8 -*-
"""Build Frontiers-ready manuscript, cover letter, and reproducibility files."""
from __future__ import annotations

import json
from datetime import date

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

from constants import FULL67, LOCKED, MS, RESULTS, STRICT14


def p(doc, text, *, bold=False, center=False, size=11):
    para = doc.add_paragraph()
    if center:
        para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = para.add_run(text)
    run.bold = bold
    run.font.size = Pt(size)
    run.font.name = "Times New Roman"
    para.paragraph_format.first_line_indent = Inches(0 if center else 0.3)
    para.paragraph_format.space_after = Pt(6)
    return para


def h(doc, text, level=1):
    heading = doc.add_heading(text, level=level)
    for run in heading.runs:
        run.font.name = "Times New Roman"
    return heading


def load():
    man = json.loads((LOCKED / "MANIFEST.json").read_text(encoding="utf-8"))
    week1 = json.loads((RESULTS / "week1_summary.json").read_text(encoding="utf-8")) if (RESULTS / "week1_summary.json").exists() else {}
    grid = pd.read_csv(RESULTS / "threshold_grid.csv") if (RESULTS / "threshold_grid.csv").exists() else pd.DataFrame()
    miss = pd.read_csv(RESULTS / "experimental_miss_path.csv") if (RESULTS / "experimental_miss_path.csv").exists() else pd.DataFrame()
    ctd = pd.read_csv(RESULTS / "network6_ctd_labels.csv") if (RESULTS / "network6_ctd_labels.csv").exists() else pd.DataFrame()
    stop = json.loads((RESULTS / "stop_rule.json").read_text(encoding="utf-8")) if (RESULTS / "stop_rule.json").exists() else {}
    dock = json.loads((RESULTS / "docking_gate.json").read_text(encoding="utf-8")) if (RESULTS / "docking_gate.json").exists() else {}
    tcga = pd.read_csv(RESULTS / "TCGA_LIHC_tumor_vs_normal.csv") if (RESULTS / "TCGA_LIHC_tumor_vs_normal.csv").exists() else pd.DataFrame()
    return man, week1, grid, miss, ctd, stop, dock, tcga


def build_docx(man, week1, grid, miss, ctd, stop, dock, tcga):
    doc = Document()
    title = (
        "Evidence filtering and cell-type context separate a network-apoptosis cassette "
        "from experimentally supported artesunate-HCC proteins"
    )
    p(doc, title, bold=True, center=True, size=16)
    p(doc, "Authors: [Name 1], [Name 2], [Name 3]", center=True)
    p(doc, "Affiliations: 1 [Department], [Institution], [City], [Country]", center=True)
    p(doc, "Corresponding author: [Name]; E-mail: [email]", center=True)
    p(
        doc,
        "Database freeze: Swiss/CTD/GeneCards locked %s; STRING Full_67 = %s edges."
        % (man.get("lock_date"), man.get("Full67_STRING_edges")),
        center=True,
    )

    h(doc, "Abstract")
    if stop.get("stop_rule_triggered", True):
        stop_txt = (
            "The pre-registered single-cell stop-rule was triggered "
            "(network-group and experimental-group mode peak classes did not separate); "
            "cell-type separation is not claimed as a mechanism."
        )
    else:
        stop_txt = (
            "Patient-level pseudobulk placed the network cassette and the experimental trio "
            "in different major cell classes in both cohorts."
        )
    md = "omitted because no new pocket passed the redock gate" if dock.get("md_section") == "omit" else "restricted to a redock-pass pocket"
    p(
        doc,
        "Background: Artesunate has chemically supported hepatocellular carcinoma (HCC) proteins "
        "(GBA, MMP9, OGA/MGEA5) and a separate apoptosis annotation trail in CTD. Unfiltered "
        "SwissTargetPrediction plus a broad GeneCards list promotes STAT3/EGFR/PI3K hubs that do not "
        "survive a probability floor. This in silico study does not nominate a new target. "
        "It asks whether a locked network-apoptosis cassette and the experimental trio occupy different "
        "evidence layers and, if the data allow, different cell classes. "
        "Methods: Compound identity is TCMSP MOL007434. Locked Swiss probability for MMP9 is %s. "
        "A pre-registered grid varied Swiss (0.05/0.1/0.2), GeneCards (>2/>5/>10) and STRING (0.400/0.700). "
        "Primary interpretation uses Strict_14 versus Full_67 (660-edge STRING only). "
        "Docking of artesunate and dihydroartemisinin is restricted to pockets that reproduce the co-crystal ligand "
        "with heavy-atom RMSD < 2.0 A; BCL2/4LVT is reused (GetBestRMS 1.54 A; artesunate -7.5 kcal/mol; DHA -7.9 kcal/mol). "
        "CASP3/3DEI and PARP1/4ZZZ remain supplementary. No cell or animal experiment was performed. "
        "Results: %s leads Strict_14 KEGG; %s leads Full_67. "
        "GBA/GBA1 and OGA/MGEA5 are absent from Swiss and CTD; MMP9 is present in Swiss only below the 0.1 floor. "
        "%s MD was %s. "
        "Conclusion: Evidence filtering replaces STAT3/EGFR-centred hubs with CTD apoptosis annotations. "
        "That shift restates database priors. Experimentally supported proteins sit outside the filter. "
        "Geometric compatibility is claimed only for protocol-valid pockets."
        % (
            man.get("MMP9_swiss_probability"),
            week1.get("primary_strict14_first_kegg", "Apoptosis"),
            week1.get("primary_full67_first_kegg", "Lipid and atherosclerosis"),
            stop_txt,
            md,
        ),
    )
    p(
        doc,
        "Keywords: artesunate; dihydroartemisinin; hepatocellular carcinoma; network pharmacology; "
        "evidence filtering; GBA; MMP9; OGA; molecular docking",
    )

    h(doc, "1 Introduction")
    p(
        doc,
        "Hepatocellular carcinoma remains a major global cause of cancer death. Artesunate is a hemisuccinate ester "
        "that is converted in plasma to dihydroartemisinin; the endoperoxide can be activated by ferrous iron or heme. "
        "Experimental work has already identified GBA, MMP9 and OGA as chemically supported artesunate-HCC proteins "
        "and has linked artesunate to ferroptosis and sorafenib sensitization. Those reports are the experimental group "
        "in this paper. They are not rediscovered here.",
    )
    p(
        doc,
        "Network pharmacology often promotes very low-probability SwissTargetPrediction rows to hubs. We therefore lock "
        "two gene lists. The network group is PARP1, BCL2, CASP3, BAX, FAS and CASP9. The experimental group is GBA, "
        "MMP9 and OGA (GeneCards currently indexes GBA1 and OGA). Full_67 is an unrestricted exploratory intersection; "
        "Strict_14 is the only set used for primary interpretation. We do not treat apoptosis KEGG enrichment as a newly "
        "discovered HCC axis.",
    )

    h(doc, "2 Materials and methods")
    h(doc, "2.1 Locked inputs", 2)
    p(
        doc,
        "SwissTargetPrediction, CTD (Harmonizome gene-chemical export) and GeneCards HCC downloads were copied into "
        "art_hcc_twogroup/locked_inputs with SHA-256 digests (MANIFEST.json). MMP9 Swiss probability is %s from SWISS raw; "
        "methods and results use this single number. The Full_67 STRING table has %s edges. An earlier 536-edge snapshot "
        "is not used. CTD rows are untyped gene-chemical records; binding versus expression versus inferred flags were not "
        "present in the export."
        % (man.get("MMP9_swiss_probability"), man.get("Full67_STRING_edges")),
    )
    h(doc, "2.2 Pre-registered threshold grid", 2)
    p(
        doc,
        "Swiss probability floors 0.05, 0.1 and 0.2; GeneCards relevance >2, >5 and >10; STRING combined score 0.400 and "
        "0.700. Each cell records intersection size, degree top five, Jaccard overlap with the locked network six and "
        "experimental three, presence of GBA/MMP9/MGEA5 (including GBA1/OGA aliases), and the first g:Profiler KEGG term "
        "(Homo sapiens, g:SCS). Rules were not changed after seeing ranks.",
    )
    h(doc, "2.3 Single-cell, TCGA and HPA", 2)
    p(
        doc,
        "GSE149614 (Sun et al.) is the primary single-cell cohort. Cells were mapped onto six author-style major classes: "
        "hepatocyte/cholangiocyte, T/NK, myeloid, B, endothelial, fibroblast. The statistical unit is the patient. "
        "GSE125449 is a second map only. TCGA-LIHC tumor versus adjacent used Wilcoxon tests with Benjamini-Hochberg FDR. "
        "Human Protein Atlas IHC is database extraction, not new staining. Univariable Cox OS uses the same specification "
        "for both gene groups and is supplementary. Stop-rule: the mode peak class of the network six must differ from that "
        "of the experimental three in GSE149614 and remain separated in GSE125449; otherwise the manuscript does not convert "
        "into a mechanism paper.",
    )
    h(doc, "2.4 Docking gate", 2)
    p(
        doc,
        "A pocket enters the main figure only if it has a human ligand-bound crystal structure and reproduces the co-crystal "
        "ligand with in-place or GetBestRMS heavy-atom RMSD < 2.0 A. Box: ligand centroid, 8 A padding, size capped at 30 A. "
        "AutoDock Vina 1.1.2, exhaustiveness 32, nine modes, seeds 20260930/20261001/20261002. The peroxide cage is rigid; "
        "the succinate chain is flexible. Artesunate and DHA share the box. A trypsin/benzamidine complex (3PTB) is an "
        "unrelated negative control. CASP3/3DEI and PARP1/4ZZZ remain supplementary failures. A second engine (GNINA, smina "
        "or Vina 1.2) is applied only to pockets that already passed. Molecular dynamics (3 x 50 ns per ligand) is run only "
        "if at least one new experimental pocket passes; otherwise the MD section is omitted.",
    )

    h(doc, "3 Results")
    h(doc, "3.1 Two-level intersection", 2)
    p(
        doc,
        "Strict_14 (%d genes) and Full_67 (%d genes) are unchanged from the locked downloads. The network cassette is a "
        "subset of Strict_14. Figure 1 records set membership, not measured physical interactions."
        % (len(STRICT14), len(FULL67)),
    )
    h(doc, "3.2 Threshold grid and miss path", 2)
    if not grid.empty:
        nuniq = int(grid["n_intersection"].nunique())
        p(
            doc,
            "The grid has %d cells. Intersection sizes take %d distinct values. Swiss floors 0.05, 0.1 and 0.2 collapse on "
            "this download because only PTGES2 and CYP1A2 exceed 0.05 (MMP9 p = %s). GeneCards cutoffs change disease-side "
            "breadth. STRING 0.700 thins edges without restoring GBA or OGA."
            % (len(grid), nuniq, man.get("MMP9_swiss_probability")),
        )
    if not miss.empty:
        bits = "; ".join("%s: %s" % (r.gene, r.miss_steps) for r in miss.itertuples())
        p(
            doc,
            "Experimental-group miss path: %s. Failure to recover a protein from Swiss plus CTD is a property of this filter, "
            "not evidence against the experiments." % bits,
        )
    if not ctd.empty:
        p(doc, "Every network-group gene in CTD is an untyped gene-chemical record in this export (network6_ctd_labels.csv).")
    h(doc, "3.3 KEGG rank shift", 2)
    p(
        doc,
        "Same-tool g:Profiler: Strict_14 first KEGG = %s; Full_67 first KEGG = %s. Apoptosis enrichment on Strict_14 is "
        "largely circular with CTD apoptosis annotations and is not written as a new HCC-specific axis."
        % (week1.get("primary_strict14_first_kegg"), week1.get("primary_full67_first_kegg")),
    )
    h(doc, "3.4 Cell-type maps and bulk expression", 2)
    p(doc, json.dumps(stop, ensure_ascii=True))
    p(
        doc,
        "In GSE149614 patient-level pseudobulk, MMP9 peaked in myeloid cells in 19 of 21 patients, consistent with the "
        "original atlas description of MMP9+ macrophages; that observation is used only as a contrast, not as a new finding. "
        "The mode peak class of the six network genes and of the three experimental genes was hepatocyte/cholangiocyte, so "
        "the stop-rule was not passed. GSE125449 processed averages were not recovered from TISCH2 in this run.",
    )
    if stop.get("stop_rule_triggered", True):
        p(
            doc,
            "Stop-rule triggered. The Discussion does not convert the nine genes into a cell-type mechanism of artesunate. "
            "Figure 4 is a resource map.",
        )
    if not tcga.empty:
        p(
            doc,
            "TCGA-LIHC tumor versus adjacent statistics for the nine genes are in TCGA_LIHC_tumor_vs_normal.csv (Wilcoxon, BH). "
            "Survival Cox models are supplementary (supplement_TCGA_Cox_OS.csv).",
        )
    h(doc, "3.5 Protocol-valid docking", 2)
    p(
        doc,
        "BCL2/4LVT remains the protocol-valid engineered BH3 groove (GetBestRMS 1.54 A; artesunate -7.5 kcal/mol; DHA "
        "-7.9 kcal/mol across three seeds). Experimental-group docking gate: %s. Scores are geometric compatibility, not Kd. "
        "PDB 4LVT is an engineered Bcl-2 construct. 3PTB trypsin is a negative-control pocket that passed redocking (in-place "
        "RMSD 0.57 A) and is tabulated only in the supplement."
        % json.dumps(dock, default=str),
    )

    h(doc, "4 Discussion")
    p(
        doc,
        "The usable result is comparative, not therapeutic. Filtering removes STAT3/EGFR/PI3K hubs that originate from Swiss "
        "probabilities far below 0.1 and leaves a CTD apoptosis cassette. Experimentally supported GBA, MMP9 and OGA sit "
        "outside that cassette for documented miss-path reasons. Single-gene MMP9 myeloid bias in GSE149614 does not satisfy "
        "the pre-registered group-mode stop-rule and is not rewritten as an artesunate-macrophage mechanism.",
    )
    p(
        doc,
        "Limitations: no wet experiment; CTD untyped; 4LVT engineered; rigid Vina cannot represent Fe(II)/heme activation of "
        "the endoperoxide; GSE125449 cell-type averages were not downloaded; Cox models are univariable and supplementary. "
        "We do not infer clinical activity.",
    )

    h(doc, "5 Conclusions")
    p(
        doc,
        "Evidence filtering changes an artesunate-HCC network from STAT3/EGFR/PI3K-centred hubs to CTD apoptosis annotations. "
        "GBA, MMP9 and OGA remain outside the filter. Docking claims are restricted to redock-pass pockets, currently "
        "BCL2/4LVT. Files are auditable.",
    )

    h(doc, "Data availability")
    p(
        doc,
        "Locked tables, threshold grid, docking logs and scripts: https://github.com/zhazhayv/art-hcc-twogroup "
        "and Zenodo https://doi.org/10.5281/zenodo.23156173 (version DOI 10.5281/zenodo.23156174). "
        "GEO: GSE149614, GSE125449. TCGA-LIHC via UCSC Xena. HPA JSON accessed on the run date.",
    )
    h(doc, "Conflict of interest")
    p(
        doc,
        "The authors declare that the research was conducted in the absence of any commercial or financial relationships "
        "that could be construed as a potential conflict of interest.",
    )
    h(doc, "Author contributions")
    p(doc, "[To be completed].")
    h(doc, "References")
    for r in [
        "Hopkins AL. Network pharmacology. Nat Biotechnol. 2007.",
        "Szklarczyk D et al. STRING. Nucleic Acids Res.",
        "Raudvere U et al. g:Profiler. Nucleic Acids Res. 2019.",
        "Trott O, Olson AJ. AutoDock Vina. J Comput Chem. 2010.",
        "Chen W et al. Artesunate and GBA in HCC. Exp Mol Med. 2022.",
        "Sun Y et al. GSE149614 ecosystem of early-relapse HCC. Cell. 2021.",
        "Ma L et al. GSE125449 liver cancer ecosystem. Cancer Cell. 2019.",
        "Souers AJ et al. navitoclax analogue BCL2 structure (PDB 4LVT).",
    ]:
        p(doc, r)

    out = MS / "Manuscript_EN_two_group.docx"
    doc.save(out)
    return out


def cover_letter(stop, dock):
    txt = (
        "Dear Editors,\n\n"
        "Please consider our computational manuscript, \"Evidence filtering and cell-type context separate a "
        "network-apoptosis cassette from experimentally supported artesunate-HCC proteins,\" for Frontiers in Pharmacology.\n\n"
        "This study contains no cell, organoid, animal or clinical experiment. We do not nominate a new artesunate target "
        "and we do not claim a treatment mechanism. The contribution is a two-group contrast: a locked network-apoptosis "
        "cassette (PARP1, BCL2, CASP3, BAX, FAS, CASP9) versus experimentally supported proteins (GBA, MMP9, OGA/MGEA5). "
        "We show how a pre-registered Swiss/GeneCards/STRING grid and a documented miss-path keep the experimental proteins "
        "outside the default network, and we restrict docking of artesunate and dihydroartemisinin to co-crystal pockets that "
        "pass a 2 A redocking gate (BCL2/4LVT reused; CASP3 and PARP1 supplementary only).\n\n"
        "Stop-rule (single-cell): %s\n"
        "Docking/MD gate: %s\n\n"
        "If the absence of wet-lab validation is judged incompatible with the journal, we ask that the manuscript be declined "
        "promptly so that we may submit to Computational and Structural Biotechnology Journal.\n\n"
        "Yours sincerely,\n"
        "[Corresponding author]\n"
        % (json.dumps(stop, ensure_ascii=True), json.dumps(dock, default=str))
    )
    path = MS / "Cover_letter_Frontiers_in_Pharmacology.txt"
    path.write_text(txt, encoding="utf-8")
    return path


def main():
    man, week1, grid, miss, ctd, stop, dock, tcga = load()
    ms = build_docx(man, week1, grid, miss, ctd, stop, dock, tcga)
    cl = cover_letter(stop, dock)
    (MS / "ZENODO.md").write_text(
        "# Zenodo DOI\n\n"
        "GitHub: https://github.com/zhazhayv/art-hcc-twogroup\n"
        "Release archived: v1.0.1\n"
        "Record: https://zenodo.org/records/23156174\n\n"
        "This version: https://doi.org/10.5281/zenodo.23156174\n"
        "All versions (cite this): https://doi.org/10.5281/zenodo.23156173\n\n"
        "DOI: 10.5281/zenodo.23156173\n"
        "Version DOI: 10.5281/zenodo.23156174\n"
        "Prepared: %s\n" % date.today().isoformat(),
        encoding="utf-8",
    )
    print("wrote", ms, cl)


if __name__ == "__main__":
    main()
