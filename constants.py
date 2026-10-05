# -*- coding: utf-8 -*-
"""Locked gene lists and pre-registered thresholds for the two-group paper."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
# art&hcc/02_two_group_pipeline/art_hcc_twogroup -> art&hcc
ROOT = HERE.parent.parent
EXTRACTS = ROOT / "03_from_HA_docking"
SOURCES = ROOT / "01_network_pharmacology_source"
HA = ROOT
LOCKED = HERE / "locked_inputs"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
MS = HERE / "manuscript"
DOCK = HERE / "docking"
for p in (LOCKED, RESULTS, FIGURES, MS, DOCK):
    p.mkdir(parents=True, exist_ok=True)

NETWORK6 = ["PARP1", "BCL2", "CASP3", "BAX", "FAS", "CASP9"]
EXPERIMENTAL3 = ["GBA", "MMP9", "MGEA5"]
EXPERIMENTAL_ALIASES = {
    "GBA": ["GBA", "GBA1"],
    "MMP9": ["MMP9"],
    "MGEA5": ["MGEA5", "OGA"],
    "OGA": ["OGA", "MGEA5"],
}
NINE = NETWORK6 + EXPERIMENTAL3
NINE_EXPR = NETWORK6 + ["GBA", "GBA1", "MMP9", "MGEA5", "OGA"]

SWISS_P = (0.05, 0.1, 0.2)
GC_SCORE = (2.0, 5.0, 10.0)
STRING_SCORE = (0.400, 0.700)

SEEDS = (20260930, 20261001, 20261002)
VINA = ROOT / "bin" / "vina.exe"
P0_BCL2 = EXTRACTS / "_p0_tmp" / "BCL2_4LVT"
ARTHCC_NEWDATA = EXTRACTS / "_arthcc_newdata"
STRICT14 = [
    "AIFM1", "BAX", "BCL2", "CASP3", "CASP9", "CAT", "CYP1A2",
    "CYP2A6", "FAS", "GSTP1", "PARP1", "SP1", "TPT1", "YY1",
]
FULL67 = [
    "ABCB1", "AIFM1", "AKR1B1", "AKR1C3", "ATP2A2", "BAX", "BCL2", "BRD4",
    "CA2", "CA9", "CASP1", "CASP3", "CASP7", "CASP9", "CAT", "CDC25B",
    "CSF1R", "CTSG", "CYP1A2", "CYP2A6", "DNM1L", "EGFR", "EIF4EBP1", "ELANE",
    "EP300", "FAS", "FGFR1", "GSTP1", "GUSB", "HDAC1", "HMGCR", "HSP90AA1",
    "IDH1", "IDO1", "JAK1", "JAK2", "KDR", "MME", "MMEL1", "MMP2", "MMP9",
    "NR1H4", "PARP1", "PFKFB3", "PIK3C2B", "PIK3CA", "PIK3CB", "PIK3CD",
    "PIK3CG", "PIN1", "PPARA", "PPARG", "PRKCA", "PRKCD", "PTGER4", "PTPN1",
    "PTPN11", "SERPINE1", "SLC22A1", "SP1", "STAT3", "TERT", "TLR2", "TLR4",
    "TLR9", "TPT1", "YY1",
]
GSE149614_CLASSES = [
    "Hepatocyte_cholangiocyte", "T_NK", "Myeloid", "B", "Endothelial", "Fibroblast",
]
NEGATIVE_CONTROL = ("3PTB", "trypsin", "BEN")
