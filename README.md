# Artesunate-HCC two-group in silico paper

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23156173.svg)](https://doi.org/10.5281/zenodo.23156173)

Home folder: `C:\Users\20917\Desktop\art&hcc` (do not write outputs into HA_docking).

Cite: https://doi.org/10.5281/zenodo.23156173  
GitHub: https://github.com/zhazhayv/art-hcc-twogroup

Locked lists:

- Network: PARP1, BCL2, CASP3, BAX, FAS, CASP9
- Experimental: GBA, MMP9, OGA (`MGEA5`)

STRING: use `locked_inputs/Full67_STRING_raw.tsv` only (660 edges). Never the 536-edge snapshot.

MMP9 Swiss probability: `locked_inputs/MANIFEST.json` field `MMP9_swiss_probability` only.

## Run

```
py -3 01_lock_inputs.py
py -3 02_threshold_grid.py
py -3 03_scrna_tcga_hpa.py
py -3 04_dock_experimental_pockets.py
py -3 05_build_manuscript.py
```

Target journal: Frontiers in Pharmacology. Backup: Computational and Structural Biotechnology Journal.
No wet experiments.
