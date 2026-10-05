# GitHub and Zenodo (this machine has no git/gh in PATH)

This folder is the paper repository. On a machine with GitHub CLI:

```
cd art_hcc_twogroup
git init
git add README.md LICENSE .zenodo.json .gitignore constants.py 0*.py plot_fig1.py
git add locked_inputs results figures manuscript docking/md
git add -u
git commit -m "Locked two-group artesunate-HCC in silico package."
gh repo create art-hcc-twogroup --private --source . --remote origin --push
```

Zenodo DOI (all versions): https://doi.org/10.5281/zenodo.23156173  
Version v1.0.1: https://doi.org/10.5281/zenodo.23156174  
GitHub: https://github.com/zhazhayv/art-hcc-twogroup

The repository must be **public** for Zenodo GitHub integration (private repos are not listed).

Do not upload `results/cohort_cache/` (GEO/Xena raw cache). Derived tables in `results/*.csv` are enough to audit the figures.
