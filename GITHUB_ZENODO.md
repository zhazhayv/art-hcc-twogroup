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

Then enable Zenodo GitHub integration, create a Release, and paste the DOI into `manuscript/ZENODO.md`.

Do not upload `results/cohort_cache/` (GEO/Xena raw cache). Derived tables in `results/*.csv` are enough to audit the figures.
