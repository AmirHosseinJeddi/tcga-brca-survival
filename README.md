# TCGA-BRCA Survival Prediction

Predicting overall survival in breast cancer from gene expression (TCGA-BRCA) using Cox models and Random Survival Forests, benchmarked against a clinical-only baseline.

## Status
- [x] Day 1: data cleaning, overall Kaplan-Meier
- [ ] Day 2: clinical baseline
- [ ] Day 3: feature selection
- [ ] Day 4: Coxnet + Random Survival Forest
- [ ] Day 5: evaluation and risk groups
- [ ] Day 6: write-up

## Data
Downloaded from UCSC Xena (GDC TCGA Breast Cancer hub). Not included in the repo. Put these in `data/`:
- `TCGA-BRCA.star_tpm.tsv.gz`
- `TCGA-BRCA.survival.tsv`

## Setup
```bash
conda env create -f environment.yml
conda activate tcga
jupyter lab
```
