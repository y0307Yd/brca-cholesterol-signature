# Cholesterol metabolism in breast cancer: an integrated multi-omics study

This repository provides the reproducible analysis pipeline for the manuscript
"Cholesterol Metabolism Gene Signatures Distinguish Estrogen-Receptor Status
and Define Immune-Relevant Molecular Subtypes of Breast Cancer".

## Study design

TCGA-BRCA (discovery) and METABRIC (external) were combined with WGCNA,
machine-learning feature selection, consensus clustering, single-cell
localisation, SMR/HEIDI Mendelian randomisation, Bayesian colocalisation
(coloc.abf and SuSiE/coloc.susie), methylation-expression correlation and
three independent Affymetrix validation cohorts (GSE21653, GSE7390, METABRIC).

## Quick start

```bash
python -m pip install -r requirements.txt
# download public data (see data/download_data.py for URLs and checksums)
python data/download_data.py --dir ../data
# run the pipeline in order (see README section "Pipeline order")
```

R scripts require R >= 4.6 with `susieR`, `coloc`, `BiocManager` and
`hgu133a.db` (used for the GSE7390 probe mapping).

## Pipeline order

1. `code/bc_01_prep_data.py` -> `bc_02_wgcna.py` -> `bc_03_ml_signature.py` ->
   `bc_04_subtypes.py` -> `bc_05_figures.py` -> `bc_06_summary.py` ->
   `bc_07_enrich_fc.py` (discovery analysis)
2. `finish_pam50.py`, `finish_pam50_expression.py`, `finish_ssgsea.py`
   (subtype characterisation)
3. `finish_geo_gse21653.py`, `finish_geo_gse7390.py` (external validation)
4. `finish_scrna_gse161529.py`, `finish_scrna_lr_deep.py`,
   `scRNA_cellxgene.py` (single-cell)
5. `make_gwas_ma.py`, `make_tnbc_ma.py`, `build_onco_ma.py`,
   `smr_breast_template.sh`, `make_supp_table_s10.py` (SMR/HEIDI)
6. `coloc_fdps.py`, `coloc_sens.py`, `fdps_conditional.py`,
   `fdps_fstat_credset.py`, `prep_susie_fdps.py`, `susie_coloc_fdps.R`
   (colocalisation)

## Key results

- 11-gene cholesterol-metabolism ER classifier: TCGA CV AUC 0.946, METABRIC
  AUC 0.922, GSE21653 AUC 0.847, GSE7390 AUC 0.903.
- FDPS is a putatively causal protective gene for breast cancer risk
  (SMR P = 9.8e-7; ER+ P = 1.5e-5 in the BCAC OncoArray sensitivity analysis).
- Single-causal-variant coloc.abf gives PP.H4 = 0.9985; SuSiE-based
  multi-signal colocalisation is less conclusive (max PP.H4 = 0.44).

Figure files in `results/figures/` are named `figNN_*.png` (and 300-dpi
`figNN_*.tiff`), where NN is the figure number in the manuscript (Figures
1-16; `Supplementary_*.png` are Supplementary Figures S1-S2).

## Citation

If you use this pipeline, please cite the manuscript (DOI to be added) and the
primary data sources listed in `DATA_LICENSES.md`.
