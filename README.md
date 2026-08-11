# Cholesterol metabolism in breast cancer: an integrated multi-omics study

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21873168.svg)](https://doi.org/10.5281/zenodo.21873168)

This repository provides the reproducible analysis pipeline for the manuscript
"Cholesterol Metabolism Gene Signatures Distinguish Estrogen-Receptor Status
and Define Immune-Relevant Molecular Subtypes of Breast Cancer".

## Study design

TCGA-BRCA (discovery) and METABRIC (external) were combined with WGCNA,
machine-learning feature selection, consensus clustering, immune
deconvolution (CIBERSORT/ESTIMATE/ssGSEA), single-cell localisation,
SMR/HEIDI Mendelian randomisation, Bayesian colocalisation (coloc.abf,
SuSiE/coloc.susie, and LIPA/FAXDC2/SREBF1 loci), methylation-expression
correlation, HPA protein evidence, pan-cancer expression and three independent
Affymetrix validation cohorts (GSE21653, GSE7390, GSE20711) plus METABRIC.

## Quick start

```bash
python -m pip install -r requirements.txt
# download public data (see data/download_data.py for URLs and checksums)
python data/download_data.py --dir ../data
# run the pipeline in order (see README section "Pipeline order")
```

R scripts require R >= 4.6 with `susieR`, `coloc`, `e1071` (CIBERSORT),
`BiocManager` and `hgu133a.db` (used for the GSE7390 probe mapping).

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
7. `bonus_analyses.py` (immune-checkpoint genes by subtype; LASSO bootstrap
   stability), `bonus_cibersort.R` (CIBERSORT nu-SVR with LM22),
   `subtype_geo_validation.py` (nearest-centroid subtype mapping in
   GSE21653/GSE7390), `coloc_bonus.py` (LIPA/FAXDC2/SREBF1 colocalisation)
8. `make_lipid_sig_comparison.py`, `make_tripod_checklist.py`,
   `make_flowchart_v20.py`, `make_graphical_abstract.py` (reporting and
   submission artefacts)
9. `finish_geo_gse20711.py` (fourth validation cohort, n = 88)

## Key results

- 11-gene cholesterol-metabolism ER classifier: TCGA CV AUC 0.946, METABRIC
  AUC 0.922, GSE21653 AUC 0.847, GSE7390 AUC 0.903, GSE20711 AUC 0.833.
- LASSO selection stability: all 11 hub genes selected in >= 73.7% of 300
  bootstrap resamples (mean 90.2%).
- Subtype ER gradient reproduced by nearest-centroid mapping in GSE21653
  (P = 4.2e-16) and GSE7390 (P = 4.7e-16); survival differences significant
  in GSE21653 (log-rank P = 0.026), not in GSE7390 (P = 0.48).
- Immune deconvolution: ESTIMATE ImmuneScore highest in C1 (P = 1.9e-3);
  CIBERSORT C1 enrichment of M1/M0 macrophages, activated DC, Tfh, monocytes
  and activated NK; checkpoint genes CD274/PDCD1/CTLA4/LAG3/IDO1 highest in
  C1 (P = 4.9e-5 to 1.7e-20).
- FDPS is a putatively causal protective gene for breast cancer risk
  (SMR P = 9.8e-7; ER+ P = 1.5e-5 in the BCAC OncoArray sensitivity analysis).
- Single-causal-variant coloc.abf gives PP.H4 = 0.9985; SuSiE-based
  multi-signal colocalisation is less conclusive (max PP.H4 = 0.44).
- LIPA/FAXDC2/SREBF1 colocalisation is eQTL-driven without strong shared
  variants (PP.H4 = 0.035/0.043/0.171), reported as a sensitivity negative.
- HPA confirms cytoplasmic FDPS protein with tissue-enhanced breast
  expression; FDPS mRNA detected in all 17 queried TCGA cancer types.

Manuscript v25 (docx + pdf) is the current content-complete version: it
restores the full methodological detail of the original draft, retains all
later additions (immune deconvolution, GSE20711, LIPA/FAXDC2/SREBF1
colocalisation, HPA/pan-cancer evidence, comparison table and TRIPOD-AI
checklist), and adds round-2/3 evidence: METABRIC immune-checkpoint validation
(Table S18), GSE20711 subtype mapping and figure (Figure S4), LASSO
selection-frequency and CIBERSORT heatmap figures (S5-S6), a complete
data inventory (Table S17) and a pathway-activity panel across subtypes
(Figure S7), which also corrected the Discussion statement about
cholesterol-synthesis-high subtypes. v24 additionally adds normal-vs-tumour
expression in a paired cohort (GSE15852; Figure S9, Table S19) and in the
larger TCGA-vs-GTEx comparison (Figure S11, Table S21), and
exploratory immunotherapy-response associations (GSE91061; Figure S10,
Table S20); TIDE scores are provided as exploratory results in
`results/v19_bonus/`. v24 is the version without the Xena confirmation;
v20 is the word-limited submission variant. The graphical
abstract, the updated study-design flowchart and the TRIPOD-AI checklist are
included under `manuscript/` and `results/figures/`; supplementary tables
S11-S21 are in `results/v19_bonus/`.

Figure files in `results/figures/` are named `figNN_*.png` (and 300-dpi
`figNN_*.tiff`), where NN is the figure number in the manuscript (Figures
1-16; `Supplementary_*.png` are Supplementary Figures S1-S3).

## Citation

If you use this pipeline, please cite the manuscript (DOI to be added) and the
primary data sources listed in `DATA_LICENSES.md`.
