# GitHub upload manifest: v25 to v28 increment

Baseline: commit `88e4e0b` (v25, 2026-08-11).

## Include in the GitHub commit

- `code/complete_v25_additions.py`: corrected TIDE mapping and adjusted models.
- `code/add_v27_robustness.py`: random forest, negative controls and AUC meta-analysis.
- `code/add_v28_stability_immune.py`: subtype stability and adjusted immune analysis.
- `results/v19_bonus/Supplementary_Table_S22_TIDE_subtype_sensitivity.csv`.
- `results/v19_bonus/Supplementary_Table_S23_random_forest_robustness.csv`.
- `results/v19_bonus/Supplementary_Table_S24_classifier_negative_controls.csv`.
- `results/v19_bonus/Supplementary_Table_S25_external_AUC_meta_analysis.csv`.
- `results/v19_bonus/Supplementary_Table_S26_subtype_stability.csv`.
- `results/v19_bonus/Supplementary_Table_S26_subtype_stability_resamples.csv`.
- `results/v19_bonus/Supplementary_Table_S27_adjusted_immune_models.csv`.
- Corrected sample-level TIDE outputs and v27/v28 machine-readable summaries.
- Corrected Supplementary Figure S8 and new Supplementary Figures S12-S14 in PNG/PDF/SVG.
- `manuscript/Manuscript_Cholesterol_Metabolism_BRCA_v28.pdf`.
- Updated repository and manuscript README files.

## Exclude from GitHub

- Editable manuscript DOCX files.
- Intermediate v26 and v27 manuscript snapshots.
- Supplementary TIFF files S12-S14 (55-101 MB each); archive these in Zenodo.
- Raw expression matrices, consensus arrays and downloaded source datasets.
- QA thumbnails, contact sheets, LibreOffice profiles and local configuration files.

## Before pushing

1. Review `git diff --stat` and `git status --short`.
2. Confirm no file exceeds the GitHub 100 MB hard limit.
3. Run `python -m py_compile` on the three new Python scripts.
4. Update the Zenodo record with the v28 PDF, production TIFF files and a source snapshot.
