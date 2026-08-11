# -*- coding: utf-8 -*-
"""Build the comparison table of published lipid-metabolism BC signatures."""
import pandas as pd

OUT = r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature"

rows = [
    {
        "Study": "Ma et al., BMC Med Genomics, 2025 (PMID 40745642)",
        "Signature": "Lipid metabolism-related prognostic signature (ER+ tamoxifen)",
        "n_genes": "Not reported in abstract (ML-derived)",
        "Cohorts": "TCGA-ER+BRCA (training); GSE17705/GSE22219/GSE42568/GSE58644 (validation)",
        "Selection": "Lipid-metabolism gene set + machine learning (Cox/LASSO)",
        "Validation": "Multiple independent cohorts; 5-y OS AUC 0.858",
        "Immune": "Immune infiltration comparison (high vs low risk)",
        "Single_cell": "Yes (TME heterogeneity)",
        "Subtyping": "High/low risk only",
        "Causal_inference": "No",
        "Protein_level": "No",
    },
    {
        "Study": "Zhao et al., Aging, 2024 (PMID 38771140)",
        "Signature": "Lipid metabolism + immune-related 9-gene signature",
        "n_genes": "9 (CALR, CCL5, CEPT, FTT3, CXCL13, FLT3, IL12B, IL18, IL24)",
        "Cohorts": "UCSC-derived DEGs (TCGA-based); RT-PCR validation of IL18",
        "Selection": "Univariate Cox + LASSO on lipid/immune DEGs",
        "Validation": "Prognostic KM/ROC; IL18 RT-PCR",
        "Immune": "Immune-related functional mining",
        "Single_cell": "No",
        "Subtyping": "High/low risk only",
        "Causal_inference": "No",
        "Protein_level": "RT-PCR (mRNA) only",
    },
    {
        "Study": "Mai et al., Sci Rep, 2025 (PMID 40715365)",
        "Signature": "Lipid metabolism-related 6-gene signature",
        "n_genes": "6 (APOC3, CEL, CPT1A, JAK2, NFKBIA, PLA2G1B)",
        "Cohorts": "TCGA (1:1 train/test) + external validation cohorts",
        "Selection": "OS-associated LMRGs + LASSO + stepwise Cox",
        "Validation": "Training/test/external cohorts, clinical subgroups",
        "Immune": "Not detailed in abstract",
        "Single_cell": "No",
        "Subtyping": "High/low risk only",
        "Causal_inference": "No",
        "Protein_level": "No",
    },
    {
        "Study": "Shen et al., Front Immunol, 2023 (PMID 37469520)",
        "Signature": "Lipid metabolism-related 9-gene signature (ER+)",
        "n_genes": "9 survival-related LMRGs",
        "Cohorts": "TCGA (training); METABRIC, GEO, own cohort (validation)",
        "Selection": "Consensus clustering of LMRG patterns + LASSO",
        "Validation": "External validation; nomogram",
        "Immune": "Immune landscapes; immunotherapy/chemotherapy response",
        "Single_cell": "No",
        "Subtyping": "2 molecular patterns",
        "Causal_inference": "No",
        "Protein_level": "IHC of hub gene(s)",
    },
    {
        "Study": "Gong et al., Int J Gen Med, 2021 (PMID 34916832)",
        "Signature": "Lipid metabolism-associated 16-gene signature",
        "n_genes": "16",
        "Cohorts": "TCGA-BRCA (n=1053)",
        "Selection": "DEG + univariate Cox + LASSO",
        "Validation": "Internal only",
        "Immune": "ssGSEA correlation",
        "Single_cell": "No",
        "Subtyping": "High/low risk only",
        "Causal_inference": "No",
        "Protein_level": "No",
    },
    {
        "Study": "Chang & Xing, Lipids Health Dis, 2022 (PMID 35562758)",
        "Signature": "FABP7 + NDUFAB1 (from 16 LRGs)",
        "n_genes": "2",
        "Cohorts": "TCGA-BRCA (training)",
        "Selection": "LASSO Cox; consensus clustering",
        "Validation": "Internal; immune/TMB association",
        "Immune": "Immune infiltration, TMB/MSI",
        "Single_cell": "No",
        "Subtyping": "Consensus clusters",
        "Causal_inference": "No",
        "Protein_level": "No",
    },
    {
        "Study": "This study (v18+)",
        "Signature": "Cholesterol metabolism 11-hub-gene signature + FDPS",
        "n_genes": "11 hub genes (+FDPS)",
        "Cohorts": "TCGA-BRCA; METABRIC, GSE21653, GSE7390 (validation)",
        "Selection": "WGCNA + cholesterol gene set intersection; LASSO + SVM-RFE",
        "Validation": "3 external expression cohorts; ER-classifier AUC 0.85-0.90; subtype external mapping",
        "Immune": "CIBERSORT, ESTIMATE, ssGSEA, immune-checkpoint genes by subtype",
        "Single_cell": "Yes (GSE161529, monocyte/macrophage localization)",
        "Subtyping": "4 molecular subtypes + external validation",
        "Causal_inference": "Yes (SMR + HEIDI + Bayesian coloc, eQTLGen x BCAC, ER-stratified)",
        "Protein_level": "HPA protein evidence (FDPS)",
    },
]

df = pd.DataFrame(rows)
path = OUT + r"\Supplementary_Table_S11_lipid_signature_comparison.csv"
df.to_csv(path, index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print("saved", path)
