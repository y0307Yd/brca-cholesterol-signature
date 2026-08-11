# -*- coding: utf-8 -*-
"""Study-design flow diagram v19 (Supplementary Figure S2, updated)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(11, 9.5))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")


def box(x, y, w, h, text, fc="#eef3fb", ec="#1f4e79", fs=7.5, bold=False):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                       fc=fc, ec=ec, lw=1.1)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight="bold" if bold else "normal")


def arrow(x1, y1, x2, y2, color="#1f4e79"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                        mutation_scale=12, lw=1.2, color=color)
    ax.add_patch(a)


# Row 1: data sources
box(2, 90, 22, 8, "TCGA-BRCA (n = 952; 122 PFI events)\nSTAR counts, MD5-verified (GDC)",
    fc="#e2f0d9", fs=7)
box(26, 90, 22, 8, "METABRIC (n = 1,979; 803 RFS events)\nIllumina microarray (cBioPortal)",
    fc="#e2f0d9", fs=7)
box(50, 90, 22, 8, "GEO validation: GSE21653 (n = 252),\nGSE7390 (n = 198)",
    fc="#e2f0d9", fs=7)
box(74, 90, 24, 8, "Cholesterol gene set (120 genes)\nKEGG hsa00100 + GO:0008203/0006695",
    fc="#fff2cc", fs=7)

# Row 2: WGCNA + pathway
box(2, 74, 22, 8, "WGCNA (8,000 most variable genes)\n15 modules; PFI module r = 0.082", fs=7)
box(26, 74, 22, 8, "Strict module \u2229 pathway intersection\n= 0 genes (documented null)",
    fc="#fbe5e5", ec="#c00000", fs=7)
box(50, 74, 22, 8, "Pathway-first adaptation\n119 measurable cholesterol genes",
    fc="#fff2cc", fs=7)
box(74, 74, 24, 8, "External ER/survival validation:\nMETABRIC + 2 GEO cohorts (AUC 0.85-0.90)",
    fc="#e2f0d9", fs=7)
arrow(13, 90, 13, 84)
arrow(37, 90, 37, 84)
arrow(61, 86, 61, 80)
arrow(31, 78, 47, 76)

# Row 3: ML
box(2, 58, 22, 8, "LASSO (C = 0.316) \u2192 69 genes\nSVM-RFE \u2192 12 genes\nintersection = 11 hub genes\n(bootstrap selection \u2265 73.7%)",
    fs=6.8)
box(26, 58, 22, 8, "ER-status classifier (11 genes)\nTCGA 5-fold CV AUC 0.946 \u00b1 0.012\nMETABRIC AUC 0.922",
    fc="#e2f0d9", fs=7)
box(50, 58, 22, 8, "Calibration: Brier 0.069 (TCGA CV)\n0.098 (METABRIC); DCA: positive net benefit",
    fc="#e2f0d9", fs=7)
box(74, 58, 24, 8, "CIBERSORT (22 immune cell types) +\nESTIMATE (immune/stromal scores)",
    fc="#dbe5f1", fs=7)
arrow(24, 62, 26, 62)
arrow(48, 62, 50, 62)
arrow(72, 62, 74, 62)
arrow(13, 74, 13, 68)

# Row 4: subtypes + scRNA
box(2, 42, 22, 8, "Consensus clustering (k = 4; PAC = 0.097)\nC1 ER-low/immune-hot; C3 ER-high/immune-cold",
    fs=6.8)
box(26, 42, 22, 8, "Subtype external validation\nMETABRIC + GSE21653 + GSE7390\nER composition reproduced; immune-checkpoint genes",
    fc="#e2f0d9", fs=6.8)
box(50, 42, 22, 8, "Single-cell (GSE176078 + GSE161529)\nhub genes in tumour epithelium, myeloid/\nmacrophage and endothelial cells",
    fs=6.8)
box(74, 42, 24, 8, "Immune deconvolution by subtype\nC1: highest CD8/Treg/NK/M1/DC signatures;\nC3: immune-cold (all P < 0.05)",
    fc="#dbe5f1", fs=6.8)
arrow(13, 58, 13, 52)
arrow(37, 58, 37, 52)
arrow(61, 58, 61, 52)
arrow(86, 58, 86, 52)

# Row 5: causal inference
box(2, 26, 22, 8, "SMR/HEIDI (eQTLGen blood \u00d7 BCAC)\nFDPS: P = 9.8e-07, P_HEIDI = 0.026\nprotective (b_SMR = -0.42)",
    fc="#fff2cc", fs=7)
box(26, 26, 22, 8, "coloc.abf: PP.H4 = 0.9985 (3248 SNPs)\n+ LIPA/FAXDC2/SREBF1 coloc\n(suggestive; SuSiE caveat reported)",
    fc="#fff2cc", fs=6.8)
box(50, 26, 22, 8, "ER-stratified SMR\nER+ P = 1.54e-05 (b = -0.44)\nTNBC sensitivity; GTEx sensitivity",
    fc="#fff2cc", fs=7)
box(74, 26, 24, 8, "Protein/methylation/pan-cancer\nHPA tissue-enhanced FDPS; HM450 methylation\nrho = -0.187 (ER-adjusted); 17 cancers expressed",
    fc="#e2f0d9", fs=6.8)
arrow(13, 42, 13, 36)
arrow(37, 42, 37, 36)
arrow(61, 42, 61, 36)
arrow(86, 42, 86, 36)

# Row 6: clinical + honesty box
box(2, 10, 30, 8, "FDPS clinical association (honest negative)\nER- high expression; not an independent prognostic factor",
    fc="#fbe5e5", ec="#c00000", fs=7)
box(36, 10, 28, 8, "GTEx v8 breast sensitivity\nLIPA/FAXDC2 nominal direction-consistent;\nFDPS lacks breast instrument (power < 0.01%)",
    fs=6.8)
box(68, 10, 30, 8, "Reporting: TRIPOD checklist; code on GitHub;\ndata/scripts archived on Zenodo (DOI 10.5281/zenodo.21873168)",
    fc="#e2f0d9", fs=6.8)
arrow(13, 26, 13, 20)
arrow(50, 26, 50, 20)
arrow(86, 26, 86, 20)

fig.suptitle("Study design: cholesterol metabolism in breast cancer "
             "(multi-cohort, multi-omics evidence chain)", fontsize=12.5, y=0.985)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(r"outputs\chol_metab_signature\figures\supp_fig2_design_flowchart_v19.png",
            dpi=300, bbox_inches="tight")
print("saved supp_fig2 v19")
