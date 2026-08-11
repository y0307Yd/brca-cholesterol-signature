# -*- coding: utf-8 -*-
"""Graphical abstract for the cholesterol-metabolism breast cancer study."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(14, 6.2))
ax.set_xlim(0, 140)
ax.set_ylim(0, 62)
ax.axis("off")


def box(x, y, w, h, title, lines, fc="#ffffff", ec="#1f4e79", tfs=10.5, lfs=8,
        tc="#1f4e79"):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5",
                       fc=fc, ec=ec, lw=1.6)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h - 4.5, title, ha="center", va="center",
            fontsize=tfs, fontweight="bold", color=tc)
    ax.text(x + w / 2, y + h - 12, lines, ha="center", va="top",
            fontsize=lfs, color="#222222", linespacing=1.45)


def arrow(x1, y1, x2, y2):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                        mutation_scale=22, lw=2.4, color="#1f4e79")
    ax.add_patch(a)


# Row 1: discovery
box(2, 44, 24, 16, "Public multi-omics",
     "TCGA-BRCA (n=952)\nMETABRIC (n=1,979)\nGEO (n=777)\nGSE176078/GSE161529 scRNA",
     fc="#eef6ee")
box(32, 44, 24, 16, "Targeted screening",
     "Cholesterol gene set (120)\nWGCNA + pathway-first\n119 measurable genes",
     fc="#fff6e5")
box(62, 44, 24, 16, "Machine learning",
     "LASSO + SVM-RFE\n11 hub genes\n(ER classifier AUC 0.946)",
     fc="#eaf1fb")
box(92, 44, 24, 16, "Patient subtypes",
     "4 consensus subtypes\nER composition + immune\npatterns externally reproduced",
     fc="#f3eafb")
arrow(26, 52, 31, 52)
arrow(56, 52, 61, 52)
arrow(86, 52, 91, 52)

# Row 2: validation layers
box(2, 20, 24, 16, "Single-cell localization",
     "Tumour epithelium,\nmyeloid/macrophage\nand endothelial cells",
     fc="#eaf1fb")
box(32, 20, 24, 16, "Immune deconvolution",
     "CIBERSORT/ESTIMATE/\nssGSEA; immune-checkpoint\ngenes; C1 immune-hot",
     fc="#eaf1fb")
box(62, 20, 24, 16, "Causal inference",
     "SMR/HEIDI + coloc.abf\nFDPS protective\n(P=9.8e-7; PP.H4=0.9985)",
     fc="#fff6e5")
box(92, 20, 24, 16, "External verification",
     "4 independent cohorts\nAUC 0.85-0.90\nHPA protein + pan-cancer",
     fc="#eef6ee")
arrow(26, 36, 31, 36)
arrow(56, 36, 61, 36)
arrow(86, 36, 91, 36)

# central arrows row1 -> row2
arrow(14, 44, 14, 37)
arrow(44, 44, 44, 37)
arrow(74, 44, 74, 37)
arrow(104, 44, 104, 37)

# bottom conclusion
box(20, 2, 100, 12, "Conclusion",
     "Cholesterol genes carry a strong cross-platform ER signature and define immune-relevant subtypes; "
     "germline-regulated FDPS expression is a putatively causal protective factor for breast cancer risk.",
     fc="#ffffff", ec="#c00000", tfs=11, lfs=9, tc="#c00000")

fig.suptitle("Cholesterol metabolism in breast cancer: a multi-cohort, multi-omics evidence chain",
             fontsize=14, fontweight="bold", y=0.98)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(r"outputs\chol_metab_signature\figures\graphical_abstract.png",
            dpi=300, bbox_inches="tight")
fig.savefig(r"outputs\chol_metab_signature\figures\graphical_abstract.pdf",
            bbox_inches="tight")
print("saved graphical abstract")
