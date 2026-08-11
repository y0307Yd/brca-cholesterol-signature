# -*- coding: utf-8 -*-
"""Round-2 bonus analyses for the content-complete manuscript (v22).

1. METABRIC immune-checkpoint genes by mapped subtype (Kruskal-Wallis).
2. GSE20711 supplementary figure (ROC, FDPS-by-ER, KM, subtype ER).
3. LASSO selection-frequency figure.
4. CIBERSORT subtype heatmap.
5. Data inventory table (Supplementary Table S17).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)

CHECK = ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "IDO1", "CD8A"]


def metabric_checkpoints():
    genes = [l.strip() for l in open(OUT / "mb_genes.txt", encoding="utf-8")
             if l.strip()]
    X = np.load(OUT / "mb_X.npy")
    mb = pd.read_csv(OUT / "metabric_subtypes_mapped.csv")
    idx = {g: i for i, g in enumerate(genes)}
    rows = []
    for g in CHECK:
        x = X[idx[g]]
        med = {}
        groups = []
        for k in range(1, 5):
            mask = (mb["subtype"] == k).to_numpy()
            med[k] = float(np.nanmedian(x[mask]))
            groups.append(x[mask])
        kw = stats.kruskal(*groups)
        rows.append({
            "gene": g, "median_C1": med[1], "median_C2": med[2],
            "median_C3": med[3], "median_C4": med[4],
            "kruskal_P": kw.pvalue,
            "highest_subtype": max(med, key=med.get),
        })
    res = pd.DataFrame(rows)
    res.to_csv(WORK / "bonus_metabric_checkpoint_by_subtype.csv",
               index=False, encoding="utf-8-sig")
    res.to_csv(OUT / "Supplementary_Table_S18_metabric_checkpoints.csv",
               index=False, encoding="utf-8-sig")
    print("=== METABRIC immune checkpoints by subtype ===")
    print(res.round(4).to_string(index=False))
    return res


def gse20711_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    from sklearn.metrics import roc_curve
    m = pd.read_csv(OUT / "geo_gse20711_patient_data.csv", index_col=0)
    sub = pd.read_csv(OUT / "geo_gse20711_subtypes_mapped.csv", index_col=0)
    m = m.join(sub[["subtype"]])
    y_er = m["ER"].dropna()
    p_er = m.loc[y_er.index, "ER_pred"]
    from online_utils import auc_ci
    auc, ci = auc_ci(y_er.astype(int).values, p_er.values)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8.4))
    fpr, tpr, _ = roc_curve(y_er.astype(int).values, p_er.values)
    axes[0, 0].plot(fpr, tpr, lw=2, color="#1f77b4")
    axes[0, 0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0, 0].set_title(f"ER classifier transfer (AUC={auc:.3f}, "
                         f"95%CI {ci[0]:.3f}-{ci[1]:.3f})")
    axes[0, 0].set_xlabel("1 - Specificity")
    axes[0, 0].set_ylabel("Sensitivity")
    axes[0, 0].grid(alpha=0.3)

    bp = axes[0, 1]
    d0 = m.loc[m["ER"] == 0, "FDPS_int"].dropna().values
    d1 = m.loc[m["ER"] == 1, "FDPS_int"].dropna().values
    bp.boxplot([d0, d1], tick_labels=["ER-", "ER+"], widths=0.5)
    bp.scatter([1] * len(d0), d0, s=8, alpha=0.35, color="#1f77b4")
    bp.scatter([2] * len(d1), d1, s=8, alpha=0.35, color="#d62728")
    from scipy import stats as st
    p = st.mannwhitneyu(d0, d1).pvalue
    bp.set_title(f"FDPS by ER status (MWU P={p:.2g})")
    bp.set_ylabel("FDPS expression (INT)")
    bp.grid(alpha=0.3, axis="y")

    med = np.median(m["FDPS_int"])
    grp = (m["FDPS_int"] > med).astype(int)
    ax = axes[1, 0]
    for name, g in [("FDPS low", grp == 0), ("FDPS high", grp == 1)]:
        kmf = KaplanMeierFitter()
        kmf.fit(m.loc[g, "time_months"], m.loc[g, "event"])
        ax.step(kmf.timeline, kmf.survival_function_["KM_estimate"],
                where="post", label=name, lw=1.8)
    from lifelines.statistics import logrank_test
    lr = logrank_test(m.loc[grp == 0, "time_months"],
                      m.loc[grp == 1, "time_months"],
                      m.loc[grp == 0, "event"], m.loc[grp == 1, "event"])
    ax.set_title(f"RFS by FDPS median (log-rank P={lr.p_value:.2f})")
    ax.set_xlabel("Months")
    ax.set_ylabel("RFS probability")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    frac = [m.loc[m["subtype"] == k, "ER"].mean() if (m["subtype"] == k).any()
            else 0 for k in [1, 2, 3, 4]]
    n = [int((m["subtype"] == k).sum()) for k in [1, 2, 3, 4]]
    ax.bar(["C1", "C2", "C3", "C4"], frac,
           color=["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"])
    for i, (f, nn) in enumerate(zip(frac, n)):
        ax.text(i, f + 0.02, f"n={nn}", ha="center", fontsize=9)
    ax.set_title("Mapped-subtype ER composition (exploratory)")
    ax.set_ylabel("Fraction ER+")
    ax.set_ylim(0, 1.1)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("GSE20711 validation (n = 88 with RFS)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "supp_fig4_gse20711_validation.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig4_gse20711_validation.png")


def lasso_frequency_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    freq = pd.read_csv(WORK / "bonus_lasso_stability.csv")
    hubs = set(pd.read_csv(OUT / "hub_genes.csv")["gene"])
    freq = freq.sort_values("selection_frequency")
    cols = ["#c00000" if g in hubs else "#6b8eb8"
            for g in freq["gene"]]
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.barh(freq["gene"], freq["selection_frequency"], color=cols, height=0.7)
    ax.axvline(0.5, color="grey", ls="--", lw=1)
    ax.set_xlabel("Selection frequency (300 LASSO bootstrap resamples)")
    ax.set_title("LASSO selection stability; red = 11 hub genes")
    ax.set_xlim(0, 1.02)
    ax.grid(alpha=0.3, axis="x")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#c00000", label="Hub gene"),
                       Patch(color="#6b8eb8", label="Other candidate")],
              frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "supp_fig5_lasso_stability.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig5_lasso_stability.png")


def cibersort_heatmap():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list
    c = pd.read_csv(WORK / "bonus_cibersort_by_subtype.csv")
    mat = c.set_index("cell")[[f"median_C{k}" for k in [1, 2, 3, 4]]].values
    cells = c["cell"].tolist()
    Z = linkage(mat, method="average")
    order = leaves_list(Z)
    mat = mat[order]
    cells = [cells[i] for i in order]
    fig, ax = plt.subplots(figsize=(7, 10))
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu")
    ax.set_yticks(range(len(cells)))
    ax.set_yticklabels(cells, fontsize=8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["C1", "C2", "C3", "C4"], fontsize=10)
    ax.set_title("CIBERSORT median fractions by subtype")
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label("median fraction")
    fig.tight_layout()
    fig.savefig(FIG / "supp_fig6_cibersort_heatmap.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig6_cibersort_heatmap.png")


def data_inventory():
    rows = [
        ("TCGA-BRCA expression/clinical/methylation", "NCI GDC portal",
         "STAR-counts, clinical XML, HM450", "MD5-verified (4,393/4,393 files)"),
        ("TCGA-CDR outcomes", "Liu et al. 2018", "PFI", "Public"),
        ("METABRIC", "cBioPortal DataHub", "Illumina expression + RFS",
         "Curtis et al. 2012"),
        ("PAM50 calls (TCGA)", "UCSC Xena GDC hub", "PAM50 labels", "Public"),
        ("FDPS HM450 methylation", "cBioPortal brca_tcga_methylation_hm450",
         "beta values", "Public"),
        ("GSE21653", "NCBI GEO (GPL570)", "series matrix + GPL570.annot",
         "Sabatier et al. 2011"),
        ("GSE7390", "ArrayExpress E-GEOD-7390 (GPL96)", "processed matrix + SDRF",
         "Desmedt et al. 2007"),
        ("GSE20711", "NCBI GEO (GPL570)", "series matrix",
         "Fackler et al. 2011"),
        ("GSE176078", "CELLxGENE", "100,064-cell atlas",
         "Wu et al. 2021"),
        ("GSE161529", "Mendeley Data mirror", "70,419-cell atlas",
         "Pal et al. 2021"),
        ("eQTLGen cis-eQTLs", "eQTLGen consortium", "SMR-format BESD",
         "Vosa et al. 2021"),
        ("GTEx v8", "Yang Lab SMR database / GTEx Portal",
         "breast-mammary + whole-blood BESD", "GTEx Consortium"),
        ("BCAC overall/TNBC GWAS", "GWAS Catalog GCST010098/GCST010100",
         "SMR-format .ma", "Zhang et al. 2020"),
        ("BCAC OncoArray ER+/ER-", "GWAS Catalog GCST004988",
         "ER-stratified summary statistics", "Michailidou et al. 2017"),
        ("1000 Genomes Phase 3 EUR", "1000 Genomes Project", "LD reference",
         "Public"),
        ("LM22 signature matrix", "Newman et al. 2015 (CIBERSORTx mirror)",
         "22 immune cell types", "Public"),
        ("HPA FDPS entry", "Human Protein Atlas ENSG00000160752",
         "protein/RNA evidence", "Uhlen et al. 2015"),
        ("Pan-cancer FDPS RNA-seq", "cBioPortal datahub", "17 TCGA studies",
         "Cerami et al. 2012"),
    ]
    df = pd.DataFrame(rows, columns=["Dataset", "Source", "Content", "Reference/use"])
    df.to_csv(OUT / "Supplementary_Table_S17_data_inventory.csv",
              index=False, encoding="utf-8-sig")
    print("saved Supplementary_Table_S17_data_inventory.csv")


if __name__ == "__main__":
    metabric_checkpoints()
    gse20711_figure()
    lasso_frequency_figure()
    cibersort_heatmap()
    data_inventory()
