# -*- coding: utf-8 -*-
"""ssGSEA immune enrichment scores (TCGA-BRCA) with the curated panels.

Implements the rank-based single-sample GSEA enrichment score (Barbie et al.
2009, alpha = 0.25) for the nine immune/stromal panels previously summarised
as mean-z marker scores, then compares subtypes and cross-validates against
the existing marker scores.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_ssgsea.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import OUT, FIG, save_patch

MARKERS = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG"],
    "CD4_T": ["CD4", "IL7R", "CCR7", "CD40LG", "ICOS"],
    "Treg": ["FOXP3", "CTLA4", "IL2RA", "IKZF2"],
    "NK": ["NKG7", "KLRD1", "KLRB1", "NCR1", "GNLY"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BLK"],
    "Macrophage": ["CD68", "CD163", "CSF1R", "LYZ", "ITGAX"],
    "DC": ["ITGAX", "CD1C", "CLEC9A", "BATF3", "FLT3"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "CEACAM8"],
    "Stroma": ["COL1A1", "COL3A1", "ACTA2", "PDGFRB", "FAP"],
}


def ssgsea_scores(X, gene_sets_idx, alpha=0.25):
    """Rank-based single-sample GSEA. X: (n_genes, n_samples) log expression.
    gene_sets_idx: dict name -> global row indices. Returns dict of 1-D arrays."""
    n_genes, n_samples = X.shape
    order = np.argsort(-X, axis=0, kind="stable")  # descending expression per sample
    scores = {}
    for name, gi in gene_sets_idx.items():
        gi = np.asarray(gi)
        weights = np.abs(X) ** alpha  # (genes, samples)
        denom = weights[gi].sum(axis=0) + 1e-12
        es = np.zeros(n_samples)
        for s in range(n_samples):
            ord_s = order[:, s]
            w_sorted = weights[ord_s, s]
            in_set = np.zeros(n_genes, dtype=bool)
            in_set[gi] = True
            in_sorted = in_set[ord_s]
            hit = np.cumsum(np.where(in_sorted, w_sorted, 0.0)) / denom[s]
            miss = np.cumsum(np.where(in_sorted, 0.0, 1.0)) / (n_genes - len(gi))
            es[s] = float(np.max(hit - miss))
        scores[name] = es
    return scores


def main():
    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")
    traits = pd.read_csv(OUT / "tcga_traits.csv")
    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    traits["sample_norm"] = traits["sample_id"].map(
        lambda s: "-".join(str(s).split("-")[:3]) + "-01" if str(s).startswith("TCGA-") else str(s))
    traits = traits.drop_duplicates(subset="sample_norm").set_index("sample_norm")

    idx = {g: i for i, g in enumerate(genes)}
    gene_sets_idx = {name: [idx[g] for g in gs if g in idx] for name, gs in MARKERS.items()}
    gene_sets_idx = {k: v for k, v in gene_sets_idx.items() if len(v) >= 3}
    print("panels:", {k: len(v) for k, v in gene_sets_idx.items()})

    scores = ssgsea_scores(X, gene_sets_idx)
    df = pd.DataFrame(scores, index=traits.index)
    df["subtype"] = sub.set_index(sub["sample_id"].map(
        lambda s: "-".join(str(s).split("-")[:3]) + "-01" if str(s).startswith("TCGA-") else str(s))
    ).reindex(df.index)["subtype"].astype(int)
    df.to_csv(OUT / "ssgsea_scores.csv", encoding="utf-8-sig")

    # subtype comparison
    ks = sorted(df["subtype"].unique())
    rows = {}
    for c in scores:
        groups = [df.loc[df["subtype"] == k, c] for k in ks]
        rows[c] = [float(stats.kruskal(*groups).pvalue)]
        rows[c] += [float(df.loc[df["subtype"] == k, c].mean()) for k in ks]
    summ = pd.DataFrame(rows, index=["kruskal_p"] + [f"C{k}" for k in ks]).T
    summ.to_csv(OUT / "ssgsea_by_subtype.csv", encoding="utf-8-sig")
    print(summ.round(4).to_string())

    # cross-check against mean-z marker scores
    try:
        old = pd.read_csv(OUT / "immune_scores.csv")
        corr = {}
        for c in scores:
            if c in old.columns:
                corr[c] = float(stats.spearmanr(df[c].to_numpy(), old[c].to_numpy())[0])
        print("marker-score vs ssGSEA spearman:", {k: round(v, 3) for k, v in corr.items()})
        with open(OUT / "ssgsea_vs_marker_corr.json", "w", encoding="utf-8") as f:
            json.dump(corr, f, indent=1)
    except Exception as e:
        print("marker cross-check skipped:", e)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mat = summ[[f"C{k}" for k in ks]].astype(float).T
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(mat.values, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(mat.index)), mat.index)
    ax.set_yticks(range(len(mat.columns)), mat.columns)
    ax.set_xlabel("Hub-gene molecular subtype")
    ax.set_title("ssGSEA immune scores by subtype (TCGA-BRCA)")
    fig.colorbar(im, ax=ax, label="mean ssGSEA ES", shrink=0.8)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig14_ssgsea_subtypes.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_patch({"kruskal_p": {k: float(v[0]) for k, v in rows.items()}},
               "ssgsea")
    print("DONE ->", FIG / "fig14_ssgsea_subtypes.png")


if __name__ == "__main__":
    main()
