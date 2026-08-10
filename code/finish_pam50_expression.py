# -*- coding: utf-8 -*-
"""Hub-gene / FDPS expression by PAM50 intrinsic subtype (TCGA-BRCA).

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_pam50_expression.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import OUT, FIG, save_patch

PAM50_ORDER = ["Basal", "HER2", "LumA", "LumB", "Normal"]
HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR", "FDPS"]


def norm_sample(s):
    s = str(s).strip()
    if s.startswith("TCGA-"):
        parts = s.split("-")
        if len(parts) >= 4:
            return "-".join(parts[:3]) + "-01"
    return s


def main():
    pam_path = Path("data/pam50/BRCA_clinicalMatrix.tsv")
    pam = pd.read_csv(pam_path, sep="\t", dtype=str)
    col = "PAM50Call_RNAseq"
    s = pd.Series(pam[col].values, index=pam["sampleID"].map(norm_sample)).dropna()
    alias = {
        "basal": "Basal", "basal-like": "Basal", "her2": "HER2",
        "her2-enriched": "HER2", "luma": "LumA", "luminal a": "LumA",
        "lumb": "LumB", "luminal b": "LumB", "normal": "Normal",
        "normal-like": "Normal",
    }
    s = s.map(lambda v: alias.get(str(v).strip().lower().replace("-", " ").replace("_", " "),
                                  str(v).strip()))
    s = s[s.isin(PAM50_ORDER)]
    s = s[~s.index.duplicated(keep="first")]
    print("PAM50 samples:", len(s))
    print(s.value_counts().to_string())

    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")
    gidx = {g: genes.index(g) for g in HUB if g in genes}
    missing = [g for g in HUB if g not in gidx]
    print("missing genes:", missing)

    traits = pd.read_csv(OUT / "tcga_traits.csv")
    traits["sample_norm"] = traits["sample_id"].map(norm_sample)
    traits = traits.drop_duplicates(subset="sample_norm").set_index("sample_norm")
    rows = []
    for g, i in gidx.items():
        e = pd.Series(X[i], index=traits.index)
        d = pd.concat([e.rename("expr"), s.rename("PAM50")], axis=1, join="inner").dropna()
        groups = [d.loc[d["PAM50"] == p, "expr"] for p in PAM50_ORDER]
        h, p_kw = stats.kruskal(*groups)
        means = {p: float(d.loc[d["PAM50"] == p, "expr"].mean()) for p in PAM50_ORDER}
        ns = {f"n_{p}": int((d["PAM50"] == p).sum()) for p in PAM50_ORDER}
        per_subtype_p = {}
        for p in PAM50_ORDER:
            other = d.loc[d["PAM50"] != p, "expr"]
            this = d.loc[d["PAM50"] == p, "expr"]
            per_subtype_p[p] = float(stats.mannwhitneyu(this, other).pvalue) \
                if len(this) > 5 and len(other) > 5 else np.nan
        rows.append({"gene": g, "kruskal_H": float(h), "kruskal_p": float(p_kw), **means, **ns})

    res = pd.DataFrame(rows).set_index("gene")
    res.to_csv(OUT / "hub_genes_pam50_expression.csv", encoding="utf-8-sig")
    print(res.round(3).to_string())

    # BH-FDR across genes x subtypes for the MWU tests
    pvals = []
    for g in gidx:
        e = pd.Series(X[gidx[g]], index=traits.index)
        d = pd.concat([e.rename("expr"), s.rename("PAM50")], axis=1, join="inner").dropna()
        for p in PAM50_ORDER:
            this = d.loc[d["PAM50"] == p, "expr"]
            other = d.loc[d["PAM50"] != p, "expr"]
            if len(this) > 5 and len(other) > 5:
                pvals.append((g, p, stats.mannwhitneyu(this, other).pvalue))
    pv = [x[2] for x in pvals]
    m = len(pv)
    order = np.argsort(pv)
    q = np.empty(m)
    running = 1.0
    for rank, i in enumerate(reversed(order)):
        running = min(1.0, pv[i] * m / (m - rank))
        q[i] = running
    q = list(q)
    fdr_df = pd.DataFrame(pvals, columns=["gene", "PAM50", "p"])
    fdr_df["FDR"] = q
    fdr_df.to_csv(OUT / "hub_genes_pam50_mwu_fdr.csv", index=False, encoding="utf-8-sig")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mat = res[[p for p in PAM50_ORDER]].astype(float)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(mat.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(PAM50_ORDER)), PAM50_ORDER)
    ax.set_yticks(range(len(mat.index)), mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Hub genes and FDPS expression by PAM50 subtype (TCGA-BRCA)")
    fig.colorbar(im, ax=ax, label="mean log1p expression", shrink=0.8)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig13_hub_genes_pam50_expression.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "n_pam50": int(len(s)),
        "per_gene_kruskal_p": {g: float(res.loc[g, "kruskal_p"]) for g in res.index},
        "fdps_by_pam50": {p: float(res.loc["FDPS", p]) for p in PAM50_ORDER},
    }
    save_patch(summary, "pam50_expr")
    print("DONE ->", FIG / "fig13_hub_genes_pam50_expression.png")


if __name__ == "__main__":
    main()
