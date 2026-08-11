# -*- coding: utf-8 -*-
"""TIDE immune-escape / immunotherapy-response features across subtypes."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)


def main():
    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8")
             if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")          # genes x samples
    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    expr = pd.DataFrame(X, index=genes,
                        columns=[f"S{i}" for i in range(X.shape[1])])
    print("expression:", expr.shape)

    from tidepy.pred import TIDE
    res = TIDE(expr, cancer="Other")
    print("TIDE output shape:", res.shape)
    print(res.head(3).to_string())
    res.to_csv(OUT / "tide_by_sample.csv", encoding="utf-8-sig")

    # merge subtype
    df = res.copy()
    df["subtype"] = sub["subtype"].to_numpy()
    ks = [1, 2, 3, 4]
    feats = ["TIDE", "IFNG", "Dysfunction", "Exclusion", "MDSC", "CAF",
             "TAM M2", "CD274", "CD8", "CTL"]
    rows = []
    for f in feats:
        groups = [df.loc[df["subtype"] == k, f].to_numpy() for k in ks]
        kw = stats.kruskal(*groups)
        rows.append({
            "feature": f,
            "kruskal_P": kw.pvalue,
            **{f"median_C{k}": float(np.median(g)) for k, g in zip(ks, groups)},
        })
    tide_sum = pd.DataFrame(rows)
    tide_sum.to_csv(OUT / "tide_by_subtype.csv", index=False,
                    encoding="utf-8-sig")
    print(tide_sum.round(4).to_string(index=False))

    # responder rates per subtype
    rr = df.groupby("subtype")["Responder"].mean()
    print("\nResponder rate by subtype:")
    print(rr.round(3).to_string())
    chi = stats.chi2_contingency(
        pd.crosstab(df["subtype"], df["Responder"]).to_numpy())
    print("chi-square P:", chi.pvalue)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    cols = {"1": "#d62728", "2": "#1f77b4", "3": "#2ca02c", "4": "#ff7f0e"}
    for ax, f in zip(axes.ravel(), ["TIDE", "IFNG", "Dysfunction",
                                    "Exclusion", "CD274", "CD8"]):
        data = [df.loc[df["subtype"] == k, f].dropna() for k in ks]
        bp = ax.boxplot(data, tick_labels=[f"C{k}" for k in ks], widths=0.5)
        for patch, k in zip(bp["boxes"], ks):
            patch.set_color(cols[str(k)])
        kw = stats.kruskal(*[d.to_numpy() for d in data])
        ax.set_title(f"{f} (Kruskal P={kw.pvalue:.2g})")
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle("TIDE immunotherapy-response features by cholesterol subtype",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "supp_fig8_tide_subtypes.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig8_tide_subtypes.png")


if __name__ == "__main__":
    main()
