# -*- coding: utf-8 -*-
"""ER-adjusted FDPS methylation-expression analysis (TCGA-BRCA).

Extends the v13 methylation result (rho = -0.203) with:
  1) rank-based partial correlation controlling for ER status;
  2) ER+ and ER- stratified Spearman correlations with bootstrap CIs;
  3) updated figure (Fig 12b) and Supplementary Table S7 columns.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_methyl_er_adjusted.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import OUT, FIG, save_patch

DL = Path(__file__).resolve().parent.parent / "data" / "methyl"


def norm_sample(s):
    s = str(s).strip()
    if s.startswith("TCGA-"):
        parts = s.split("-")
        if len(parts) >= 4:
            return "-".join(parts[:3]) + "-01"
    return s


def spearman_ci(x, y, seed=20260810, n_boot=2000):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = []
    n = len(x)
    for _ in range(n_boot):
        ii = rng.integers(0, n, n)
        boots.append(stats.spearmanr(x[ii], y[ii])[0])
    ci = np.percentile(boots, [2.5, 97.5])
    return float(rho), float(p), float(ci[0]), float(ci[1])


def partial_spearman(x, y, z):
    """Partial Spearman correlation of x,y controlling for z (rank-transform)."""
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    # residuals of ranks on z
    def resid(u):
        A = np.vstack([np.ones(len(u)), rz]).T
        beta, *_ = np.linalg.lstsq(A, u, rcond=None)
        return u - A @ beta

    ex, ey = resid(rx), resid(ry)
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p)


def main():
    path = DL / "cbioportal_fdps_methylation.json"
    arr = json.loads(path.read_text(encoding="utf-8"))
    meth = pd.DataFrame(
        [{"sample_id": norm_sample(x["sampleId"]), "beta": float(x["value"])}
         for x in arr if x.get("value") is not None]
    ).drop_duplicates(subset="sample_id").dropna()

    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")
    traits = pd.read_csv(OUT / "tcga_traits.csv")
    gidx = genes.index("FDPS")
    expr = pd.DataFrame({
        "sample_id": traits["sample_id"].map(norm_sample),
        "expr": X[gidx],
        "ER": traits["ER"].to_numpy(),
    })
    df = meth.merge(expr, on="sample_id").dropna(subset=["ER"])
    print("matched samples with ER status:", len(df))
    if len(df) < 200:
        raise SystemExit("too few matched samples")

    # 1) overall
    rho_all, p_all, lo_all, hi_all = spearman_ci(df["expr"], df["beta"])
    # 2) partial controlling ER
    rho_part, p_part = partial_spearman(df["expr"], df["beta"], df["ER"])
    # 3) strata
    rho_pos, p_pos, lo_pos, hi_pos = spearman_ci(
        df.loc[df["ER"] == 1, "expr"], df.loc[df["ER"] == 1, "beta"])
    rho_neg, p_neg, lo_neg, hi_neg = spearman_ci(
        df.loc[df["ER"] == 0, "expr"], df.loc[df["ER"] == 0, "beta"])

    print(f"overall rho={rho_all:.3f} ({lo_all:.3f}-{hi_all:.3f}) P={p_all:.3e} n={len(df)}")
    print(f"partial (ER-adjusted) rho={rho_part:.3f} P={p_part:.3e}")
    print(f"ER+ rho={rho_pos:.3f} ({lo_pos:.3f}-{hi_pos:.3f}) P={p_pos:.3e} "
          f"n={(df['ER']==1).sum()}")
    print(f"ER- rho={rho_neg:.3f} ({lo_neg:.3f}-{hi_neg:.3f}) P={p_neg:.3e} "
          f"n={(df['ER']==0).sum()}")

    summary = {
        "n": int(len(df)),
        "n_er_pos": int((df["ER"] == 1).sum()),
        "n_er_neg": int((df["ER"] == 0).sum()),
        "rho_overall": rho_all, "p_overall": p_all,
        "ci_overall": [lo_all, hi_all],
        "rho_partial_er": rho_part, "p_partial_er": p_part,
        "rho_er_pos": rho_pos, "p_er_pos": p_pos, "ci_er_pos": [lo_pos, hi_pos],
        "rho_er_neg": rho_neg, "p_er_neg": p_neg, "ci_er_neg": [lo_neg, hi_neg],
    }
    with open(OUT / "fdps_methylation_er_adjusted.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)

    df_out = df.copy()
    df_out["ER_label"] = np.where(df_out["ER"] == 1, "ER+", "ER-")
    df_out.to_csv(OUT / "Supplementary_Table_S7_fdps_methylation_er.csv",
                  index=False, encoding="utf-8-sig")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    panels = [
        ("All (n=%d)" % len(df), df["expr"], df["beta"], rho_all, p_all, "#1f77b4"),
        ("ER+ (n=%d)" % (df["ER"] == 1).sum(), df.loc[df["ER"] == 1, "expr"],
         df.loc[df["ER"] == 1, "beta"], rho_pos, p_pos, "#d62728"),
        ("ER- (n=%d)" % (df["ER"] == 0).sum(), df.loc[df["ER"] == 0, "expr"],
         df.loc[df["ER"] == 0, "beta"], rho_neg, p_neg, "#2ca02c"),
    ]
    for ax, (title, x, y, r, p, c) in zip(axes, panels):
        ax.scatter(x, y, s=10, alpha=0.45, color=c)
        b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, np.polyval(b, xs), color="black", lw=1.5)
        ax.set_xlabel("FDPS expression (log1p)")
        ax.set_ylabel("FDPS methylation (HM450 beta)")
        ax.set_title(f"{title}: rho={r:.3f}, P={p:.2e}", fontsize=10)
        ax.grid(alpha=0.3)
    fig.suptitle("TCGA-BRCA FDPS methylation vs expression, overall and by ER status",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig12b_fdps_methylation_er.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_patch({"er_adjusted": summary}, "methyl_er")
    print("DONE ->", FIG / "fig12b_fdps_methylation_er.png")


if __name__ == "__main__":
    main()
