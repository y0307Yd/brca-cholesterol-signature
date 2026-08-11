# -*- coding: utf-8 -*-
"""Bonus analyses: immune-checkpoint genes by subtype and LASSO stability."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")

RNG = np.random.default_rng(20260804)


def load_tcga():
    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8")
             if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")          # genes x samples
    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    return genes, X, sub


def checkpoint_by_subtype():
    genes, X, sub = load_tcga()
    idx = {g: i for i, g in enumerate(genes)}
    check = ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "IDO1", "CD8A"]
    rows = []
    for g in check:
        x = X[idx[g]]
        med = {}
        groups = []
        for k in range(1, 5):
            mask = (sub["subtype"] == k).to_numpy()
            med[k] = float(np.nanmedian(x[mask]))
            groups.append(x[mask])
        kw = stats.kruskal(*groups)
        rows.append({
            "gene": g,
            "median_C1": med[1], "median_C2": med[2],
            "median_C3": med[3], "median_C4": med[4],
            "kruskal_P": kw.pvalue,
            "highest_subtype": max(med, key=med.get),
            "lowest_subtype": min(med, key=med.get),
        })
    res = pd.DataFrame(rows)
    res.to_csv(WORK / "bonus_checkpoint_by_subtype.csv", index=False)
    print(res.to_string(index=False))
    return res


def lasso_stability(n_boot=300, C=0.316):
    genes, X, sub = load_tcga()
    pool = [g for g in pd.read_csv(OUT / "ml_candidate_pool.csv")["gene"]]
    idx = {g: i for i, g in enumerate(genes)}
    Xp = np.vstack([X[idx[g]] for g in pool]).T          # samples x genes
    Xz = (Xp - Xp.mean(axis=0)) / Xp.std(axis=0, ddof=1)
    y = sub["ER"].to_numpy()
    ok = ~np.isnan(y)
    Xz, y = Xz[ok], y[ok].astype(int)

    counts = np.zeros(len(pool))
    for b in range(n_boot):
        ridx = RNG.integers(0, len(y), size=len(y))
        Xb, yb = Xz[ridx], y[ridx]
        if len(np.unique(yb)) < 2:
            continue
        clf = LogisticRegression(penalty="l1", C=C, solver="liblinear",
                                 max_iter=5000, random_state=0)
        clf.fit(Xb, yb)
        counts += np.abs(clf.coef_[0]) > 1e-8
    freq = pd.DataFrame({"gene": pool, "selection_frequency": counts / n_boot})
    freq = freq.sort_values("selection_frequency", ascending=False)
    freq.to_csv(WORK / "bonus_lasso_stability.csv", index=False)
    hubs = pd.read_csv(OUT / "hub_genes.csv")["gene"].tolist()
    hub_freq = freq[freq["gene"].isin(hubs)]
    print(f"top 15 by selection frequency (n_boot={n_boot}, C={C}):")
    print(freq.head(15).to_string(index=False))
    print("\nhub-gene frequencies:")
    print(hub_freq.sort_values("selection_frequency", ascending=False).to_string(index=False))
    print("\nmean hub frequency:", hub_freq["selection_frequency"].mean().round(3))
    return freq


if __name__ == "__main__":
    print("===== immune checkpoint genes by subtype =====")
    checkpoint_by_subtype()
    print("\n===== LASSO bootstrap stability =====")
    lasso_stability()
