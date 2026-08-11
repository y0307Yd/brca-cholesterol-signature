# -*- coding: utf-8 -*-
"""Normal vs tumour expression of cholesterol genes in GSE15852 (43 pairs)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")
OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = OUT / "figures"
sys.path.insert(0, str(WORK))
from finish_geo_gse21653 import collapse_probe, parse_series_matrix

MATRIX = (r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE15852"
          r"\GSE15852_series_matrix.txt.gz")
PROBE_MAP = OUT / "hgu133a_probe_map.csv"

HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]


def main():
    clin, expr = parse_series_matrix(MATRIX)
    is_normal = clin["histopathological exam"].str.lower().str.contains("normal")
    normal_ids = clin.index[is_normal].tolist()
    tumor_ids = clin.index[~is_normal].tolist()
    print("normal:", len(normal_ids), "tumour:", len(tumor_ids))
    # paired by order of appearance
    pairs = list(zip(normal_ids, tumor_ids))
    print("pairs:", len(pairs))

    pmap = pd.read_csv(PROBE_MAP)
    by_sym = pmap.groupby("SYMBOL")["PROBEID"].apply(list).to_dict()

    def collapse(ps):
        sub = expr[ps]
        v = sub.var(axis=0)
        return sub[v.idxmax()]

    targets = HUB + ["FDPS"]
    rows = []
    normal_mat, tumor_mat = [], []
    for g in targets:
        ps = [p for p in by_sym.get(g, []) if p in expr.columns]
        if not ps:
            print("no probe for", g)
            continue
        x = collapse(ps)
        nv = x[normal_ids].to_numpy(dtype=float)
        tv = x[tumor_ids].to_numpy(dtype=float)
        normal_mat.append(nv)
        tumor_mat.append(tv)
        w = stats.wilcoxon(nv, tv)
        rows.append({
            "gene": g,
            "normal_mean": float(np.mean(nv)),
            "tumor_mean": float(np.mean(tv)),
            "log2FC_tumor_vs_normal": float(np.log2(np.mean(tv) / np.mean(nv))),
            "paired_wilcoxon_P": w.pvalue,
            "direction": "up in tumour" if np.mean(tv) > np.mean(nv) else "down in tumour",
        })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "gse15852_normal_tumor_cholesterol.csv", index=False,
               encoding="utf-8-sig")
    res.to_csv(OUT / "Supplementary_Table_S19_normal_tumor_cholesterol.csv",
               index=False, encoding="utf-8-sig")
    print(res.round(4).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5))
    gs = res["gene"].tolist()
    fc = res["log2FC_tumor_vs_normal"].to_numpy()
    cols = ["#c00000" if g in HUB else "#1f77b4" for g in gs]
    order = np.argsort(fc)
    ax.barh([gs[i] for i in order], fc[order], color=[cols[i] for i in order])
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("log2 fold change (tumour / normal), GSE15852 (43 pairs)")
    ax.set_title("Cholesterol genes: tumour vs normal breast (paired)")
    ax.grid(alpha=0.3, axis="x")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#c00000", label="Hub gene"),
                       Patch(color="#1f77b4", label="FDPS")],
              frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "supp_fig9_normal_tumor_cholesterol.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig9_normal_tumor_cholesterol.png")


if __name__ == "__main__":
    main()
