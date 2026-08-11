# -*- coding: utf-8 -*-
"""Immunotherapy-response association in GSE91061 (melanoma anti-PD-1).

Pre-treatment rlog expression of the 11 hub genes, FDPS and immune
features is compared between responders (CR/PR) and non-responders (PD).
Cross-cancer exploratory validation of the cholesterol-signature immune axis.
"""
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = OUT / "figures"
MATRIX = (r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE91061"
          r"\GSE91061_series_matrix.txt.gz")
RLD = (r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE91061"
       r"\GSE91061_rld.csv.gz")

ENTREZ = {
    "ABCG1": 9619, "DHCR24": 1718, "DHCR7": 1717, "FDXR": 2232,
    "G6PD": 2539, "HMGCS2": 3158, "HSD17B7": 51478, "LIMA1": 51474,
    "NSDHL": 50814, "PRKAA1": 5562, "VLDLR": 7436, "FDPS": 2224,
    "GZMA": 3001, "PRF1": 5551, "CD8A": 925, "CD274": 29126,
    "IFNG": 3458,
}


def parse_clinical():
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")]
    title = None
    char_rows = []
    for r in rows:
        if r and r[0] == "!Sample_title":
            title = [x.strip('"') for x in r[1:] if x]
        if r and r[0].startswith("!Sample_characteristics_ch1"):
            char_rows.append([x.strip('"') for x in r[1:] if x])
    clin = pd.DataFrame({"title": title})
    for cr in char_rows:
        keys = set()
        for v in cr:
            if ":" in v:
                keys.add(v.split(":", 1)[0].strip())
        if not keys:
            continue
        key = sorted(keys)[0]
        vals = []
        for v in cr:
            if ":" in v and v.split(":", 1)[0].strip() == key:
                vals.append(v.split(":", 1)[1].strip())
            else:
                vals.append(np.nan)
        if "visit" in key.lower():
            clin["visit"] = vals
        elif "response" in key.lower():
            clin["response"] = vals
        elif "tissue" in key.lower():
            clin["tissue"] = vals
    return clin


def read_rld():
    df = pd.read_csv(RLD, index_col=0)
    df.index = df.index.astype(str)
    return df


def main():
    clin = parse_clinical()
    print("samples:", len(clin))
    print(clin["response"].value_counts(dropna=False).to_dict())
    pre = clin[clin["visit"] == "Pre"].copy()
    print("Pre samples:", len(pre))

    rld = read_rld()
    print("rld:", rld.shape)
    cols = [c for c in pre["title"] if c in rld.columns]
    print("matched Pre samples:", len(cols))
    pre = pre[pre["title"].isin(cols)].copy()
    X = rld[pre["title"].tolist()].astype(float)

    features = {}
    for sym, eid in ENTREZ.items():
        key = str(eid)
        if key in X.index:
            features[sym] = X.loc[key].to_numpy()
        else:
            print("missing Entrez", sym, eid)
    features["CYT"] = (features["GZMA"] + features["PRF1"]) / 2.0

    resp = pre["response"].map({"PRCR": 1, "PD": 0, "SD": np.nan,
                                "UNK": np.nan})
    rows = []
    for f, v in features.items():
        r = resp.to_numpy(dtype=float)
        ok = ~np.isnan(r) & ~np.isnan(v)
        if ok.sum() < 10:
            continue
        r_ok, v_ok = r[ok].astype(int), v[ok]
        mw = stats.mannwhitneyu(v_ok[r_ok == 1], v_ok[r_ok == 0])
        rows.append({
            "feature": f,
            "n_responder": int((r_ok == 1).sum()),
            "n_nonresponder": int((r_ok == 0).sum()),
            "mean_responder": float(v_ok[r_ok == 1].mean()),
            "mean_nonresponder": float(v_ok[r_ok == 0].mean()),
            "mannwhitney_P": mw.pvalue,
            "direction": ("up in responders" if v_ok[r_ok == 1].mean() >
                           v_ok[r_ok == 0].mean() else "down in responders"),
        })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "gse91061_immunotherapy_response.csv", index=False,
               encoding="utf-8-sig")
    res.to_csv(OUT / "Supplementary_Table_S20_immunotherapy_response.csv",
               index=False, encoding="utf-8-sig")
    print(res.round(4).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_feats = ["FDPS", "CD274", "CD8A", "CYT", "IFNG", "DHCR7",
                  "DHCR24", "ABCG1", "VLDLR"]
    fig, axes = plt.subplots(3, 3, figsize=(11, 9))
    for ax, f in zip(axes.ravel(), plot_feats):
        if f not in features:
            ax.axis("off")
            continue
        v = features[f]
        r = resp.to_numpy(dtype=float)
        ok = ~np.isnan(r) & ~np.isnan(v)
        d0 = v[ok & (r == 0)]
        d1 = v[ok & (r == 1)]
        ax.boxplot([d1, d0], tick_labels=["R", "NR"], widths=0.5)
        ax.scatter([1] * len(d1), d1, s=6, alpha=0.4, color="#2ca02c")
        ax.scatter([2] * len(d0), d0, s=6, alpha=0.4, color="#d62728")
        p = stats.mannwhitneyu(d1, d0).pvalue
        ax.set_title(f"{f} (P={p:.2g})")
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle("GSE91061 pre-treatment expression vs anti-PD-1 response "
                 "(melanoma; R=CR/PR, NR=PD)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "supp_fig10_immunotherapy_response.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig10_immunotherapy_response.png")


if __name__ == "__main__":
    main()
