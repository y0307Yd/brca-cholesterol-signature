# -*- coding: utf-8 -*-
"""Remaining task 1: PAM50 composition of the four hub-gene molecular subtypes.

Downloads PAM50 calls (UCSC Xena S3, cBioPortal, or Xena clinical matrix
fallback), merges with outputs/chol_metab_signature/tcga_subtypes.csv and
writes cross-tabulation, chi-square statistics, a figure and
Supplementary Table S4 plus v13 patch numbers.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_pam50.py
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_pam50.py --data data\\pam50.tsv
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import DATA, FIG, OUT, WORK, download, read_text_maybe_gzip, save_patch

URL_CLINICAL = ("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
                "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix")
URL_PAM50_OBJ = ("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
                 "TCGA.BRCA.sampleMap/PAM50")
URL_CBIO = ("https://www.cbioportal.org/api/studies/"
            "brca_tcga_pan_can_atlas_2018/clinical-data"
            "?attributeId=SUBTYPE&clinicalDataType=SAMPLE&projection=SUMMARY")

PAM50_ORDER = ["Basal", "HER2", "LumA", "LumB", "Normal"]
SUBTYPE_ORDER = ["C1", "C2", "C3", "C4"]

PAM50_ALIAS = {
    "basal": "Basal", "basal-like": "Basal", "basal_like": "Basal",
    "her2": "HER2", "her2-enriched": "HER2", "her2_enriched": "HER2",
    "her2+": "HER2",
    "luma": "LumA", "luminal a": "LumA", "luminala": "LumA",
    "lumb": "LumB", "luminal b": "LumB", "luminalb": "LumB",
    "normal": "Normal", "normal-like": "Normal", "normal_like": "Normal",
}


def norm_pam50(v):
    if pd.isna(v):
        return np.nan
    return PAM50_ALIAS.get(str(v).strip().lower().replace("-", " ").replace("_", " "),
                           str(v).strip())


def norm_sample(s):
    s = str(s).strip()
    if s.startswith("TCGA-"):
        parts = s.split("-")
        if len(parts) >= 4:
            return "-".join(parts[:3]) + "-01"
    return s


def load_pam50(path, source):
    txt = read_text_maybe_gzip(path)
    if source == "cbioportal":
        arr = json.loads(txt)
        rows = {norm_sample(x.get("sampleId")): x.get("value") for x in arr
                if x.get("value")}
        return pd.Series(rows, name="PAM50")

    def pick_pam50_col(cols):
        for pat in ("pam50call", "pam50_call"):
            hits = [c for c in cols if pat in c.lower()]
            if hits:
                return hits[0]
        hits = [c for c in cols if c.lower() == "pam50"]
        if hits:
            return hits[0]
        hits = [c for c in cols if "pam50" in c.lower()]
        if hits:
            return hits[0]
        raise ValueError("PAM50 column not found")

    df = pd.read_csv(io_string(txt), sep="\t", dtype=str)
    first = df.columns[0]
    if first.lower() in ("sample", "sampleid", ""):
        col = pick_pam50_col(df.columns)
        s = pd.Series(df[col].values, index=df[first].map(norm_sample))
    else:
        # clinical matrix fallback: rows are samples, look for PAM50Call_RNAseq
        col = pick_pam50_col(df.columns)
        s = pd.Series(df[col].values, index=df[first].map(norm_sample))
    s = s.dropna()
    return s.map(norm_pam50).dropna()


def io_string(txt):
    import io
    return io.StringIO(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="optional local PAM50 file")
    args = ap.parse_args()

    local = args.data
    source = "local"
    if local is None:
        dl_dir = DATA / "pam50"
        local = dl_dir / "BRCA_clinicalMatrix.tsv"
        if not local.exists():
            print("== downloading PAM50 ==")
            try:
                download(URL_CLINICAL, local)
                source = "clinical"
            except Exception as e:
                print("Xena clinicalMatrix failed:", e)
                try:
                    local = dl_dir / "PAM50.tsv"
                    download(URL_PAM50_OBJ, local)
                    source = "xena"
                except Exception as e2:
                    print("Xena PAM50 object failed:", e2)
                    local = dl_dir / "cbioportal_pam50.json"
                    download(URL_CBIO, local)
                    source = "cbioportal"
    pam = load_pam50(local, source)
    print(f"PAM50 loaded: n={len(pam)}, source={source}")
    print(pam.value_counts().to_string())

    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    sub["sample_norm"] = sub["sample_id"].map(norm_sample)
    merged = sub[["sample_norm", "subtype"]].merge(
        pam.rename("PAM50"), left_on="sample_norm", right_index=True, how="inner")
    merged["subtype"] = "C" + merged["subtype"].astype(int).astype(str)
    print(f"matched: {len(merged)}/{len(sub)} TCGA samples")
    if len(merged) < 100:
        raise SystemExit("too few PAM50 matches")

    ct = pd.crosstab(merged["subtype"], merged["PAM50"])
    ct = ct.reindex(index=[s for s in SUBTYPE_ORDER if s in ct.index],
                    columns=[p for p in PAM50_ORDER if p in ct.columns],
                    fill_value=0)
    chi2, p_chi, dof, exp = stats.chi2_contingency(ct.values)
    n = int(ct.values.sum())
    cramer = float(np.sqrt(chi2 / (n * (min(ct.shape) - 1)))) if n and min(ct.shape) > 1 else np.nan
    pct = ct.div(ct.sum(axis=1), axis=0) * 100

    ct_out = ct.copy()
    ct_out["n"] = ct.sum(axis=1)
    ct_out.to_csv(OUT / "pam50_subtype_crosstab.csv", encoding="utf-8-sig")
    pct.to_csv(OUT / "pam50_subtype_percent.csv", encoding="utf-8-sig")
    print("chi2 =", chi2, "p =", p_chi, "dof =", dof, "Cramer V =", cramer)
    print(ct.to_string())
    print(pct.round(1).to_string())

    # figure: stacked bars + heatmap of -log10 expected-over/under enrichment
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    colors = {"Basal": "#d62728", "HER2": "#ff7f0e", "LumA": "#2ca02c",
              "LumB": "#1f77b4", "Normal": "#7f7f7f"}
    bottom = np.zeros(len(ct.index))
    for j, col in enumerate(ct.columns):
        vals = ct[col].values
        ax[0].bar(ct.index, vals, bottom=bottom, label=col,
                  color=colors.get(col, "#999999"), width=0.62)
        for i, sub in enumerate(ct.index):
            v = pct.loc[sub, col]
            if v >= 5:
                ypos = bottom[i] + vals[i] / 2
                frac = vals[i] / max(ct.sum(axis=1)[sub], 1)
                ax[0].text(sub, ypos, f"{v:.0f}%",
                           ha="center", va="center", fontsize=7.5,
                           color="white" if frac > 0.4 else "black")
        bottom += vals
    ax[0].set_ylabel("Patients")
    ax[0].set_xlabel("Hub-gene molecular subtype")
    ax[0].set_title("PAM50 composition (n=%d)" % n)
    ax[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
    ax[0].tick_params(axis="x", rotation=0)

    im = ax[1].imshow(pct.values, cmap="YlOrRd", aspect="auto")
    ax[1].set_xticks(range(len(ct.columns)), ct.columns, rotation=30, ha="right")
    ax[1].set_yticks(range(len(ct.index)), ct.index)
    for i in range(len(ct.index)):
        for j in range(len(ct.columns)):
            ax[1].text(j, i, f"{pct.values[i, j]:.1f}", ha="center", va="center",
                       fontsize=8)
    ax[1].set_title("Row percentages")
    fig.colorbar(im, ax=ax[1], label="% of subtype", shrink=0.8)
    fig.suptitle("PAM50 calls by hub-gene molecular subtype (TCGA-BRCA)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig9_pam50_subtype_composition.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)

    supp = ct.copy().reset_index().rename(columns={"index": "subtype"})
    supp.to_csv(OUT / "Supplementary_Table_S4_pam50_subtype_composition.csv",
                index=False, encoding="utf-8-sig")

    dom = {}
    for s in ct.index:
        top = pct.loc[s].idxmax()
        dom[s] = {"top": str(top), "pct": float(pct.loc[s, top])}
    patch = {
        "n_matched": int(len(merged)),
        "chi2": float(chi2),
        "p_chi2": float(p_chi),
        "cramer_v": cramer,
        "dominant": dom,
        "table_path": str(OUT / "Supplementary_Table_S4_pam50_subtype_composition.csv"),
        "figure_path": str(FIG / "fig9_pam50_subtype_composition.png"),
    }
    save_patch(patch, "pam50")
    print("DONE ->", FIG / "fig9_pam50_subtype_composition.png")


if __name__ == "__main__":
    main()
