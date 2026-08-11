# -*- coding: utf-8 -*-
"""TCGA-BRCA tumour vs GTEx normal breast for cholesterol genes (Xena)."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = OUT / "figures"
DATA = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\xena")

ENSEMBL = {
    "ABCG1": "ENSG00000160179", "DHCR24": "ENSG00000116133",
    "DHCR7": "ENSG00000172893", "FDXR": "ENSG00000173040",
    "G6PD": "ENSG00000160211", "HMGCS2": "ENSG00000134240",
    "HSD17B7": "ENSG00000032113", "LIMA1": "ENSG00000050405",
    "NSDHL": "ENSG00000147383", "PRKAA1": "ENSG00000132356",
    "VLDLR": "ENSG00000147852", "FDPS": "ENSG00000160752",
}


def main():
    pheno = pd.read_csv(DATA / "TcgaTargetGTEX_phenotype.txt.gz", sep="\t",
                        encoding="latin-1")
    brca = pheno[(pheno["_study"] == "TCGA") &
                 (pheno["_primary_site"] == "Breast") &
                 (pheno["_sample_type"] == "Primary Tumor")]
    gtex = pheno[(pheno["_study"] == "GTEX") &
                 (pheno["_primary_site"] == "Breast")]
    print("TCGA-BRCA tumour:", len(brca), "| GTEx breast normal:", len(gtex))
    brca_ids = set(brca["sample"])
    gtex_ids = set(gtex["sample"])

    want = {v: k for k, v in ENSEMBL.items()}
    rows = {}
    header = None
    with gzip.open(DATA / "TcgaTargetGtex_rsem_gene_tpm.gz", "rt",
                   encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            gid = parts[0].split(".")[0]
            if gid in want:
                vals = np.asarray(pd.to_numeric(parts[1:], errors="coerce"))
                rows[want[gid]] = vals
    print("extracted genes:", sorted(rows))
    X = pd.DataFrame(rows, index=header[1:]).astype(float)
    X = np.log2(X + 1.0)
    brca_ids &= set(X.index)
    gtex_ids &= set(X.index)
    print("matched TCGA:", len(brca_ids), "| GTEx:", len(gtex_ids))

    res = []
    for g in sorted(rows):
        t = X.loc[list(brca_ids), g].to_numpy()
        n = X.loc[list(gtex_ids), g].to_numpy()
        tt, nn = t[~np.isnan(t)], n[~np.isnan(n)]
        mw = stats.mannwhitneyu(tt, nn)
        res.append({
            "gene": g, "n_tumor": len(t), "n_normal": len(n),
            "n_tumor_measured": int(len(tt)), "n_normal_measured": int(len(nn)),
            "mean_log2tpm_tumor": float(np.nanmean(t)),
            "mean_log2tpm_normal": float(np.nanmean(n)),
            "mean_diff_tumor_minus_normal": float(np.nanmean(t) - np.nanmean(n)),
            "mannwhitney_P": mw.pvalue,
            "direction": ("up in tumour" if np.nanmean(t) > np.nanmean(n)
                          else "down in tumour"),
        })
    df = pd.DataFrame(res)
    df.to_csv(OUT / "xena_tcga_gtex_normal_tumor_cholesterol.csv",
              index=False, encoding="utf-8-sig")
    df.to_csv(OUT / "Supplementary_Table_S21_tcga_gtex_normal_tumor.csv",
              index=False, encoding="utf-8-sig")
    print(df.round(4).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    gs = df["gene"].tolist()
    diff = df["mean_diff_tumor_minus_normal"].to_numpy()
    order = np.argsort(diff)
    fig, ax = plt.subplots(figsize=(9, 5))
    cols = ["#c00000" if g in ("ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD",
                               "HMGCS2", "HSD17B7", "LIMA1", "NSDHL",
                               "PRKAA1", "VLDLR") else "#1f77b4" for g in gs]
    ax.barh([gs[i] for i in order], diff[order],
            color=[cols[i] for i in order])
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Mean log2(TPM+1) difference (TCGA-BRCA tumour - GTEx normal breast)")
    ax.set_title("Cholesterol genes in tumour vs normal breast (Xena TCGA+GTEx)")
    ax.grid(alpha=0.3, axis="x")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#c00000", label="Hub gene"),
                       Patch(color="#1f77b4", label="FDPS")],
              frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "supp_fig11_tcga_gtex_normal_tumor.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig11_tcga_gtex_normal_tumor.png")


if __name__ == "__main__":
    main()
