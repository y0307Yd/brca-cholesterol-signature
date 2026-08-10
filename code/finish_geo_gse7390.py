# -*- coding: utf-8 -*-
"""Third independent validation cohort GSE7390 (Desmedt et al. 2007).

Uses the EBI ArrayExpress processed matrix (E-GEOD-7390, Affymetrix U133A,
198 node-negative patients), SDRF clinical annotations (ER, RFS) and the
hgu133a.db probe mapping; applies the TCGA-trained 11-hub-gene ER classifier.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import OUT, FIG, auc_ci, coxph, fmt_p, int_transform, \
    logrank, save_patch

MATRIX = Path(r"data\geo\E-GEOD-7390-processed-data-1631276164.txt")
SDRF = Path(r"data\geo\GSE7390_sdrf.txt")
PROBE_MAP = OUT / "hgu133a_probe_map.csv"
HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]


def parse_sdrf(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    rows = {}
    for _, r in df.iterrows():
        desc = str(r.get("Description", ""))
        kv = {}
        for part in desc.split(";"):
            if ":" in part:
                k, v = part.split(":", 1)
                kv[k.strip()] = v.strip()
        rows[str(r["Scan Name"]).strip()] = {
            "er": kv.get("er", np.nan),
            "t_rfs": kv.get("t.rfs", np.nan),
            "e_rfs": kv.get("e.rfs", np.nan),
            "age": kv.get("age", np.nan),
            "grade": kv.get("grade", np.nan),
            "node": kv.get("node", np.nan),
        }
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    cols = pd.read_csv(MATRIX, sep="\t", header=None, nrows=1).iloc[0].tolist()
    expr = pd.read_csv(MATRIX, sep="\t", header=None, skiprows=2)
    expr.columns = cols
    expr = expr.set_index(cols[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    print("matrix:", expr.shape)
    clin = parse_sdrf(SDRF)
    common = [s for s in expr.columns if s in clin.index]
    print("samples with clinical:", len(common))
    clin = clin.loc[common]
    expr = expr[common]

    for c in ["er", "t_rfs", "e_rfs", "age", "grade", "node"]:
        clin[c] = pd.to_numeric(clin[c], errors="coerce")
    clin["er"] = np.where(clin["er"] == 1, 1, np.where(clin["er"] == 0, 0, np.nan))
    meta = clin.dropna(subset=["t_rfs", "e_rfs"])
    print("samples with RFS:", len(meta), "events:", int(meta["e_rfs"].sum()))
    print("ER available:", int(meta["er"].notna().sum()))
    if len(meta) < 100 or meta["e_rfs"].sum() < 20:
        raise SystemExit("clinical parsing too sparse")

    pmap = pd.read_csv(PROBE_MAP)
    by_sym = pmap.groupby("SYMBOL")["PROBEID"].apply(list).to_dict()
    hub_probes = {g: [p for p in by_sym.get(g, []) if p in expr.index]
                  for g in HUB}
    n_hub = sum(bool(v) for v in hub_probes.values())
    print(f"hub genes with U133A probes: {n_hub}/11")
    if n_hub < 8:
        raise SystemExit("too few hub genes mapped")

    def collapse(row_df):
        v = row_df.var(axis=1)
        return row_df.loc[v.idxmax()]

    hub_expr = pd.DataFrame({g: collapse(expr.loc[ps])
                             for g, ps in hub_probes.items() if ps})
    hub_expr = hub_expr.reindex(index=meta.index)

    fdps_probes = [p for p in by_sym.get("FDPS", []) if p in expr.index]
    if not fdps_probes:
        raise SystemExit("no FDPS probe")
    fdps = collapse(expr.loc[fdps_probes]).reindex(index=meta.index)
    meta["FDPS_int"] = int_transform(fdps.values)

    hub_coef = pd.read_csv(OUT / "hub_genes.csv")
    coef = dict(zip(hub_coef["gene"], hub_coef["coef"]))
    hub_int = pd.DataFrame({g: int_transform(hub_expr[g].values)
                            for g in hub_expr.columns}, index=hub_expr.index)
    lp = sum(hub_int[g] * coef[g] for g in hub_int.columns)

    tcga_genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8")
                  if l.strip()]
    Xtcga = np.load(OUT / "tcga_Xlog.npy")
    tt = pd.read_csv(OUT / "tcga_traits.csv")
    gene_idx = {g: i for i, g in enumerate(tcga_genes)}
    ztcga = {g: (Xtcga[gene_idx[g]] - Xtcga[gene_idx[g]].mean()) /
             (Xtcga[gene_idx[g]].std(ddof=1) or 1) for g in hub_expr.columns}
    lp_tcga = sum(ztcga[g] * coef[g] for g in hub_expr.columns)
    y_tcga = tt["ER"].to_numpy()
    ok = ~np.isnan(y_tcga)
    from sklearn.linear_model import LogisticRegression
    cal = LogisticRegression(max_iter=2000)
    cal.fit(lp_tcga[ok].reshape(-1, 1), y_tcga[ok].astype(int))
    logit = float(cal.intercept_[0]) + float(cal.coef_[0, 0]) * lp
    meta["ER_pred"] = 1.0 / (1.0 + np.exp(-np.clip(logit, -30, 30)))

    y_er = meta["er"].dropna()
    p_er_v = meta.loc[y_er.index, "ER_pred"]
    auc, ci = auc_ci(y_er.astype(int).values, p_er_v.values)
    brier = float(np.mean((p_er_v - y_er.to_numpy()) ** 2))
    er_yes = meta.loc[meta["er"] == 1, "FDPS_int"]
    er_no = meta.loc[meta["er"] == 0, "FDPS_int"]
    p_er_fdps = float(stats.mannwhitneyu(er_yes, er_no).pvalue) \
        if len(er_yes) > 5 and len(er_no) > 5 else np.nan
    d = float((er_yes.mean() - er_no.mean()) / np.sqrt(
        (er_yes.var(ddof=1) + er_no.var(ddof=1)) / 2)) \
        if len(er_yes) > 5 and len(er_no) > 5 else np.nan
    cox_univ = coxph(meta["t_rfs"].to_numpy(), meta["e_rfs"].to_numpy(),
                     meta[["FDPS_int"]].to_numpy(), ["FDPS"])
    med = np.median(meta["FDPS_int"])
    grp = (meta["FDPS_int"] > med).astype(int)
    lr_fdps = logrank(meta["t_rfs"].to_numpy(), meta["e_rfs"].to_numpy(),
                      grp.to_numpy())

    rows = {
        "cohort": "GSE7390", "n_survival": int(len(meta)),
        "n_events": int(meta["e_rfs"].sum()), "n_er": int(meta["er"].notna().sum()),
        "er_pos": int((meta["er"] == 1).sum()), "er_neg": int((meta["er"] == 0).sum()),
        "platform": "GPL96", "hub_genes_mapped": n_hub,
        "er_classifier_auc": auc, "er_classifier_auc_ci_low": ci[0],
        "er_classifier_auc_ci_high": ci[1], "er_classifier_brier": brier,
        "fdps_er_mwu_p": p_er_fdps, "fdps_er_cohens_d": d,
        "fdps_hr_univ": cox_univ["FDPS"]["HR"],
        "fdps_p_univ": cox_univ["FDPS"]["p"],
        "fdps_logrank_p": lr_fdps["p"],
    }
    res = pd.DataFrame([rows])
    res.to_csv(OUT / "geo_gse7390_validation.csv", index=False, encoding="utf-8-sig")
    meta.to_csv(OUT / "geo_gse7390_patient_data.csv", encoding="utf-8-sig")
    print(res.T.to_string())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    from lifelines import KaplanMeierFitter
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    fpr, tpr, _ = roc_curve(y_er.astype(int).to_numpy(), p_er_v.to_numpy())
    axes[0].plot(fpr, tpr, lw=2, color="#1f77b4")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].set_title(f"ER classifier transfer (AUC={auc:.3f}, "
                      f"95%CI {ci[0]:.3f}-{ci[1]:.3f})")
    axes[0].set_xlabel("1 - Specificity"); axes[0].set_ylabel("Sensitivity")
    axes[0].grid(alpha=0.3)
    kmf = KaplanMeierFitter()
    for name, g in [("FDPS low", grp == 0), ("FDPS high", grp == 1)]:
        kmf.fit(meta.loc[g, "t_rfs"], meta.loc[g, "e_rfs"])
        axes[1].step(kmf.timeline, kmf.survival_function_["KM_estimate"],
                     where="post", label=name, lw=1.8)
    axes[1].set_title(f"FDPS expression (log-rank P={fmt_p(lr_fdps['p'])})")
    axes[1].set_xlabel("Days"); axes[1].set_ylabel("RFS probability")
    axes[1].legend(frameon=False); axes[1].grid(alpha=0.3)
    fig.suptitle(f"GSE7390 independent validation (n={len(meta)})")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig10b_geo_gse7390_validation.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)

    save_patch(rows, "gse7390")
    print("DONE ->", OUT / "geo_gse7390_validation.csv")


if __name__ == "__main__":
    main()
