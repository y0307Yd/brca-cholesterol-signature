# -*- coding: utf-8 -*-
"""Independent validation cohort GSE20711 (Affymetrix GPL570).

Node-negative, systemically untreated early breast cancers. Applies the
TCGA-trained 11-hub-gene ER classifier with per-gene inverse-normal
transformation and reports AUC/Brier, FDPS-ER associations, and subtype
nearest-centroid mapping.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from online_utils import OUT, FIG, auc_ci, coxph, fmt_p, int_transform, \
    logrank, save_patch

WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")
sys.path.insert(0, str(WORK))

from finish_geo_gse21653 import collapse_probe, load_gpl_annot, parse_series_matrix

SERIES = "GSE20711"
DL_DIR = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo") / SERIES
HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]


def main():
    matrix = DL_DIR / "GSE20711_series_matrix.txt.gz"
    if not matrix.exists() or matrix.stat().st_size < 10_000_000:
        raise SystemExit("series matrix missing/incomplete; run the download "
                         "first (data/geo/GSE20711/)")
    annot = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE21653\GPL570.annot.gz")
    clin, expr = parse_series_matrix(matrix)
    print("samples:", expr.shape[0], "probes:", expr.shape[1])
    print("clinical fields:", list(clin.columns))

    def field(s, *keys):
        for k in keys:
            if k in clin.columns:
                v = clin.loc[s, k]
                if pd.notna(v) and str(v).lower() not in ("na", "null", ""):
                    return str(v)
        return np.nan

    er = clin.index.map(lambda s: pd.to_numeric(
        field(s, "er status", "er ihc", "er", "erihc"), errors="coerce"))
    evt = clin.index.map(lambda s: pd.to_numeric(
        field(s, "e.rfs", "dfs evt", "event", "rfs event", "e_rfs",
              "relapse"), errors="coerce"))
    time_years = clin.index.map(lambda s: pd.to_numeric(
        field(s, "t.rfs", "dfs time (months)", "dfs time", "rfs time",
              "t_rfs", "time to relapse (months)", "follow up (months)"),
        errors="coerce"))
    # GSE20711 reports RFS time in years; convert to months for consistency
    time = np.where(time_years.notna() & (time_years <= 20),
                    time_years * 12.0, time_years)
    meta = pd.DataFrame({"ER": er.values, "event": evt.values,
                         "time_months": np.asarray(time, dtype=float)},
                        index=clin.index)
    meta = meta[meta["time_months"].notna() & meta["event"].notna()]
    print("samples with survival:", len(meta), "events:",
          int(meta["event"].sum()))
    print("ER available:", int(meta["ER"].notna().sum()))

    mapping = load_gpl_annot(annot)
    hub_probes = {g: [p for p in mapping.get(g, []) if p in expr.columns]
                  for g in HUB}
    n_hub = sum(bool(v) for v in hub_probes.values())
    print(f"hub genes with GPL570 probes: {n_hub}/11")
    if n_hub < 8:
        raise SystemExit("too few hub genes mapped")
    hub_expr = pd.DataFrame({g: collapse_probe(expr, ps)
                             for g, ps in hub_probes.items() if ps})
    hub_expr = hub_expr.reindex(index=meta.index)

    fdps_probes = [p for p in mapping.get("FDPS", []) if p in expr.columns]
    if not fdps_probes:
        raise SystemExit("no FDPS probe")
    fdps = collapse_probe(expr, fdps_probes).reindex(index=meta.index)
    meta["FDPS_int"] = int_transform(fdps.values)
    print("FDPS probe(s):", fdps_probes)

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

    y_er = meta["ER"].dropna()
    p_er_v = meta.loc[y_er.index, "ER_pred"]
    auc, ci = auc_ci(y_er.astype(int).values, p_er_v.values)
    brier = float(np.mean((p_er_v - y_er.to_numpy()) ** 2))
    er_yes = meta.loc[meta["ER"] == 1, "FDPS_int"]
    er_no = meta.loc[meta["ER"] == 0, "FDPS_int"]
    from scipy import stats
    p_er_fdps = float(stats.mannwhitneyu(er_yes, er_no).pvalue) \
        if len(er_yes) > 5 and len(er_no) > 5 else np.nan
    d = (er_yes.mean() - er_no.mean()) / np.sqrt(
        (er_yes.var(ddof=1) + er_no.var(ddof=1)) / 2) \
        if len(er_yes) > 5 and len(er_no) > 5 else np.nan
    cox_univ = coxph(meta["time_months"].to_numpy(), meta["event"].to_numpy(),
                     meta[["FDPS_int"]].to_numpy(), ["FDPS"])
    med = np.median(meta["FDPS_int"])
    grp = (meta["FDPS_int"] > med).astype(int)
    lr_fdps = logrank(meta["time_months"].to_numpy(), meta["event"].to_numpy(),
                      grp.to_numpy())

    rows = {
        "cohort": SERIES, "n_survival": int(len(meta)),
        "n_events": int(meta["event"].sum()),
        "n_er": int(meta["ER"].notna().sum()),
        "er_pos": int((meta["ER"] == 1).sum()),
        "er_neg": int((meta["ER"] == 0).sum()),
        "platform": "GPL570",
        "hub_genes_mapped": n_hub,
        "er_classifier_auc": auc,
        "er_classifier_auc_ci_low": ci[0],
        "er_classifier_auc_ci_high": ci[1],
        "er_classifier_brier": brier,
        "fdps_probe": ";".join(fdps_probes),
        "fdps_er_mwu_p": p_er_fdps, "fdps_er_cohens_d": float(d),
        "fdps_hr_univ": cox_univ["FDPS"]["HR"],
        "fdps_hr_univ_ci_low": cox_univ["FDPS"]["CI95"][0],
        "fdps_hr_univ_ci_high": cox_univ["FDPS"]["CI95"][1],
        "fdps_p_univ": cox_univ["FDPS"]["p"],
        "fdps_logrank_p": lr_fdps["p"],
    }
    res = pd.DataFrame([rows])
    res.to_csv(OUT / f"geo_{SERIES.lower()}_validation.csv", index=False,
               encoding="utf-8-sig")
    meta.to_csv(OUT / f"geo_{SERIES.lower()}_patient_data.csv",
                encoding="utf-8-sig")
    meta.to_csv(OUT / f"Supplementary_Table_S16_{SERIES.lower()}_validation.csv",
                encoding="utf-8-sig")
    print(res.T.to_string())

    # subtype nearest-centroid mapping
    sys.path.insert(0, str(WORK))
    from subtype_geo_validation import load_tcga, map_subtypes
    cent = load_tcga()
    sub, _ = map_subtypes(hub_expr.T.to_numpy(), cent)
    meta["subtype"] = sub
    meta.to_csv(OUT / "geo_gse20711_subtypes_mapped.csv", encoding="utf-8-sig")
    print("\nsubtype counts:", meta["subtype"].value_counts().sort_index().to_dict())
    print("ER fraction by subtype:",
          meta.groupby("subtype")["ER"].mean().round(3).to_dict())
    from lifelines.statistics import multivariate_logrank_test
    glr = multivariate_logrank_test(meta["time_months"].values,
                                    meta["subtype"].values,
                                    meta["event"].values)
    print("global log-rank P:", glr.p_value)
    save_patch(dict(rows), "geo20711")
    print("DONE")


if __name__ == "__main__":
    main()
