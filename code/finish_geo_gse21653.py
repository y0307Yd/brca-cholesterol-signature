# -*- coding: utf-8 -*-
"""Remaining task 2: independent GEO validation cohort GSE21653.

Downloads GSE21653 (Affymetrix GPL570, 266 early breast cancers, DFS),
applies the TCGA-trained 11-hub-gene ER classifier (hub_genes.csv) using
per-gene inverse-normal transformation, validates FDPS expression
associations (ER, DFS) and produces figures/Supplementary Table S5.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_geo_gse21653.py
"""
import argparse
import csv
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from online_utils import DATA, FIG, OUT, WORK, auc_ci, coxph, download, fmt_p, \
    int_transform, logrank, read_text_maybe_gzip, save_patch

SERIES = "GSE21653"
DL_DIR = DATA / "geo" / SERIES
URL_MATRIX = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21653/matrix/GSE21653_series_matrix.txt.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE21653&format=file&file=GSE21653_series_matrix.txt.gz",
]
URL_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL570/annot/GPL570.annot.gz"
HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]


def parse_series_matrix(path):
    """Return (sample_frame, expr_frame) from a GEO series matrix."""
    if str(path).endswith(".gz"):
        fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        fh = open(path, "r", encoding="utf-8", errors="replace")
    with fh:
        rows = [r for r in csv.reader(fh, delimiter="\t")]
    header, data = [], []
    in_data = False
    for r in rows:
        if not r:
            continue
        if r[0].startswith("!series_matrix_table_begin"):
            in_data = True
            continue
        if r[0].startswith("!series_matrix_table_end"):
            break
        if in_data:
            if r[0].strip('"').upper().startswith("ID_"):
                continue
            data.append(r)
        elif r[0].startswith("!"):
            header.append(r)

    samples = None
    for r in header:
        key = r[0]
        if "=" in key:
            key = key.split("=", 1)[0].strip()
        if key == "!Sample_geo_accession":
            vals = [r[0].split("=", 1)[1].strip(' "')] if "=" in r[0] else []
            vals += [x.strip('"') for x in r[1:]]
            samples = [x for x in vals if x]
    if not samples:
        raise ValueError("no sample IDs found")

    # clinical fields from characteristics rows
    clin = {s: {} for s in samples}
    for r in header:
        key = r[0]
        if "=" in key:
            key = key.split("=", 1)[0].strip()
        if key.startswith("!Sample_characteristics_ch"):
            vals = [r[0].split("=", 1)[1].strip(' "')] if "=" in r[0] else []
            vals += [x.strip('"') for x in r[1:]]
            for s, v in zip(samples, vals):
                v = v.strip('"')
                if ":" in v:
                    k, val = v.split(":", 1)
                    clin[s][k.strip().lower()] = val.strip()
    clin_df = pd.DataFrame.from_dict(clin, orient="index")

    # expression table
    probes = [r[0].strip('"') for r in data]
    exp = np.array([[float(x) if x.strip('"') not in ("", "NA", "null") else np.nan
                     for x in r[1:]] for r in data], dtype=float)
    expr = pd.DataFrame(exp.T, columns=probes, index=samples)
    return clin_df, expr


def load_gpl_annot(path):
    txt = read_text_maybe_gzip(path)
    lines = [l for l in txt.splitlines()
             if l and not l.startswith(("#", "!", "^"))]
    rows = list(csv.reader(lines, delimiter="\t"))
    header = rows[0]
    sym_idx = [i for i, h in enumerate(header) if h.strip().lower() == "gene symbol"]
    if not sym_idx:
        raise ValueError("Gene Symbol column not found in GPL annotation")
    j = sym_idx[0]
    mapping = {}
    for r in rows[1:]:
        if len(r) <= j or not r[j]:
            continue
        for g in r[j].split("///"):
            g = g.strip()
            if g:
                mapping.setdefault(g, []).append(r[0])
    return mapping


def collapse_probe(expr, probes_for_gene, mode="var"):
    sub = expr[probes_for_gene]
    if mode == "var":
        v = sub.var(axis=0, ddof=1)
        pick = v.idxmax()
    else:
        pick = sub.mean(axis=0).idxmax()
    return sub[pick]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default=None)
    ap.add_argument("--annot", default=None)
    args = ap.parse_args()

    DL_DIR.mkdir(parents=True, exist_ok=True)
    matrix = args.matrix
    if matrix is None:
        matrix = DL_DIR / "GSE21653_series_matrix.txt.gz"
        if not matrix.exists():
            for u in URL_MATRIX:
                try:
                    matrix = Path(download(u, matrix))
                    break
                except Exception as e:
                    print("matrix download failed:", e)
            else:
                raise SystemExit("could not download series matrix; "
                                 "place it at data/geo/GSE21653/")
    annot = args.annot
    if annot is None:
        annot = DL_DIR / "GPL570.annot.gz"
        if not annot.exists():
            annot = Path(download(URL_ANNOT, annot))

    clin, expr = parse_series_matrix(matrix)
    print("samples:", expr.shape[0], "probes:", expr.shape[1])
    print("clinical fields:", list(clin.columns))

    # survival / clinical parsing
    def field(s, *keys):
        for k in keys:
            if k in clin.columns:
                v = clin.loc[s, k]
                if pd.notna(v) and str(v).lower() not in ("na", "null", ""):
                    return str(v)
        return np.nan

    er = clin.index.map(lambda s: pd.to_numeric(field(s, "er ihc", "er"), errors="coerce"))
    evt = clin.index.map(lambda s: pd.to_numeric(field(s, "dfs evt", "event"), errors="coerce"))
    time = clin.index.map(lambda s: pd.to_numeric(field(s, "dfs time (months)",
                                                         "dfs time", "rfs time"),
                                                  errors="coerce"))
    age = clin.index.map(lambda s: pd.to_numeric(field(s, "age at diagnosis", "age"), errors="coerce"))
    grade = clin.index.map(lambda s: pd.to_numeric(field(s, "sbr grade", "grade"), errors="coerce"))
    subtype = clin.index.map(lambda s: field(s, "molecular subtype", "pam50"))
    meta = pd.DataFrame({"ER": er.values, "event": evt.values,
                         "time_months": time.values, "age": age.values,
                         "grade": grade.values, "subtype": subtype.values},
                        index=clin.index)
    meta = meta[meta["time_months"].notna() & meta["event"].notna()]
    print("samples with survival:", len(meta),
          "events:", int(meta["event"].sum()))
    print("ER available:", int(meta["ER"].notna().sum()))
    if len(meta) < 150 or meta["event"].sum() < 20:
        raise SystemExit("survival parsing too sparse; check clinical fields")

    mapping = load_gpl_annot(annot)
    hub_probes = {g: [p for p in mapping.get(g, []) if p in expr.columns]
                  for g in HUB}
    n_hub = sum(bool(v) for v in hub_probes.values())
    print(f"hub genes with GPL570 probes: {n_hub}/11")
    if n_hub < 8:
        raise SystemExit("too few hub genes mapped")
    hub_expr = pd.DataFrame({g: collapse_probe(expr, ps) for g, ps in hub_probes.items() if ps})
    hub_expr = hub_expr.reindex(index=meta.index)

    # FDPS (not a hub, added separately)
    fdps_probes = [p for p in mapping.get("FDPS", []) if p in expr.columns]
    if not fdps_probes:
        raise SystemExit("no FDPS probe on GPL570")
    fdps = collapse_probe(expr, fdps_probes).reindex(index=meta.index)
    fdps_int = int_transform(fdps.values)
    meta["FDPS_int"] = fdps_int
    print("FDPS probe(s):", fdps_probes)

    # --- TCGA classifier transfer -------------------------------------
    hub_coef = pd.read_csv(OUT / "hub_genes.csv")
    coef = dict(zip(hub_coef["gene"], hub_coef["coef"]))
    tcga_genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
    Xtcga = np.load(OUT / "tcga_Xlog.npy")
    tt = pd.read_csv(OUT / "tcga_traits.csv")
    gene_idx = {g: i for i, g in enumerate(tcga_genes)}
    ztcga = {}
    for g in hub_expr.columns:
        x = Xtcga[gene_idx[g]]
        ztcga[g] = (x - x.mean()) / (x.std(ddof=1) if x.std(ddof=1) > 0 else 1)
    lp_tcga = sum(ztcga[g] * coef[g] for g in hub_expr.columns)
    y_tcga = tt["ER"].to_numpy()
    ok = ~np.isnan(y_tcga)
    lp_tcga = lp_tcga[ok]
    y_tcga = y_tcga[ok].astype(int)
    # recalibrate intercept/slope on TCGA (AUC unchanged by monotone transform)
    from sklearn.linear_model import LogisticRegression
    cal = LogisticRegression(max_iter=2000)
    cal.fit(lp_tcga.reshape(-1, 1), y_tcga)
    b0, b1 = float(cal.intercept_[0]), float(cal.coef_[0, 0])

    hub_int = pd.DataFrame({g: int_transform(hub_expr[g].values)
                            for g in hub_expr.columns}, index=hub_expr.index)
    lp_gse = sum(hub_int[g] * coef[g] for g in hub_int.columns)
    logit = b0 + b1 * lp_gse
    p_er = 1.0 / (1.0 + np.exp(-np.clip(logit, -30, 30)))
    meta["ER_pred"] = p_er

    y_er = meta["ER"].dropna()
    p_er_v = p_er.reindex(y_er.index)
    auc, ci = auc_ci(y_er.astype(int).values, p_er_v.values)
    brier = float(np.mean((p_er_v - y_er.to_numpy()) ** 2))
    # calibration intercept/slope via logistic recalibration in GSE21653
    lpl = np.log(p_er_v / (1 - p_er_v)).replace([np.inf, -np.inf], np.nan).dropna()
    yy = y_er.loc[lpl.index]
    if len(lpl) > 50:
        cal2 = LogisticRegression(max_iter=2000).fit(lpl.to_numpy().reshape(-1, 1),
                                                     yy.astype(int).to_numpy())
        cal_int, cal_slope = float(cal2.intercept_[0]), float(cal2.coef_[0, 0])
    else:
        cal_int = cal_slope = np.nan

    # --- FDPS clinical associations ------------------------------------
    m = meta.copy()
    er_yes = m.loc[m["ER"] == 1, "FDPS_int"]
    er_no = m.loc[m["ER"] == 0, "FDPS_int"]
    if len(er_yes) > 5 and len(er_no) > 5:
        p_er_fdps = float(stats.mannwhitneyu(er_yes, er_no).pvalue)
        d = (er_yes.mean() - er_no.mean()) / np.sqrt(
            (er_yes.var(ddof=1) + er_no.var(ddof=1)) / 2)
    else:
        p_er_fdps, d = np.nan, np.nan
    cox_univ = coxph(m["time_months"].to_numpy(), m["event"].to_numpy(),
                     m[["FDPS_int"]].to_numpy(), ["FDPS"])
    adj = m.dropna(subset=["FDPS_int", "ER", "age", "grade"])
    if len(adj) > 50:
        cox_multi = coxph(adj["time_months"].to_numpy(), adj["event"].to_numpy(),
                          adj[["FDPS_int", "ER", "age", "grade"]].to_numpy(),
                          ["FDPS", "ER", "age", "grade"])
    else:
        cox_multi = {}
    med = np.median(m["FDPS_int"])
    grp = (m["FDPS_int"] > med).astype(int)
    lr_fdps = logrank(m["time_months"].to_numpy(), m["event"].to_numpy(), grp.to_numpy())

    tert = pd.qcut(m["ER_pred"], 3, labels=["T1", "T2", "T3"], duplicates="drop")
    lr_tert = logrank(m["time_months"].to_numpy(), m["event"].to_numpy(),
                      (tert == "T3").astype(int).to_numpy())

    # survival tables for KM plots
    from lifelines import KaplanMeierFitter
    km_fdps = []
    for name, g in [("FDPS low", grp == 0), ("FDPS high", grp == 1)]:
        kmf = KaplanMeierFitter()
        kmf.fit(m.loc[g, "time_months"], m.loc[g, "event"])
        km_fdps.append((name, kmf.timeline, kmf.survival_function_["KM_estimate"]))
    km_cls = []
    for name, g in [("ER-classifier T1", tert == "T1"), ("ER-classifier T3", tert == "T3")]:
        kmf = KaplanMeierFitter()
        kmf.fit(m.loc[g, "time_months"], m.loc[g, "event"])
        km_cls.append((name, kmf.timeline, kmf.survival_function_["KM_estimate"]))

    # --- outputs ---------------------------------------------------------
    rows = {
        "cohort": SERIES, "n_survival": int(len(m)),
        "n_events": int(m["event"].sum()), "n_er": int(m["ER"].notna().sum()),
        "er_pos": int((m["ER"] == 1).sum()), "er_neg": int((m["ER"] == 0).sum()),
        "platform": "GPL570",
        "hub_genes_mapped": n_hub,
        "er_classifier_auc": auc, "er_classifier_auc_ci_low": ci[0],
        "er_classifier_auc_ci_high": ci[1], "er_classifier_brier": brier,
        "er_classifier_cal_int": cal_int, "er_classifier_cal_slope": cal_slope,
        "fdps_probe": ";".join(fdps_probes),
        "fdps_er_mwu_p": p_er_fdps, "fdps_er_cohens_d": float(d),
        "fdps_er_pos_mean": float(er_yes.mean()) if len(er_yes) else np.nan,
        "fdps_er_neg_mean": float(er_no.mean()) if len(er_no) else np.nan,
        "fdps_hr_univ": cox_univ["FDPS"]["HR"],
        "fdps_hr_univ_ci_low": cox_univ["FDPS"]["CI95"][0],
        "fdps_hr_univ_ci_high": cox_univ["FDPS"]["CI95"][1],
        "fdps_p_univ": cox_univ["FDPS"]["p"],
        "fdps_hr_multi": cox_multi.get("FDPS", {}).get("HR", np.nan),
        "fdps_hr_multi_ci_low": cox_multi.get("FDPS", {}).get("CI95", [np.nan, np.nan])[0],
        "fdps_hr_multi_ci_high": cox_multi.get("FDPS", {}).get("CI95", [np.nan, np.nan])[1],
        "fdps_p_multi": cox_multi.get("FDPS", {}).get("p", np.nan),
        "fdps_logrank_p": lr_fdps["p"],
        "classifier_t3_vs_t1_logrank_p": lr_tert["p"],
    }
    res = pd.DataFrame([rows])
    res.to_csv(OUT / "geo_gse21653_validation.csv", index=False, encoding="utf-8-sig")
    m.to_csv(OUT / "geo_gse21653_patient_data.csv", encoding="utf-8-sig")
    print(res.T.to_string())

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.6))
    fpr, tpr, _ = roc_curve(y_er.astype(int).to_numpy(), p_er_v.to_numpy())
    ax = axes[0, 0]
    ax.plot(fpr, tpr, color="#1f77b4", lw=2)
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title(f"ER classifier transfer (AUC={auc:.3f}, 95%CI {ci[0]:.3f}-{ci[1]:.3f})")
    ax.grid(alpha=0.3)
    for name, t, s in km_fdps:
        axes[0, 1].step(t, s, where="post", label=name, lw=1.8)
    axes[0, 1].set_xlabel("Months")
    axes[0, 1].set_ylabel("DFS probability")
    axes[0, 1].set_title(f"FDPS expression (log-rank P={fmt_p(lr_fdps['p'])}")
    axes[0, 1].legend(frameon=False)
    axes[0, 1].grid(alpha=0.3)
    for name, t, s in km_cls:
        axes[1, 0].step(t, s, where="post", label=name, lw=1.8)
    axes[1, 0].set_xlabel("Months")
    axes[1, 0].set_ylabel("DFS probability")
    axes[1, 0].set_title(f"ER-classifier tertiles (log-rank P={fmt_p(lr_tert['p'])}")
    axes[1, 0].legend(frameon=False)
    axes[1, 0].grid(alpha=0.3)
    bp = axes[1, 1]
    data_box = [m.loc[m["ER"] == 0, "FDPS_int"].dropna().values,
                m.loc[m["ER"] == 1, "FDPS_int"].dropna().values]
    bp.boxplot(data_box, tick_labels=["ER-", "ER+"], widths=0.55)
    bp.scatter([1] * len(data_box[0]), data_box[0], s=8, alpha=0.35, color="#1f77b4")
    bp.scatter([2] * len(data_box[1]), data_box[1], s=8, alpha=0.35, color="#d62728")
    bp.set_title(f"FDPS by ER status (MWU P={fmt_p(p_er_fdps)})")
    bp.set_ylabel("FDPS expression (INT)")
    bp.grid(alpha=0.3, axis="y")
    fig.suptitle(f"Independent validation: {SERIES} (n={len(m)}, events={int(m['event'].sum())})",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig10_geo_gse21653_validation.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)

    supp = m.copy()
    supp.to_csv(OUT / "Supplementary_Table_S5_geo_gse21653_validation.csv",
                encoding="utf-8-sig")
    patch = dict(rows)
    patch["table_path"] = str(OUT / "Supplementary_Table_S5_geo_gse21653_validation.csv")
    patch["figure_path"] = str(FIG / "fig10_geo_gse21653_validation.png")
    save_patch(patch, "geo")
    print("DONE ->", FIG / "fig10_geo_gse21653_validation.png")


if __name__ == "__main__":
    main()
