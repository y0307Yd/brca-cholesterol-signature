# -*- coding: utf-8 -*-
"""External validation of TCGA-defined cholesterol subtypes in GEO cohorts.

Applies the TCGA nearest-centroid mapping (hub-gene within-cohort z-scores)
to GSE21653 (GPL570) and GSE7390 (U133A), then compares ER composition and
RFS/DFS survival across mapped subtypes.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")
OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature\figures")
sys.path.insert(0, str(WORK))

from finish_geo_gse21653 import collapse_probe, load_gpl_annot, parse_series_matrix
from online_utils import int_transform

HUB = ["ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
       "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]


def load_tcga():
    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")
    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    idx = {g: i for i, g in enumerate(genes)}
    Xh = np.vstack([X[idx[g]] for g in HUB])
    Zt = ((Xh - Xh.mean(axis=1)[:, None]) / Xh.std(axis=1, ddof=1)[:, None]).T
    lab = sub["subtype"].to_numpy()
    cent = np.vstack([Zt[lab == k].mean(axis=0) for k in [1, 2, 3, 4]])
    return cent


def map_subtypes(expr_hub, cent):
    Zh = ((expr_hub - expr_hub.mean(axis=1)[:, None]) /
          expr_hub.std(axis=1, ddof=1)[:, None]).T
    from scipy.spatial.distance import cdist
    dist = cdist(Zh, cent, metric="euclidean")
    return np.argmin(dist, axis=1) + 1, Zh


def gse21653(cent):
    dl = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE21653")
    clin, expr = parse_series_matrix(dl / "GSE21653_series_matrix.txt.gz")
    mapping = load_gpl_annot(dl / "GPL570.annot.gz")
    hub_probes = {g: [p for p in mapping.get(g, []) if p in expr.columns]
                  for g in HUB}
    hub_expr = pd.DataFrame({g: collapse_probe(expr, ps)
                             for g, ps in hub_probes.items() if ps})
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
                                                        "dfs time", "rfs time"), errors="coerce"))
    meta = pd.DataFrame({"ER": er.values, "event": evt.values,
                         "time": time.values}, index=clin.index)
    meta = meta[meta["time"].notna() & meta["event"].notna()]
    meta["subtype"], _ = map_subtypes(hub_expr.reindex(index=meta.index).T.to_numpy(), cent)
    meta["cohort"] = "GSE21653"
    return meta


def gse7390(cent):
    from finish_geo_gse7390 import parse_sdrf
    matrix = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\E-GEOD-7390-processed-data-1631276164.txt")
    sdrf = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\data\geo\GSE7390_sdrf.txt")
    cols = pd.read_csv(matrix, sep="\t", header=None, nrows=1).iloc[0].tolist()
    expr = pd.read_csv(matrix, sep="\t", header=None, skiprows=2)
    expr.columns = cols
    expr = expr.set_index(cols[0]).apply(pd.to_numeric, errors="coerce")
    clin = parse_sdrf(sdrf)
    common = [s for s in expr.columns if s in clin.index]
    clin = clin.loc[common]
    expr = expr[common]
    for c in ["er", "t_rfs", "e_rfs"]:
        clin[c] = pd.to_numeric(clin[c], errors="coerce")
    clin["er"] = np.where(clin["er"] == 1, 1, np.where(clin["er"] == 0, 0, np.nan))
    pmap = pd.read_csv(OUT / "hgu133a_probe_map.csv")
    by_sym = pmap.groupby("SYMBOL")["PROBEID"].apply(list).to_dict()
    hub_probes = {g: [p for p in by_sym.get(g, []) if p in expr.index]
                  for g in HUB}
    def collapse_probes_rowwise(probe_rows):
        v = probe_rows.var(axis=1)
        return probe_rows.loc[v.idxmax()]

    hub_expr = pd.DataFrame({g: collapse_probes_rowwise(expr.loc[ps])
                             for g, ps in hub_probes.items() if ps})
    meta = pd.DataFrame({"ER": clin["er"], "event": clin["e_rfs"], "time": clin["t_rfs"]},
                        index=clin.index)
    meta = meta[meta["time"].notna() & meta["event"].notna()]
    meta["subtype"], _ = map_subtypes(hub_expr.reindex(index=meta.index).T.to_numpy(), cent)
    meta["cohort"] = "GSE7390"
    return meta


def summarize(meta):
    from scipy import stats
    from lifelines.statistics import multivariate_logrank_test, logrank_test
    from online_utils import coxph
    rows = []
    for k in [1, 2, 3, 4]:
        m = meta[meta["subtype"] == k]
        rows.append({
            "cohort": meta["cohort"].iloc[0], "subtype": k, "n": len(m),
            "ER_pos_frac": m["ER"].mean() if m["ER"].notna().any() else np.nan,
            "events": int(m["event"].sum()),
        })
    res = pd.DataFrame(rows)
    glr = multivariate_logrank_test(meta["time"].values, meta["subtype"].values,
                                    meta["event"].values)
    res.attrs["global_logrank_p"] = float(glr.p_value)
    pair = {}
    for a in [1, 2, 3, 4]:
        for b in range(a + 1, 5):
            lr = logrank_test(meta.loc[meta["subtype"] == a, "time"],
                              meta.loc[meta["subtype"] == b, "time"],
                              meta.loc[meta["subtype"] == a, "event"],
                              meta.loc[meta["subtype"] == b, "event"])
            pair[f"C{a}_vs_C{b}"] = float(lr.p_value)
    res.attrs["pairwise_p"] = pair
    er_med = meta.groupby("subtype")["ER"].mean()
    res.attrs["er_kruskal_p"] = float(stats.kruskal(
        *[meta.loc[meta["subtype"] == k, "ER"].dropna() for k in [1, 2, 3, 4]]).pvalue)
    return res


def main():
    cent = load_tcga()
    m1 = gse21653(cent)
    m2 = gse7390(cent)
    both = pd.concat([m1, m2])
    both.to_csv(OUT / "geo_subtypes_mapped.csv", index=True, encoding="utf-8-sig")
    r1 = summarize(m1)
    r2 = summarize(m2)
    r1.to_csv(OUT / "geo_gse21653_subtype_stats.csv", index=False)
    r2.to_csv(OUT / "geo_gse7390_subtype_stats.csv", index=False)
    print("=== GSE21653 ===")
    print(r1.to_string(index=False))
    print("global log-rank P:", r1.attrs["global_logrank_p"])
    print("pairwise:", {k: round(v, 4) for k, v in r1.attrs["pairwise_p"].items()})
    print("ER Kruskal P:", r1.attrs["er_kruskal_p"])
    print("\n=== GSE7390 ===")
    print(r2.to_string(index=False))
    print("global log-rank P:", r2.attrs["global_logrank_p"])
    print("pairwise:", {k: round(v, 4) for k, v in r2.attrs["pairwise_p"].items()})
    print("ER Kruskal P:", r2.attrs["er_kruskal_p"])

    # combined figure: ER composition + KM per cohort
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    colors = {1: "#d62728", 2: "#1f77b4", 3: "#2ca02c", 4: "#ff7f0e"}
    for row, (m, r, title) in enumerate([(m1, r1, "GSE21653"), (m2, r2, "GSE7390")]):
        ax = axes[row, 0]
        frac = [m.loc[m["subtype"] == k, "ER"].mean() if (m["subtype"] == k).any() else 0 for k in [1, 2, 3, 4]]
        ax.bar(["C1", "C2", "C3", "C4"], frac, color=[colors[k] for k in [1, 2, 3, 4]])
        ax.set_ylabel("Fraction ER+")
        ax.set_title(f"{title} ER composition")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3, axis="y")
        ax = axes[row, 1]
        for k in [1, 2, 3, 4]:
            kmf = KaplanMeierFitter()
            mm = m[m["subtype"] == k]
            if len(mm):
                kmf.fit(mm["time"], mm["event"])
                ax.step(kmf.timeline, kmf.survival_function_["KM_estimate"],
                        where="post", label=f"C{k} (n={len(mm)})", color=colors[k])
        ax.set_xlabel("Months")
        ax.set_ylabel("RFS/DFS probability")
        ax.set_title(f"{title} KM by subtype (log-rank P={r.attrs['global_logrank_p']:.3g})")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.3)
        ax = axes[row, 2]
        ax.axis("off")
        tbl = r.copy()
        tbl["ER_pos_frac"] = tbl["ER_pos_frac"].round(3)
        txt = tbl.to_string(index=False)
        pair = "\n".join(f"{k}: {v:.3g}" for k, v in list(r.attrs["pairwise_p"].items())[:3])
        ax.text(0.02, 0.95, f"{title}\n{txt}\n\nGlobal log-rank P: "
                f"{r.attrs['global_logrank_p']:.3g}\nPairwise P:\n{pair}",
                va="top", fontsize=9, family="monospace")
    fig.suptitle("External validation of cholesterol subtypes in GEO cohorts", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig_subtype_geo_validation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved ->", FIG / "fig_subtype_geo_validation.png")


if __name__ == "__main__":
    main()
