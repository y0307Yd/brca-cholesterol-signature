# -*- coding: utf-8 -*-
"""Remaining task 3: second single-cell atlas GSE161529 (Pal et al. 2021).

Downloads GSE161529_RAW.tar (2.2 GB, 69 10x profiles, ~421k cells) and
`GSE161529_features.tsv.gz`, extracts per-sample 10x matrices, then runs a
scanpy-free preprocessing + marker-based cell-type annotation pipeline
(numpy/scipy/sklearn/matplotlib), localises FDPS and the 11 hub genes, and
writes Supplementary Table S6 + figures + v13 patch numbers.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_scrna_gse161529.py
"""
import argparse
import gzip
import shutil
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

from online_utils import DATA, FIG, OUT, download, save_patch

SERIES = "GSE161529"
DL_DIR = DATA / "scRNA" / SERIES
URL_TAR = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/GSE161529_RAW.tar"
URL_FEAT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/GSE161529_features.tsv.gz"

GENES = ["FDPS", "ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2",
         "HSD17B7", "LIMA1", "NSDHL", "PRKAA1", "VLDLR"]

CELL_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "CD2"],
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB"],
    "CD4_T": ["CD4", "IL7R", "CCR7"],
    "Treg": ["FOXP3", "CTLA4", "IL2RA"],
    "NK": ["NKG7", "KLRD1", "GNLY"],
    "B_cell": ["MS4A1", "CD79A", "CD79B"],
    "Plasma": ["MZB1", "JCHAIN", "TNFRSF17"],
    "Macrophage": ["CD68", "CD163", "CSF1R"],
    "Monocyte": ["LYZ", "FCN1", "S100A8", "S100A9"],
    "DC": ["ITGAX", "CD1C", "CLEC9A", "FLT3"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Fibroblast": ["COL1A1", "COL3A1", "DCN", "PDGFRB"],
    "Epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Myoepithelial": ["ACTA2", "MYLK", "KRT14", "KRT5"],
}
LIGAND_RECEPTOR = [("CD74", "MIF"), ("LGALS9", "HAVCR2"), ("SPP1", "CD44"),
                   ("CCL18", "CCR8"), ("CXCL16", "CXCR6")]


def download_data():
    DL_DIR.mkdir(parents=True, exist_ok=True)
    tar = DL_DIR / "GSE161529_RAW.tar"
    feat = DL_DIR / "GSE161529_features.tsv.gz"
    if not tar.exists() or tar.stat().st_size < 100_000_000:
        print("== downloading GSE161529_RAW.tar (2.2 GB, resumable) ==")
        download(URL_TAR, tar, max_try=5, timeout=300)
    if not feat.exists():
        download(URL_FEAT, feat)
    return tar, feat


def extract_tar(tar_path, out_dir):
    if (out_dir / "_extract_done.flag").exists():
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    print("== extracting (may take several minutes) ==")
    with tarfile.open(tar_path, "r:*") as tf:
        members = tf.getmembers()
        for i, m in enumerate(members):
            if m.isfile():
                try:
                    tf.extract(m, out_dir, filter="data")
                except TypeError:
                    tf.extract(m, out_dir)
            if (i + 1) % 500 == 0:
                print(f"  extracted {i + 1}/{len(members)}")
    (out_dir / "_extract_done.flag").write_text("ok", encoding="utf-8")
    return out_dir


def organize_10x(root):
    """Group matrix/barcodes/features into per-sample dirs; handle flat names."""
    samples = defaultdict(dict)
    for f in root.rglob("*.gz"):
        name = f.name.lower()
        parent = str(f.parent)
        if "matrix.mtx" in name:
            samples[parent]["matrix"] = f
        elif "barcodes.tsv" in name:
            samples[parent]["barcodes"] = f
        elif "features.tsv" in name or "genes.tsv" in name:
            samples[parent]["features"] = f
    ready = [v for v in samples.values()
             if all(k in v for k in ("matrix", "barcodes", "features"))]
    if ready:
        return sorted(ready, key=lambda d: str(d["matrix"]))
    # flat fallback: files carry GSM prefixes
    flat = {}
    for f in root.rglob("*.gz"):
        n = f.name
        gsm = None
        for part in n.split("_"):
            if part.upper().startswith("GSM"):
                gsm = part.upper()
                break
        if gsm is None:
            continue
        flat.setdefault(gsm, {})
        nl = n.lower()
        if "matrix.mtx" in nl:
            flat[gsm]["matrix"] = f
        elif "barcodes.tsv" in nl:
            flat[gsm]["barcodes"] = f
        elif "features.tsv" in nl or "genes.tsv" in nl:
            flat[gsm]["features"] = f
    ready = [v for v in flat.values()
             if all(k in v for k in ("matrix", "barcodes", "features"))]
    return sorted(ready, key=lambda d: str(d["matrix"]))


def load_sample(files, feature_fallback=None):
    import scipy.io
    mtx = scipy.io.mmread(files["matrix"])
    mtx = mtx.tocsr()
    with gzip.open(files["features"], "rt", errors="replace") as f:
        frows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    if frows and len(frows[0]) >= 2:
        symbols = [r[1] if len(r) > 1 else r[0] for r in frows]
    elif feature_fallback is not None:
        symbols = feature_fallback
    else:
        symbols = [str(i) for i in range(mtx.shape[0])]
    with gzip.open(files["barcodes"], "rt", errors="replace") as f:
        barcodes = [l.strip() for l in f if l.strip()]
    if mtx.shape[1] != len(barcodes):
        barcodes = barcodes[:mtx.shape[1]]
    return mtx, np.asarray(symbols, dtype=object), np.asarray(barcodes, dtype=object)


def concat_samples(entries):
    gene_sets = [set(e[1]) for e in entries]
    universe = sorted(set().union(*gene_sets))
    gene_idx = {g: i for i, g in enumerate(universe)}
    cols, rows = [], []
    for k, (mtx, symbols, barcodes) in enumerate(entries):
        idx = np.array([gene_idx.get(s, -1) for s in symbols])
        keep = idx >= 0
        sub = mtx[keep]
        j = idx[keep]
        rows.append(sparse.coo_matrix((sub.data, (np.repeat(np.arange(sub.shape[0]), np.diff(sub.indptr)),
                                                  j[sub.indices])),
                                      shape=(sub.shape[0], len(universe))).tocsr())
        cols.append(np.asarray([f"GSE161529_{k}_{b}" for b in barcodes], dtype=object))
    X = sparse.vstack(rows).tocsr()
    return X, np.asarray(universe, dtype=object), np.concatenate(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", default=None,
                    help="already-extracted dir (optional)")
    ap.add_argument("--max_cells", type=int, default=0,
                    help="subsample for testing (0 = all cells)")
    args = ap.parse_args()

    if args.data_dir:
        ex_dir = Path(args.data_dir)
    elif (DL_DIR / "extracted" / "_extract_done.flag").exists():
        ex_dir = DL_DIR / "extracted"
        print("using existing extraction:", ex_dir)
    else:
        tar, feat = download_data()
        ex_dir = DL_DIR / "extracted"
        extract_tar(tar, ex_dir)
    samples = organize_10x(ex_dir)
    print("10x sample groups found:", len(samples))
    if not samples:
        raise SystemExit("no complete 10x matrix/barcodes/features groups")
    entries = [load_sample(s) for s in samples]
    print("loaded samples:", len(entries))

    X, genes, barcodes = concat_samples(entries)
    n0 = X.shape[0]
    print(f"raw cells: {n0}, genes: {len(genes)}")

    # QC
    ngenes = (X > 0).sum(axis=1).A1
    mito = [i for i, g in enumerate(genes) if g.startswith("MT-")]
    if mito:
        mito_frac = (X[:, mito].sum(axis=1).A1 / X.sum(axis=1).A1)
    else:
        mito_frac = np.zeros(n0)
    keep = (ngenes >= 200) & (mito_frac < 0.20)
    X = X[keep]
    genes_keep = np.asarray((X > 0).sum(axis=0)).ravel() >= 3
    X = X[:, genes_keep]
    genes = genes[genes_keep]
    n1 = X.shape[0]
    print(f"after QC: cells {n0} -> {n1}, genes {len(genes)}")
    if args.max_cells and n1 > args.max_cells:
        rng = np.random.default_rng(20260809)
        sub = rng.choice(n1, args.max_cells, replace=False)
        X = X[sub]

    # normalize CPM + log1p (sparse)
    cs = X.sum(axis=1).A1
    cs[cs == 0] = 1
    Xn = X.multiply(1.0 / cs[:, None]).multiply(1e4).tocsr()
    Xn.data = np.log1p(Xn.data)
    Xn.eliminate_zeros()

    # HVG by variance
    v = np.asarray(Xn.power(2).mean(axis=0)).ravel() - np.asarray(Xn.mean(axis=0)).ravel() ** 2
    hv_idx = np.argsort(v)[::-1][:3000]
    Xh = Xn[:, hv_idx]
    hvg = genes[hv_idx]

    # PCA via TruncatedSVD on sparse data
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import StandardScaler
    svd = TruncatedSVD(n_components=30, random_state=20260809)
    P = svd.fit_transform(Xh)
    P = StandardScaler().fit_transform(P)
    print("PCA done, explained variance:", round(float(svd.explained_variance_ratio_.sum()), 3))

    # marker scores (mean log1p expression of markers per cell, sparse-safe)
    gmap = {g: i for i, g in enumerate(genes)}
    scores = {}
    used = {}
    for ct, ms in CELL_MARKERS.items():
        idx = [gmap[m] for m in ms if m in gmap]
        if len(idx) >= 3:
            scores[ct] = np.asarray(Xn[:, idx].mean(axis=1)).ravel()
            used[ct] = ms
    S = pd.DataFrame(scores)
    celltype = S.idxmax(axis=1).to_numpy()
    maxscore = S.max(axis=1).to_numpy()
    celltype[maxscore < 0.05] = "Other"
    print("cell types:", pd.Series(celltype).value_counts().to_dict())

    # gene-by-celltype enrichment for FDPS + hubs
    rows = []
    gidx = {g: i for i, g in enumerate(genes)}
    for g in GENES:
        if g not in gidx:
            continue
        vals = np.asarray(Xn[:, gidx[g]].todense()).ravel()
        for ct in sorted(set(celltype)):
            m = celltype == ct
            if m.sum() < 3 or (~m).sum() < 3:
                continue
            p = stats.mannwhitneyu(vals[m], vals[~m]).pvalue
            rows.append({"gene": g, "celltype": ct, "n": int(m.sum()),
                         "mean_expr": float(vals[m].mean()),
                         "pct_expr": float((vals[m] > 0).mean() * 100),
                         "p": float(p)})
    res = pd.DataFrame(rows)
    pv = res["p"].dropna().to_numpy()
    order = np.argsort(pv)
    adj = np.empty(len(pv))
    adj[order] = np.minimum.accumulate(
        (pv[order] * len(pv) / np.arange(len(pv), 0, -1))[::-1])[::-1]
    res["padj"] = np.nan
    res.loc[res["p"].notna(), "padj"] = adj
    res = res.sort_values("padj")
    out_dir = OUT / "scRNA_GSE161529"
    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "gene_celltype_enrichment.csv", index=False, encoding="utf-8-sig")
    mean_by_ct = pd.DataFrame(
        {g: np.asarray(Xn[:, gidx[g]].mean(axis=1)).ravel() for g in GENES if g in gidx}
    ).groupby(celltype).mean()
    mean_by_ct.to_csv(out_dir / "gene_mean_by_celltype.csv", encoding="utf-8-sig")

    # T cell - macrophage ligand-receptor co-expression
    tcell = (celltype == "T_cell") | (celltype == "CD8_T") | (celltype == "CD4_T") | (celltype == "Treg")
    mac = (celltype == "Macrophage") | (celltype == "Monocyte")
    lr_rows = []
    for lig, rec in LIGAND_RECEPTOR:
        if lig not in gidx or rec not in gidx:
            continue
        lv = np.asarray(Xn[:, gidx[lig]].todense()).ravel()
        rv = np.asarray(Xn[:, gidx[rec]].todense()).ravel()
        lr_rows.append({
            "ligand": lig, "receptor": rec,
            "ligand_in_T_pct": float((lv[tcell] > 0).mean() * 100),
            "receptor_in_Mac_pct": float((rv[mac] > 0).mean() * 100),
            "coexpressed_cell_pairs_frac": float(((lv > 0) & (rv > 0))[tcell | mac].mean() * 100),
        })
    lr = pd.DataFrame(lr_rows)
    lr.to_csv(out_dir / "ligand_receptor_t_mac.csv", index=False, encoding="utf-8-sig")
    print(lr.to_string(index=False))

    # t-SNE on subsample for figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.manifold import TSNE
    rng = np.random.default_rng(1)
    n_plot = min(30000, len(celltype))
    ii = rng.choice(len(celltype), n_plot, replace=False)
    emb = TSNE(n_components=2, perplexity=40, random_state=20260809,
               n_jobs=-1, init="pca").fit_transform(P[ii])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    cmap = plt.get_cmap("tab20")
    cts = sorted(set(celltype[ii]))
    for k, ct in enumerate(cts):
        m = celltype[ii] == ct
        axes[0].scatter(emb[m, 0], emb[m, 1], s=1.5, alpha=0.5,
                        color=cmap(k / max(len(cts), 1)), label=ct)
    axes[0].legend(frameon=False, fontsize=6.5, markerscale=4,
                   loc="center left", bbox_to_anchor=(1, 0.5))
    axes[0].set_title("Cell types (marker-based)")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    if "FDPS" in gidx:
        fv = np.asarray(Xn[ii, gidx["FDPS"]].todense()).ravel()
        sc = axes[1].scatter(emb[:, 0], emb[:, 1], s=1.5, c=fv, cmap="magma",
                             vmin=0, vmax=np.quantile(fv, 0.98))
        axes[1].set_title("FDPS expression (log1p CPM)")
        fig.colorbar(sc, ax=axes[1], shrink=0.8)
    axes[1].set_xticks([]); axes[1].set_yticks([])
    fig.suptitle(f"GSE161529 second single-cell atlas (n={n1:,} cells after QC)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig11_scrna_gse161529.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    top = res[res["gene"] == "FDPS"].head(5)
    mac_enr = res[(res["celltype"].isin(["Macrophage", "Monocyte"])) &
                  (res["gene"] == "FDPS")]
    t_enr = res[(res["celltype"].isin(["T_cell", "CD8_T", "CD4_T"])) &
                (res["gene"] == "FDPS")]
    mac_means = dict(zip(mac_enr["celltype"], mac_enr["mean_expr"]))
    patch = {
        "n_samples": len(entries),
        "n_cells_raw": int(n0),
        "n_cells_qc": int(n1),
        "n_genes": int(len(genes)),
        "fdps_macrophage_mean": float(mac_means.get("Macrophage", np.nan)),
        "fdps_monocyte_mean": float(mac_means.get("Monocyte", np.nan)),
        "fdps_tcell_mean": float(t_enr["mean_expr"].mean()) if len(t_enr) else np.nan,
        "top_fdps_celltypes": top[["celltype", "mean_expr", "padj"]].to_dict("records"),
        "ligand_receptor": lr.to_dict("records"),
        "table_dir": str(out_dir),
        "figure_path": str(FIG / "fig11_scrna_gse161529.png"),
    }
    save_patch(patch, "scrna2")
    print(res.head(25).to_string(index=False))
    print("DONE ->", out_dir, "|", FIG / "fig11_scrna_gse161529.png")


if __name__ == "__main__":
    main()
