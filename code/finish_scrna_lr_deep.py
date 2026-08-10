# -*- coding: utf-8 -*-
"""Deeper T-myeloid ligand-receptor analysis in GSE161529.

Reuses the extracted 10x matrices, QC/normalisation and marker-based
annotation from finish_scrna_gse161529.py but skips t-SNE and expands the
ligand-receptor list (20 pairs), adding FDPS co-expression context.

Usage:
  C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_scrna_lr_deep.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from finish_scrna_gse161529 import CELL_MARKERS, concat_samples, load_sample, organize_10x
from online_utils import OUT, save_patch

EX_DIR = Path(r"data\scRNA\GSE161529\extracted")
OUT_DIR = OUT / "scRNA_GSE161529"

LIGAND_RECEPTOR = [
    ("SPP1", "CD44"), ("CD74", "MIF"), ("LGALS9", "HAVCR2"),
    ("CCL18", "CCR8"), ("CXCL16", "CXCR6"),
    ("CCL5", "CCR5"), ("CCL4", "CCR5"), ("CCL3", "CCR5"),
    ("CXCL9", "CXCR3"), ("CXCL10", "CXCR3"), ("CXCL11", "CXCR3"),
    ("CCL2", "CCR2"), ("CCL13", "CCR2"), ("CCL7", "CCR2"),
    ("CX3CL1", "CX3CR1"), ("IL15", "IL2RB"), ("IL18", "IL18R1"),
    ("CSF1", "CSF1R"), ("IL1B", "IL1R2"), ("TNFSF13B", "TNFRSF17"),
]


def main():
    samples = organize_10x(EX_DIR)
    print("10x sample groups found:", len(samples))
    entries = [load_sample(s) for s in samples]
    X, genes, barcodes = concat_samples(entries)
    n0 = X.shape[0]
    ngenes = (X > 0).sum(axis=1).A1
    mito = [i for i, g in enumerate(genes) if g.startswith("MT-")]
    mito_frac = (X[:, mito].sum(axis=1).A1 / X.sum(axis=1).A1) if mito else np.zeros(n0)
    keep = (ngenes >= 200) & (mito_frac < 0.20)
    X = X[keep]
    genes_keep = np.asarray((X > 0).sum(axis=0)).ravel() >= 3
    X = X[:, genes_keep]
    genes = genes[genes_keep]
    n1 = X.shape[0]
    cs = X.sum(axis=1).A1
    cs[cs == 0] = 1
    Xn = X.multiply(1.0 / cs[:, None]).multiply(1e4).tocsr()
    Xn.data = np.log1p(Xn.data)
    Xn.eliminate_zeros()
    print(f"after QC: {n1} cells, {len(genes)} genes")

    gmap = {g: i for i, g in enumerate(genes)}
    scores = {}
    for ct, ms in CELL_MARKERS.items():
        idx = [gmap[m] for m in ms if m in gmap]
        if len(idx) >= 3:
            scores[ct] = np.asarray(Xn[:, idx].mean(axis=1)).ravel()
    S = pd.DataFrame(scores)
    celltype = S.idxmax(axis=1).to_numpy()
    maxscore = S.max(axis=1).to_numpy()
    celltype[maxscore < 0.05] = "Other"
    print("cell types:", pd.Series(celltype).value_counts().to_dict())

    tcell = np.isin(celltype, ["T_cell", "CD8_T", "CD4_T", "Treg"])
    mac = np.isin(celltype, ["Macrophage", "Monocyte"])
    rows = []
    for lig, rec in LIGAND_RECEPTOR:
        if lig not in gmap or rec not in gmap:
            continue
        lv = np.asarray(Xn[:, gmap[lig]].todense()).ravel()
        rv = np.asarray(Xn[:, gmap[rec]].todense()).ravel()
        rows.append({
            "ligand": lig, "receptor": rec,
            "ligand_in_T_pct": float((lv[tcell] > 0).mean() * 100),
            "receptor_in_Mac_pct": float((rv[mac] > 0).mean() * 100),
            "coexpressed_cell_pairs_frac": float(((lv > 0) & (rv > 0))[tcell | mac].mean() * 100),
            "ligand_in_Mac_pct": float((lv[mac] > 0).mean() * 100),
            "receptor_in_T_pct": float((rv[tcell] > 0).mean() * 100),
        })
    lr = pd.DataFrame(rows).sort_values("coexpressed_cell_pairs_frac", ascending=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lr.to_csv(OUT_DIR / "ligand_receptor_t_mac_expanded.csv", index=False,
              encoding="utf-8-sig")
    print(lr.to_string(index=False))

    # FDPS co-expression context: % FDPS+ cells within T and myeloid
    fdps = np.asarray(Xn[:, gmap["FDPS"]].todense()).ravel()
    fdps_high = fdps > np.quantile(fdps[fdps > 0], 0.75) if (fdps > 0).any() else fdps > 0
    summary = {
        "n_cells_qc": int(n1),
        "n_t_cells": int(tcell.sum()),
        "n_myeloid": int(mac.sum()),
        "fdps_pct_t": float((fdps[tcell] > 0).mean() * 100),
        "fdps_pct_myeloid": float((fdps[mac] > 0).mean() * 100),
        "fdps_mean_t": float(fdps[tcell].mean()),
        "fdps_mean_myeloid": float(fdps[mac].mean()),
        "top_pairs": lr.head(10).to_dict("records"),
    }
    save_patch(summary, "scrna_lr_deep")
    print("DONE ->", OUT_DIR / "ligand_receptor_t_mac_expanded.csv")


if __name__ == "__main__":
    main()
