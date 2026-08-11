# -*- coding: utf-8 -*-
"""Round-3 bonus: pathway-activity panel across cholesterol subtypes (v22+).

Rank-based single-sample pathway activity (ssGSEA-style, alpha = 0.25) for
cholesterol biosynthesis/metabolism, fatty-acid oxidation, glycolysis,
oxidative phosphorylation, cell cycle, estrogen response and immune
cytokine programs, compared across the four hub-gene subtypes.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)

# curated gene sets (all resolvable from the TCGA gene list)
FAO = ["CPT1A", "CPT1B", "CPT1C", "CPT2", "SLC25A20", "ACADM", "ACADL",
       "ACADS", "ACADVL", "HADHA", "HADHB", "ACOX1", "ACOX2", "ACOX3",
       "EHHADH", "DECR1", "ECH1", "ECHS1", "HADH", "ACADSB", "ACAD11",
       "ACAD10", "ACAD9", "ACAD8", "ACAA1", "ACAA2"]
GLY = ["HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA",
       "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2", "PGAM1",
       "PGAM2", "ENO1", "ENO2", "ENO3", "PKLR", "PKM", "LDHA", "LDHB",
       "PDK1", "PDK2", "PDK3", "PDK4", "SLC2A1", "SLC2A3"]
ER_RESP = ["ESR1", "ESR2", "PGR", "TFF1", "TFF3", "GREB1", "CA12",
           "CCND1", "MYC", "STC2", "AGR2", "GATA3", "FOXA1", "XBP1",
           "MUC1", "AREG", "ERBB2", "NRIP1", "SCUBE2"]
CYTOKINE = ["CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL2", "CCL3", "CCL4",
            "CCL5", "CCL19", "CCL20", "CX3CL1", "IL1A", "IL1B", "IL6",
            "IL8", "IL10", "IL12A", "IL12B", "IL15", "IL18", "TNF",
            "IFNG", "CSF1", "CSF2", "CSF3", "CXCL8"]
OXPHOS = ["NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6",
          "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12",
          "NDUFA13", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5",
          "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11",
          "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6",
          "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3", "SDHA",
          "SDHB", "SDHC", "SDHD", "UQCRC1", "UQCRC2", "UQCRB", "UQCRQ",
          "UQCR10", "UQCR11", "UQCRFS1", "COX4I1", "COX5A", "COX5B",
          "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX7C",
          "COX8A", "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D",
          "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PB",
          "ATP5PD", "ATP5PO", "ATP5PF"]
CELLCYCLE = ["CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CCNA1", "CCNA2",
             "CCNB1", "CCNB2", "CCNC", "CCND1", "CCND2", "CCND3",
             "CCNE1", "CCNE2", "CCNF", "MCM2", "MCM3", "MCM4", "MCM5",
             "MCM6", "MCM7", "MCM10", "CDC20", "CDC25A", "CDC25B",
             "CDC25C", "CDC6", "CDT1", "BUB1", "BUB1B", "BUB3", "AURKA",
             "AURKB", "CHEK1", "CHEK2", "E2F1", "E2F2", "RB1", "PLK1",
             "PLK4", "CKS1B", "CKS2"]


def build_sets(genes):
    gene_set = {
        "Cholesterol biosynthesis (KEGG)": None,   # filled from csv below
        "Cholesterol metabolism (GO)": None,
        "Fatty acid oxidation": FAO,
        "Glycolysis": GLY,
        "Estrogen response": ER_RESP,
        "Immune cytokines": CYTOKINE,
        "Oxidative phosphorylation": OXPHOS,
        "Cell cycle": CELLCYCLE,
    }
    chol = pd.read_csv(OUT / "cholesterol_genes.csv")
    long_assoc = pd.read_csv(OUT / "cholesterol_genes_long.csv")
    kegg = set()
    go_metab = set()
    for _, r in long_assoc.iterrows():
        src = str(r.get("source", "")).strip('"')
        sym = str(r.get("SYMBOL", "")).strip('"')
        if "KEGG" in src:
            kegg.add(sym)
        if "0008203" in src:
            go_metab.add(sym)
    gene_set["Cholesterol biosynthesis (KEGG)"] = sorted(kegg)
    gene_set["Cholesterol metabolism (GO)"] = sorted(go_metab)
    idx = {g: i for i, g in enumerate(genes)}
    out = {}
    for name, gs in gene_set.items():
        gi = [idx[g] for g in gs if g in idx]
        if len(gi) < 3:
            print(f"SKIP {name}: only {len(gi)} genes")
            continue
        out[name] = gi
        print(f"{name}: {len(gi)}/{len(gs)} genes mapped")
    return out


def ssgsea_scores(X, gene_sets_idx, alpha=0.25):
    n_genes, n_samples = X.shape
    order = np.argsort(-X, axis=0, kind="stable")
    scores = {}
    for name, gi in gene_sets_idx.items():
        gi = np.asarray(gi)
        weights = np.abs(X) ** alpha
        denom = weights[gi].sum(axis=0) + 1e-12
        es = np.zeros(n_samples)
        for s in range(n_samples):
            ord_s = order[:, s]
            w_sorted = weights[ord_s, s]
            in_set = np.zeros(n_genes, dtype=bool)
            in_set[gi] = True
            in_sorted = in_set[ord_s]
            hit = np.cumsum(np.where(in_sorted, w_sorted, 0.0)) / denom[s]
            miss = np.cumsum(np.where(in_sorted, 0.0, 1.0)) / (n_genes - len(gi))
            es[s] = float(np.max(hit - miss))
        scores[name] = es
    return scores


def main():
    genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8")
             if l.strip()]
    X = np.load(OUT / "tcga_Xlog.npy")
    sub = pd.read_csv(OUT / "tcga_subtypes.csv")
    sets = build_sets(genes)
    scores = ssgsea_scores(X, sets)
    df = pd.DataFrame(scores)
    df["subtype"] = sub["subtype"].to_numpy()

    ks = [1, 2, 3, 4]
    rows = []
    for name in scores:
        groups = [df.loc[df["subtype"] == k, name].to_numpy() for k in ks]
        kw = stats.kruskal(*groups)
        rows.append({
            "pathway": name,
            "kruskal_P": kw.pvalue,
            **{f"mean_C{k}": float(df.loc[df["subtype"] == k, name].mean())
               for k in ks},
        })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "pathway_activity_by_subtype.csv", index=False,
               encoding="utf-8-sig")
    print(res.round(4).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mat = res[[f"mean_C{k}" for k in ks]].to_numpy()
    names = res["pathway"].tolist()
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto",
                   vmin=-np.abs(mat).max(), vmax=np.abs(mat).max())
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"C{k}" for k in ks], fontsize=10)
    ax.set_title("Pathway-activity ssGSEA scores by cholesterol subtype")
    fig.colorbar(im, ax=ax, label="mean enrichment score", shrink=0.8)
    fig.tight_layout()
    fig.savefig(FIG / "supp_fig7_pathway_panel.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print("saved supp_fig7_pathway_panel.png")


if __name__ == "__main__":
    main()
