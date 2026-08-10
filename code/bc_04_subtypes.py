# -*- coding: utf-8 -*-
"""Step 3: consensus clustering of hub genes; subtype prognosis/immune/pathway characterization."""
import os
os.environ["LOKY_MAX_CPU_COUNT"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"
import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")
OUT = Path(r".\outputs\chol_metab_signature")
RNG = 20260804
rng = np.random.default_rng(RNG)

X = np.load(OUT / "tcga_Xlog.npy")
genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
tr = pd.read_csv(OUT / "tcga_traits.csv")
hub = pd.read_csv(OUT / "hub_genes.csv")["gene"].tolist()
idx = {g: i for i, g in enumerate(genes)}
Z = np.vstack([(X[idx[g]] - X[idx[g]].mean()) / X[idx[g]].std() for g in hub]).T  # 952 x hub

# ---------------- consensus clustering ----------------
n = len(tr)
pac = {}
cons = {}
for k in range(2, 7):
    C = np.zeros((n, n)); W = np.zeros((n, n))
    for b in range(100):
        s = rng.choice(n, size=int(0.8 * n), replace=False)
        km = KMeans(n_clusters=k, n_init=10, random_state=RNG + b)
        lab = km.fit_predict(Z[s])
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                if lab[i] == lab[j]:
                    C[s[i], s[j]] += 1
                W[s[i], s[j]] += 1
    M = C / np.maximum(W, 1)
    M = (M + M.T) / 2
    np.fill_diagonal(M, 1)
    cons[k] = M
    off = M[np.triu_indices(n, 1)]
    pac[k] = float(((off > 0.1) & (off < 0.9)).mean())
k_best = min(pac, key=pac.get)
M = cons[k_best]
lab = fcluster(linkage(squareform(1 - M), method="average"), t=k_best, criterion="maxclust")
print("PAC:", {k: round(v, 4) for k, v in pac.items()}, "| chosen k:", k_best)
print("cluster sizes:", {c: int((lab == c).sum()) for c in np.unique(lab)})

tr = tr.copy()
tr["subtype"] = lab
tr.to_csv(OUT / "tcga_subtypes.csv", index=False, encoding="utf-8-sig")
json.dump({"pac": pac, "k_best": k_best},
          open(OUT / "consensus_summary.json", "w", encoding="utf-8"), indent=1)

# ---------------- survival by subtype ----------------
ks = np.unique(lab)
logp = None
if len(ks) == 2:
    lr = logrank_test(tr["time"][lab == ks[0]], tr["time"][lab == ks[1]],
                      tr["event"][lab == ks[0]], tr["event"][lab == ks[1]])
    logp = float(lr.p_value)
    print("log-rank p:", logp)

# ---------------- immune marker scores ----------------
markers = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG"],
    "CD4_T": ["CD4", "IL7R", "CCR7", "CD40LG", "ICOS"],
    "Treg": ["FOXP3", "CTLA4", "IL2RA", "IKZF2"],
    "NK": ["NKG7", "KLRD1", "KLRB1", "NCR1", "GNLY"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BLK"],
    "Macrophage": ["CD68", "CD163", "CSF1R", "LYZ", "ITGAX"],
    "DC": ["ITGAX", "CD1C", "CLEC9A", "BATF3", "FLT3"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "CEACAM8"],
    "Stroma": ["COL1A1", "COL3A1", "ACTA2", "PDGFRB", "FAP"],
}
score = {}
for cell, gs in markers.items():
    gs = [g for g in gs if g in idx]
    if len(gs) < 3:
        continue
    e = np.vstack([(X[idx[g]] - X[idx[g]].mean()) / X[idx[g]].std() for g in gs])
    score[cell] = e.mean(axis=0)
imm = pd.DataFrame(score, index=tr["sample_id"])
imm["subtype"] = lab
imm.to_csv(OUT / "immune_scores.csv", index=False, encoding="utf-8-sig")

imm_sum = imm.groupby("subtype")[list(score)].mean().T
krus = {}
for c in score:
    groups = [imm[c][lab == k] for k in ks]
    krus[c] = stats.kruskal(*groups).pvalue
imm_sum["kruskal_p"] = [krus[c] for c in imm_sum.index]
imm_sum.to_csv(OUT / "immune_by_subtype.csv", encoding="utf-8-sig")
print(imm_sum.round(3).to_string())

# cholesterol signature / pathway scores by subtype
sig = np.load(OUT / "signature_tcga.npy")
chol_genes = pd.read_csv(OUT / "cholesterol_genes.csv")["symbol"].tolist()
ci = [g for g in chol_genes if g in idx]
ce = np.vstack([(X[idx[g]] - X[idx[g]].mean()) / X[idx[g]].std() for g in ci])
path_score = ce.mean(axis=0)
sub_stats = pd.DataFrame({"subtype": lab, "signature": sig, "chol_pathway": path_score,
                          "ER": tr["ER"].values, "event": tr["event"].values, "time": tr["time"].values})
grp = sub_stats.groupby("subtype")
tbl = pd.DataFrame({
    "n": grp.size(),
    "ER_pos_frac": grp["ER"].mean(),
    "signature_mean": grp["signature"].mean(),
    "chol_pathway_mean": grp["chol_pathway"].mean(),
})
tbl.to_csv(OUT / "subtype_stats.csv", encoding="utf-8-sig")
print(tbl.round(3).to_string())

# ---------------- enrichment (cluster A vs B) ----------------
def bh(pv):
    pv = np.asarray(pv, float); order = np.argsort(pv); ranked = pv[order]
    adj = np.minimum.accumulate(ranked[::-1] * len(ranked) / np.arange(len(ranked), 0, -1))[::-1]
    return np.clip(adj[np.argsort(order)], 0, 1)

a = lab == ks[0]; b = lab == ks[1]
pv = np.full(len(genes), np.nan); fc = np.full(len(genes), np.nan)
for i in range(len(genes)):
    pv[i], _ = stats.mannwhitneyu(X[i, a], X[i, b], alternative="two-sided")
    fc[i] = X[i, a].mean() - X[i, b].mean()
padj = bh(pv[~np.isnan(pv)])
de = pd.DataFrame({"gene": genes, "fc_A_minus_B": fc, "p": pv, "padj": np.nan})
de.loc[~np.isnan(pv), "padj"] = padj
de["direction"] = np.where(de.fc_A_minus_B > 0, "up_in_A", "up_in_B")
de.to_csv(OUT / "subtype_de_genes.csv", index=False, encoding="utf-8-sig")
sig_genes = de[de.padj < 0.05].copy()
print("DE genes padj<0.05:", len(sig_genes))

def enrich(gene_sets, anno, bg, direction, label):
    rows = []
    for term, gs in gene_sets.items():
        gs = set(gs)
        ov = len(gs & set(sig_genes.gene))
        if ov == 0:
            continue
        pop = len(bg); drawn = len(sig_genes); hits_in_pop = len(gs & bg)
        p = stats.hypergeom.sf(ov - 1, pop, hits_in_pop, drawn)
        rows.append({"term": term, "label": label, "direction": direction, "overlap": ov,
                     "pathway_size": hits_in_pop, "p": float(p)})
    return rows

kg = pd.read_csv(OUT / "kegg_annotations.csv")
go = pd.read_csv(OUT / "go_annotations.csv")
bg = set(genes)
kegg_sets = {t: set(g.symbol) for t, g in kg.groupby("kegg_id") if len(g) >= 3}
go_sets = {t: set(g.symbol) for t, g in go[go.ontology == "BP"].groupby("go_id") if len(g) >= 5}
rows = []
for d in ["up_in_A", "up_in_B"]:
    sub = sig_genes[sig_genes.direction == d]
    if len(sub) == 0:
        continue
    rows += enrich(kegg_sets, kg, bg, d, "KEGG")
    rows += enrich(go_sets, go, bg, d, "GO_BP")
enr = pd.DataFrame(rows)
if len(enr):
    enr["padj"] = bh(enr.p.values)
    enr = enr.sort_values("padj")
    enr.to_csv(OUT / "subtype_enrichment.csv", index=False, encoding="utf-8-sig")
    print(enr.head(25).to_string(index=False))
json.dump({"k_best": k_best, "logrank_p": logp, "n_de": int(len(sig_genes)),
           "n_enriched_terms": int(len(enr)) if len(enr) else 0},
          open(OUT / "subtype_summary.json", "w", encoding="utf-8"), indent=1)
print("saved step3 outputs")