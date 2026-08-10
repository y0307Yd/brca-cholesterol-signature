# -*- coding: utf-8 -*-
"""Step 1: WGCNA-style module detection + cholesterol-gene-set intersection."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats, cluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster

OUT = Path(r".\outputs\chol_metab_signature")
RNG = 20260804
np.random.seed(RNG)

X = np.load(OUT / "tcga_Xlog.npy")          # genes x samples
genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
tr = pd.read_csv(OUT / "tcga_traits.csv")

# ---- filter & variable-gene selection ----
present = (X >= 1.0).mean(axis=1) >= 0.10
keep = np.where(present)[0]
mad = np.median(np.abs(X[keep] - np.median(X[keep], axis=1, keepdims=True)), axis=1)
top = np.argsort(mad)[::-1][:8000]
idx = keep[top]
Xg = X[idx].astype(np.float64)
gsel = [genes[i] for i in idx]
print("network genes:", len(gsel))

# ---- soft threshold ----
def sf_fit(beta):
    A = np.abs(np.corrcoef(Xg)) ** beta
    k = A.sum(axis=1) - 1.0
    k = k[k > 0]
    counts, edges = np.histogram(k, bins=100)
    centers = (edges[:-1] + edges[1:]) / 2
    m = counts > 0
    if m.sum() < 5:
        return 0.0, np.nan
    lk = np.log10(centers[m]); lp = np.log10(counts[m] / counts.sum())
    slope, intercept, r, p, se = stats.linregress(lk, lp)
    return r ** 2, slope

sf_rows = []
for beta in range(1, 21):
    r2, sl = sf_fit(beta)
    sf_rows.append({"beta": beta, "scale_free_R2": round(r2, 3), "slope": round(sl, 3) if sl == sl else np.nan})
sf = pd.DataFrame(sf_rows)
sf.to_csv(OUT / "soft_threshold.csv", index=False, encoding="utf-8-sig")
ok = sf[(sf.scale_free_R2 >= 0.80) & (sf.slope <= -0.5)]
BETA = int(ok.beta.min()) if len(ok) else 6
if len(ok) == 0:
    BETA = int(sf.loc[sf.scale_free_R2.idxmax(), "beta"]) if sf.scale_free_R2.max() >= 0.5 else 6
print(sf.to_string(index=False))
print("chosen beta:", BETA)

# ---- adjacency & TOM ----
tom_path = OUT / "dissTOM.npy"
if tom_path.exists():
    diss = np.load(tom_path)
    print("loaded saved dissTOM")
else:
    corr = np.corrcoef(Xg)
    A = np.abs(corr) ** BETA
    np.fill_diagonal(A, 0)
    k = A.sum(axis=1)
    denom = np.minimum.outer(k, k) + 1.0 - A
    num = A + A @ A
    np.fill_diagonal(num, 1)
    TOM = num / denom
    np.fill_diagonal(TOM, 1)
    diss = (1.0 - TOM).astype(np.float32)
    np.save(tom_path, diss)
print("TOM ready")

# ---- hierarchical clustering & module cutting ----
sq = squareform(diss, checks=False)
link = cluster.hierarchy.average(sq)
# height scan: prefer a partition with 10-40 modules of size >= 20
heights = np.linspace(0.3, 1.5, 121)
cand = []
for h in heights:
    lab = fcluster(link, t=float(h), criterion="distance")
    sizes = np.bincount(lab)
    nmod = int((sizes >= 20).sum())
    if 8 <= nmod <= 40:
        cand.append((float(h), nmod, lab))
if cand:
    cand.sort(key=lambda x: (abs(x[1] - 20), x[0]))
    CUT, nmod0, lab0 = cand[0]
else:
    CUT, nmod0, lab0 = 0.9, 0, fcluster(link, t=0.9, criterion="distance")
print("cut height:", CUT, "initial modules (size>=20):", nmod0, "total clusters:", len(set(lab0)))

# ---- merge modules by eigengene correlation >= 0.75 ----
def eigengene(Xm):
    Xc = Xm - Xm.mean(axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(Xc, full_matrices=False)
    return vt[0] * s[0]

mods = {m: np.where(lab0 == m)[0] for m in np.unique(lab0)}
mods = {m: v for m, v in mods.items() if len(v) >= 30}
changed = True
while changed:
    changed = False
    keys = sorted(mods)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            if a not in mods or b not in mods: continue
            ea = eigengene(Xg[mods[a]])
            eb = eigengene(Xg[mods[b]])
            r = np.corrcoef(ea, eb)[0, 1]
            if r >= 0.80:
                mods[a] = np.concatenate([mods[a], mods[b]])
                del mods[b]; changed = True; break
        if changed: break

modules = {}
for m, ix in mods.items():
    modules["M%d" % m] = ix
print("final modules:", len(modules), "sizes:", {m: len(v) for m, v in modules.items()})

# ---- module eigengenes & trait association ----
traits = {
    "event": tr["event"].values.astype(float),
    "time": tr["time"].values,
    "ER": tr["ER"].values,
    "stage": tr["stage"].values,
}
me_rows, mtraits = [], {}
for m, ix in modules.items():
    e = eigengene(Xg[ix])
    mtraits[m] = e
    kme = np.corrcoef(Xg[ix], e)[:-1, -1]
    me_rows.append({"module": m, "n_genes": len(ix),
                    "kME_mean": float(np.mean(kme)), "kME_min": float(np.min(kme))})
ME = pd.DataFrame(me_rows)
ME.to_csv(OUT / "modules_summary.csv", index=False, encoding="utf-8-sig")
np.save(OUT / "module_eigengenes.npy",
        np.vstack([mtraits[m] for m in sorted(modules)]))

mod_trait = []
for m in sorted(modules):
    e = mtraits[m]
    for tn, tv in traits.items():
        mask = ~np.isnan(tv)
        if mask.sum() < 50: continue
        if tn == "time" or tn == "stage":
            r, p = stats.spearmanr(e[mask], tv[mask])
        else:
            r, p = stats.pearsonr(e[mask], tv[mask])
        mod_trait.append({"module": m, "trait": tn, "n": int(mask.sum()),
                          "cor": round(float(r), 3), "p": float(p)})
mt = pd.DataFrame(mod_trait)
# BH within trait
for tn in mt.trait.unique():
    ix = mt.trait == tn
    pv = mt.loc[ix, "p"].values
    order = np.argsort(pv); ranked = pv[order]
    adj = np.minimum.accumulate(ranked[::-1] * len(ranked) / np.arange(len(ranked), 0, -1))[::-1]
    mt.loc[ix, "padj"] = np.clip(adj[np.argsort(order)], 0, 1)
mt.to_csv(OUT / "module_trait.csv", index=False, encoding="utf-8-sig")
print(mt.sort_values("p").head(20).to_string(index=False))

# ---- gene-module map ----
map_rows = []
for m, ix in modules.items():
    for i in ix:
        map_rows.append((gsel[i], m))
gene_mod = pd.DataFrame(map_rows, columns=["gene", "module"])
gene_mod.to_csv(OUT / "gene_module_map.csv", index=False, encoding="utf-8-sig")
np.save(OUT / "linkage.npy", link)

# ---- pick disease module (event) and ER module (fallback phenotype) ----
ev = mt[mt.trait == "event"].sort_values("p")
er = mt[mt.trait == "ER"].sort_values("p")
ev_mod = ev.iloc[0]["module"] if len(ev) else None
er_mod = er.iloc[0]["module"] if len(er) else None

chol = pd.read_csv(OUT / "cholesterol_genes.csv")["symbol"].tolist()
chols = set(chol)
mod_genes = {m: [gsel[i] for i in ix] for m, ix in modules.items()}
cand_ev = sorted(set(mod_genes[ev_mod]) & chols) if ev_mod else []
cand_er = sorted(set(mod_genes[er_mod]) & chols) if er_mod else []
use_mod = ev_mod if (ev_mod and cand_ev) else er_mod
cand = cand_ev if cand_ev else cand_er
print("cholesterol genes in network:", len(chols & set(gsel)))
print("event module:", ev_mod, "candidates:", cand_ev)
print("ER module:", er_mod, "candidates:", cand_er)
print("used module:", use_mod, "| candidates:", cand)

pd.DataFrame({"module": ev_mod, "gene": mod_genes[ev_mod]}).to_csv(
    OUT / "event_module_genes.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"module": er_mod, "gene": mod_genes[er_mod]}).to_csv(
    OUT / "er_module_genes.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"gene": cand}).to_csv(OUT / "candidate_genes.csv", index=False, encoding="utf-8-sig")
(OUT / "disease_module_genes.csv").unlink(missing_ok=True)
json.dump({"beta": BETA, "cut": CUT, "network_genes": len(gsel), "n_modules": len(modules),
           "event_module": ev_mod,
           "event_module_padj": float(ev.iloc[0]["padj"]) if len(ev) else None,
           "event_module_genes": len(mod_genes[ev_mod]) if ev_mod else 0,
           "er_module": er_mod,
           "er_module_padj": float(er.iloc[0]["padj"]) if len(er) else None,
           "candidates_event_module": cand_ev, "candidates_er_module": cand_er,
           "used_module": use_mod, "candidates": cand},
          open(OUT / "wgcna_summary.json", "w", encoding="utf-8"), indent=1)
print("saved wgcna outputs")