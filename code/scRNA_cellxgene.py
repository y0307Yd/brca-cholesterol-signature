# -*- coding: utf-8 -*-
"""Hub-gene cell-type localization using CELLxGENE GSE176078 h5ad (author annotations)."""
import matplotlib; matplotlib.use('Agg')
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")

HUB = ["ABCG1","DHCR24","DHCR7","FDXR","G6PD","HMGCS2","HSD17B7","LIMA1","NSDHL","PRKAA1","VLDLR"]
OUT = Path(r"outputs\chol_metab_signature\scRNA"); OUT.mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(r"data\scRNA\GSE176078\GSE176078_cellxgene.h5ad")

# map gene names
var_map = dict(zip(adata.var_names, adata.var["feature_name"].astype(str)))
def find_gene(sym):
    hits=[g for g,v in var_map.items() if v==sym]
    return hits[0] if hits else None
hub_idx = {sym: find_gene(sym) for sym in HUB}
hub_idx = {k:v for k,v in hub_idx.items() if v is not None}
print("matched hub genes:", len(hub_idx), list(hub_idx.keys()))
if len(hub_idx) < 8: raise SystemExit("too few hub genes matched")

X = adata[:, list(hub_idx.values())].X
if hasattr(X, "toarray"): X = X.toarray()
expr = pd.DataFrame(X, index=adata.obs_names, columns=list(hub_idx.keys()))

for ctcol in ["celltype_major","celltype_minor"]:
    ct = adata.obs[ctcol].astype(str)
    rows=[]
    for g in hub_idx:
        v = expr[g].values
        for c in sorted(ct.unique()):
            m = (ct==c).values
            if m.sum()<3 or (~m).sum()<3:
                rows.append({"gene":g,"celltype":c,"mean":np.nan,"pct":np.nan,"p":np.nan}); continue
            p = stats.mannwhitneyu(v[m], v[~m]).pvalue
            rows.append({"gene":g,"celltype":c,"mean":v[m].mean(),"pct":(v[m]>0).mean()*100,"p":p})
    res = pd.DataFrame(rows)
    pv = res.p.dropna().values; order = np.argsort(pv)
    adj = np.empty(len(pv)); adj[order] = np.minimum.accumulate((pv[order]*len(pv)/np.arange(len(pv),0,-1))[::-1])[::-1]
    res["padj"] = np.nan; res.loc[res.p.notna(),"padj"] = adj
    res.sort_values("padj").to_csv(OUT/f"hub_gene_{ctcol}_enrichment.csv", index=False, encoding="utf-8-sig")
    expr.assign(celltype=ct.values).groupby("celltype")[list(hub_idx.keys())].mean().to_csv(OUT/f"hub_gene_{ctcol}_by_celltype.csv", encoding="utf-8-sig")
    print(f"== {ctcol}: top 12 ==")
    print(res.sort_values("padj").head(12).to_string(index=False))

# figures: dotplot + umap on major celltypes
adata.obs["celltype_major"] = adata.obs["celltype_major"].astype(str)
sc.settings.figdir = str(OUT); sc.settings.dpi = 300
sc.pl.dotplot(adata, var_names=list(hub_idx.values()), groupby="celltype_major", save="_hub_dotplot.png", show=False, standard_scale="var")
sc.pl.umap(adata, color=list(hub_idx.values())[:6], save="_hub_umap.png", show=False, ncols=3)
# per-celltype UMAP with hub signature (mean of matched hub genes)
sig = expr[list(hub_idx.keys())].mean(axis=1)
adata.obs["hub_signature"] = sig.values
sc.pl.umap(adata, color=["hub_signature","celltype_major"], save="_hub_sig.png", show=False, ncols=2)
print("DONE ->", OUT)