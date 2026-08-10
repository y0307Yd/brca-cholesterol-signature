# -*- coding: utf-8 -*-
"""Enrichment v2: DE genes with |mean diff| >= 1 (log2 scale) per subtype vs rest."""
import os
os.environ["LOKY_MAX_CPU_COUNT"] = "4"
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(r".\outputs\chol_metab_signature")
de = pd.read_csv(OUT / "subtype_de_genes.csv")
kg = pd.read_csv(OUT / "kegg_annotations.csv")
go = pd.read_csv(OUT / "go_annotations.csv")
genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
bg = set(genes)
kegg_sets = {t: set(g.symbol) for t, g in kg.groupby("kegg_id") if len(g) >= 3}
go_sets = {t: set(g.symbol) for t, g in go[go.ontology == "BP"].groupby("go_id") if 10 <= len(g) <= 500}

def bh(pv):
    pv = np.asarray(pv, float); order = np.argsort(pv); ranked = pv[order]
    adj = np.minimum.accumulate(ranked[::-1] * len(ranked) / np.arange(len(ranked), 0, -1))[::-1]
    return np.clip(adj[np.argsort(order)], 0, 1)

rows = []
for c in sorted(de.cluster.unique()):
    d = de[de.cluster == c]
    sig = d[(d.padj < 0.05) & (d.fc_cluster_minus_rest.abs() >= 1.0)]
    print("cluster", int(c), "DE with |fc|>=1:", len(sig))
    for direction, sign in [("up", 1), ("down", -1)]:
        sset = set(sig[sig.direction == direction].gene)
        if len(sset) < 5:
            continue
        for label, sets in [("KEGG_C%d" % int(c), kegg_sets), ("GOBP_C%d" % int(c), go_sets)]:
            for term, gs in sets.items():
                ov = len(gs & sset)
                if ov < 3:
                    continue
                p = stats.hypergeom.sf(ov - 1, len(bg), len(gs & bg), len(sset))
                rows.append({"label": label, "direction": direction, "term": term,
                             "overlap": ov, "pathway_size": len(gs & bg), "n_de": len(sset),
                             "p": float(p)})
enr = pd.DataFrame(rows)
enr["padj"] = bh(enr.p.values)
enr = enr.sort_values("padj")
enr.to_csv(OUT / "subtype_enrichment_fc1.csv", index=False, encoding="utf-8-sig")
print("total terms:", len(enr))
print(enr.head(40).to_string(index=False))