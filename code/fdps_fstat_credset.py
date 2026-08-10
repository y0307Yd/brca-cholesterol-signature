#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""FDPS instrument strength (F) and 99% credible set from coloc SNP-PP.H4."""
import json
import numpy as np
import pandas as pd

OUT = r"outputs\chol_metab_signature"
df = pd.read_csv(OUT + r"\coloc_fdps_region_ld.csv")
df["z_eQTL"] = df["b_eQTL"] / df["se_eQTL"]

# strongest instrument SNP for FDPS expression
top = df.loc[df["z_eQTL"].abs().idxmax()]
F = top["z_eQTL"] ** 2
print("strongest eQTL instrument:", top["SNP"], "z =", round(top["z_eQTL"], 2),
      "F =", round(F, 1))

# credible set from SNP-level PP.H4 (eQTL x GWAS shared signal)
pp = df["snp_PP_H4"].to_numpy()
ppn = pp / pp.sum()
order = np.argsort(-ppn)
cs = []
acc = 0.0
for i in order:
    cs.append(df["SNP"].iloc[i])
    acc += ppn[i]
    if acc >= 0.99:
        break
cs_df = df.iloc[order][["SNP", "BP", "b_eQTL", "se_eQTL", "z_eQTL",
                        "b_GWAS", "se_GWAS", "snp_PP_H4"]].head(len(cs)).copy()
cs_df["snp_PP_H4_norm"] = ppn[order[:len(cs)]]
cs_df["cum_PP_H4"] = np.cumsum(ppn[order[:len(cs)]])
print("99% credible set size:", len(cs),
      "| cumulative PP.H4:", round(float(cs_df["cum_PP_H4"].iloc[-1]), 4))
print(cs_df.to_string(index=False))

cs_df.to_csv(OUT + r"\fdps_credible_set.csv", index=False,
             encoding="utf-8-sig")
json.dump({
    "top_eqtl_snp": str(top["SNP"]),
    "top_eqtl_z": float(top["z_eQTL"]),
    "instrument_F": float(F),
    "credible_set_99_size": len(cs),
    "credible_set_snps": cs,
    "cumulative_pp": float(cs_df["cum_PP_H4"].iloc[-1]),
}, open(OUT + r"\fdps_fstat_credset.json", "w", encoding="utf-8"), indent=2)
print("saved")
