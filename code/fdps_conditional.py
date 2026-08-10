#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Conditional (COJO-style) analysis at the FDPS locus.

Conditional z-scores for a SNP given a conditioning SNP:
    z_cond = (z - r * z_cond_snp) / sqrt(1 - r^2)
with r the signed genotype correlation in 1000 Genomes EUR.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"work")
from coloc_locus_plot import read_fam, read_bed_snps, read_bim_region

BFILE = r"data\smr\g1000_eur\g1000_eur"
REGION = r"outputs\chol_metab_signature\coloc_fdps_region_ld.csv"
OUT_CSV = r"outputs\chol_metab_signature\fdps_conditional_analysis.csv"

df = pd.read_csv(REGION)
print("region SNPs:", len(df))
df["z_GWAS"] = df["b_GWAS"] / df["se_GWAS"]
df["z_eQTL"] = df["b_eQTL"] / df["se_eQTL"]

# load genotypes for the region SNPs from 1000G EUR
fam = read_fam(BFILE + ".fam")
n = len(fam)
hits = read_bim_region(BFILE + ".bim",
                       start=int(df["BP"].min()) - 50000,
                       end=int(df["BP"].max()) + 50000)
bim = pd.DataFrame(hits, columns=["row", "SNP", "BP", "A1", "A2"])
bim = bim.drop_duplicates(subset="SNP", keep="first")
bim_sel = bim[bim["SNP"].isin(set(df["SNP"]))]
dos = read_bed_snps(BFILE + ".bed", bim_sel["row"].tolist(), n)
rsids = bim_sel["SNP"].tolist()
print("genotypes loaded for", len(rsids), "SNPs")

# signed pairwise correlation with conditioning SNPs (complete cases)
def signed_r_with(dos, rsids, target):
    i = rsids.index(target)
    x = dos[i]
    out = np.full(dos.shape[0], np.nan)
    for j in range(dos.shape[0]):
        ok = ~np.isnan(x) & ~np.isnan(dos[j])
        if ok.sum() < 50:
            continue
        out[j] = np.corrcoef(x[ok], dos[j, ok])[0, 1]
    return out


def cond_z(z, z_cond_snp, r):
    zc = (z - r * z_cond_snp) / np.sqrt(np.maximum(1 - r * r, 1e-6))
    return zc


def analyze(df, zcol, label, cond_snps):
    z = df[zcol].to_numpy()
    rows = []
    for target in cond_snps:
        if target not in rsids:
            print("skip", target, "not in genotype data")
            continue
        r = signed_r_with(dos, rsids, target)
        rmap = dict(zip(rsids, r))
        rv = df["SNP"].map(rmap).to_numpy()
        zt = z[df["SNP"] == target][0]
        zc = cond_z(z, zt, rv)
        pc = 2 * stats.norm.sf(np.abs(zc))
        tmp = df[["SNP", "BP"]].copy()
        tmp["z_uncond"] = z
        tmp["z_cond_on_" + target] = zc
        tmp["p_cond_on_" + target] = pc
        tmp["r_with_" + target] = rv
        rows.append(tmp)
        keep = tmp[pc < 5e-8]
        print(f"[{label}] conditioned on {target}: "
              f"{len(keep)} SNPs remain P<5e-8; "
              f"top: {keep.sort_values('p_cond_on_'+target)['SNP'].iloc[0] if len(keep) else '-'} "
              f"(P={keep['p_cond_on_'+target].min():.2g} if any)")
    merged = rows[0]
    for t in rows[1:]:
        merged = merged.merge(t.drop(columns=["BP", "z_uncond"]), on="SNP",
                              how="outer")
    return merged


gwas = analyze(df, "z_GWAS", "GWAS",
               ["rs4971059", "rs12091730"])
eqtl = analyze(df, "z_eQTL", "eQTL",
               ["rs12091730", "rs6677385"])

out = gwas.merge(eqtl.drop(columns=["BP"]), on="SNP", how="outer",
                 suffixes=("_GWAS", "_eQTL"))
out.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

# key summaries
def summ(df, zcol, cond, label):
    keep = df[df["p_cond_on_" + cond] < 5e-8]
    print(f"\n{label}: {len(keep)} genome-wide SNPs after conditioning on {cond}")
    print(keep.sort_values("p_cond_on_" + cond)[
        ["SNP", "BP", zcol, "r_with_" + cond, "p_cond_on_" + cond]
    ].head(8).to_string(index=False))

summ(gwas, "z_uncond", "rs4971059", "GWAS")
summ(gwas, "z_uncond", "rs12091730", "GWAS")
summ(eqtl, "z_uncond", "rs12091730", "eQTL")
summ(eqtl, "z_uncond", "rs6677385", "eQTL")
print("saved", OUT_CSV)
