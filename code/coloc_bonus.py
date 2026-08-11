# -*- coding: utf-8 -*-
"""Colocalisation for LIPA / FAXDC2 / SREBF1 (eQTLGen x BCAC overall).

Same pipeline as coloc_fdps.py: parse SMR --plot output, merge eQTL and GWAS
summary statistics by rsID, then coloc.abf-style posterior probabilities.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from coloc_fdps import coloc_abf, load_smr_plot

ROOT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat")
OUT = ROOT / "outputs" / "chol_metab_signature"
PLOT = ROOT / "plot"


def logsumexp(x):
    x = np.asarray(x, dtype=float)
    m = np.max(x)
    if not np.isfinite(m):
        return m
    return m + np.log(np.sum(np.exp(x - m)))


def log_subtract(x, y):
    """log(exp(x) - exp(y)) with x > y; -inf if x <= y."""
    x = float(x)
    y = float(y)
    if x <= y:
        return -np.inf
    return x + np.log1p(-np.exp(y - x))


def coloc_abf_stable(rsid, beta1, se1, beta2, se2,
                     n1, n2, type1="quant", type2="cc",
                     s2=None, sdY1=1.0, corrections=True,
                     p1=1e-4, p2=1e-4, p12=1e-5, w=0.2):
    """Numerically stable coloc.abf (all operations in log space)."""
    from coloc_fdps import wakefield_abf
    rsid = np.asarray(rsid, dtype=str)
    beta1 = np.asarray(beta1, dtype=float)
    se1 = np.asarray(se1, dtype=float)
    beta2 = np.asarray(beta2, dtype=float)
    se2 = np.asarray(se2, dtype=float)
    ok = (~np.isnan(beta1)) & (~np.isnan(se1)) & (~np.isnan(beta2)) & (~np.isnan(se2))
    ok &= (se1 > 0) & (se2 > 0)
    rsid, beta1, se1, beta2, se2 = rsid[ok], beta1[ok], se1[ok], beta2[ok], se2[ok]
    n = len(rsid)
    if n < 2:
        raise ValueError("fewer than 2 usable SNPs")
    if corrections:
        if type1 == "quant" and n1 is not None:
            se1 = np.sqrt(se1**2 + sdY1**2 / (2.0 * n1))
        if type2 == "cc" and s2 is not None and n2 is not None:
            se2 = np.sqrt(se2**2 + 1.0 / (2.0 * n2 * s2 * (1.0 - s2)))
    l1 = wakefield_abf(beta1, se1, w=w)
    l2 = wakefield_abf(beta2, se2, w=w)
    lA = logsumexp(l1)
    lB = logsumexp(l2)
    lC = logsumexp(l1 + l2)
    lH0 = 0.0
    lH1 = lA + np.log(p1) - np.log(n)
    lH2 = lB + np.log(p2) - np.log(n)
    lH3 = log_subtract(lA + lB, lC) + np.log(p1) + np.log(p2) - 2.0 * np.log(n)
    lH4 = lC + np.log(p12) - np.log(n)
    lp0 = np.log(max(1.0 - p1 - p2 - p12, np.finfo(float).tiny))
    logs = np.array([lp0, lH1, lH2, lH3, lH4])
    denom = logsumexp(logs)
    pp = np.exp(logs - denom)
    snp_pp_h4 = np.exp(l1 + l2 - logsumexp(l1 + l2))
    lead = int(np.argmax(l1 + l2))
    return {
        "n_snps": n,
        "PP.H0.abf": float(pp[0]),
        "PP.H1.abf": float(pp[1]),
        "PP.H2.abf": float(pp[2]),
        "PP.H3.abf": float(pp[3]),
        "PP.H4.abf": float(pp[4]),
        "lead_snp": rsid[lead],
        "lead_snp_pp_h4": float(snp_pp_h4[lead]),
        "n1": n1, "n2": n2, "s2": s2,
        "p1": p1, "p2": p2, "p12": p12, "w": w,
    }

GENES = {
    "LIPA":   {"ensg": "ENSG00000107798", "bp": 91_073_000},
    "FAXDC2": {"ensg": "ENSG00000170271", "bp": 154_218_000},
    "SREBF1": {"ensg": "ENSG00000072310", "bp": 17_727_000},
}


def run_gene(name, ensg, probe_bp):
    plot = PLOT / f"{name.lower()}.{ensg}.txt"
    if not plot.exists():
        print(f"[{name}] plot file missing: {plot}")
        return None
    snp_info, gwas, eqtl_by_probe, meta = load_smr_plot(plot)
    if ensg not in eqtl_by_probe:
        print(f"[{name}] {ensg} not in eQTL blocks: {list(eqtl_by_probe)[:8]}")
        return None
    eqtl = eqtl_by_probe[ensg]
    m = eqtl.merge(gwas, on="SNP", how="inner")
    m = m.merge(snp_info[["SNP", "Chr", "BP", "A1", "A2"]], on="SNP", how="left")
    m["BP"] = pd.to_numeric(m["BP"], errors="coerce")
    m = m.dropna(subset=["b_eQTL", "se_eQTL", "b_GWAS", "se_GWAS", "BP"])
    m = m[(m["se_eQTL"] > 0) & (m["se_GWAS"] > 0)]
    m = m.sort_values("BP").reset_index(drop=True)
    print(f"[{name}] eQTL SNPs: {len(eqtl)}, GWAS: {len(gwas)}, "
          f"merged usable: {len(m)}")
    if len(m) < 2:
        print(f"[{name}] too few SNPs")
        return None
    m["z_eQTL"] = m["b_eQTL"] / m["se_eQTL"]
    m["z_GWAS"] = m["b_GWAS"] / m["se_GWAS"]

    def run_one(frac, label):
        sub = m if frac is None else m[m["BP"].between(probe_bp - frac / 2.0,
                                                       probe_bp + frac / 2.0)]
        if len(sub) < 2:
            return None
        res = coloc_abf_stable(
            rsid=sub["SNP"].values,
            beta1=sub["b_eQTL"].values,
            se1=sub["se_eQTL"].values,
            beta2=sub["b_GWAS"].values,
            se2=sub["se_GWAS"].values,
            n1=31684, n2=247173, type1="quant", type2="cc",
            s2=0.2929, corrections=True, p1=1e-4, p2=1e-4,
            p12=1e-5, w=0.2,
        )
        res["window"] = label
        res["gene"] = name
        res["probe"] = ensg
        res["region_chr"] = str(sub["Chr"].iloc[0])
        res["region_start"] = int(sub["BP"].min())
        res["region_end"] = int(sub["BP"].max())
        lead = sub[sub["SNP"] == res["lead_snp"]].iloc[0]
        res["lead_b_eQTL"] = float(lead["b_eQTL"])
        res["lead_se_eQTL"] = float(lead["se_eQTL"])
        res["lead_z_eQTL"] = float(lead["z_eQTL"])
        res["lead_b_GWAS"] = float(lead["b_GWAS"])
        res["lead_se_GWAS"] = float(lead["se_GWAS"])
        res["lead_z_GWAS"] = float(lead["z_GWAS"])
        return res

    all_res = [run_one(None, "full_4Mb")]
    for kb in [1000, 500, 250]:
        all_res.append(run_one(kb * 1000, f"pm{kb}kb"))
    all_res = [r for r in all_res if r]
    if not all_res:
        return None
    out_json = OUT / f"coloc_{name.lower()}_eqtlgen_bcac.json"
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(all_res, fh, indent=2, ensure_ascii=False)
    print(f"[{name}] saved {out_json}")
    return all_res


def main():
    summary_rows = []
    for name, info in GENES.items():
        res = run_gene(name, info["ensg"], info["bp"])
        if res:
            for r in res:
                summary_rows.append({
                    k: r.get(k) for k in [
                        "gene", "window", "n_snps", "PP.H0.abf", "PP.H1.abf",
                        "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "lead_snp",
                        "lead_snp_pp_h4", "lead_z_eQTL", "lead_z_GWAS",
                        "region_chr", "region_start", "region_end",
                    ]
                })
    if summary_rows:
        df = pd.DataFrame(summary_rows)
        df.to_csv(OUT / "coloc_bonus_summary.csv", index=False,
                  encoding="utf-8-sig")
        print("\n=== SUMMARY ===")
        print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
