#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FDPS locus colocalisation (coloc.abf-style, implemented from the published
Giambartolomei et al. 2014 model using Wakefield approximate Bayes factors).

Data:
  eQTL side  : eQTLGen whole-blood cis-eQTLs (SMR-format BESD)
  GWAS side  : BCAC overall breast cancer GWAS
  LD/region  : SMR --plot output for ENSG00000160752 (FDPS) +/- 2 Mb

The per-SNP Wakefield ABF is
  ABF = sqrt(v/(v+w)) * exp( w*beta^2 / (2*v*(v+w)) )
with v = se^2 and w = 0.2 (prior variance, coloc default).

Posterior probabilities are computed from first principles with priors
p1 = p2 = 1e-4, p12 = 1e-5 (coloc.abf defaults / recommended settings).
"""

import sys
import numpy as np
import pandas as pd


def logsumexp(x):
    x = np.asarray(x, dtype=float)
    m = np.max(x)
    if not np.isfinite(m):
        return m
    return m + np.log(np.sum(np.exp(x - m)))


def wakefield_abf(beta, se, w=0.2):
    """Return log(ABF) for a single SNP (Wakefield 2007)."""
    v = np.asarray(se, dtype=float) ** 2
    b = np.asarray(beta, dtype=float)
    v = np.maximum(v, 1e-300)
    l = 0.5 * (np.log(v) - np.log(v + w)) + (w * b * b) / (2.0 * v * (v + w))
    return l


def coloc_abf(rsid, beta1, se1, beta2, se2,
              n1, n2, type1="quant", type2="cc",
              s2=None, sdY1=1.0, corrections=True,
              p1=1e-4, p2=1e-4, p12=1e-5, w=0.2):
    """
    Bayesian colocalisation for one genomic region.

    Returns dict with PP.H0-H4, SNP-level PP.H4, and the leading SNP.
    """
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
        raise ValueError("fewer than 2 usable SNPs in region")

    # The R coloc.abf implementation adds a small variance term to account
    # for uncertainty in the phenotype variance:
    #   quant: varbeta + sdY^2/(2*N)      (sdY defaults to 1)
    #   cc   : varbeta + 1/(2*N*s*(1-s))
    # We reproduce this to stay faithful to the reference implementation.
    if corrections:
        if type1 == "quant" and n1 is not None:
            se1 = np.sqrt(se1**2 + sdY1**2 / (2.0 * n1))
        if type2 == "cc" and s2 is not None and n2 is not None:
            se2 = np.sqrt(se2**2 + 1.0 / (2.0 * n2 * s2 * (1.0 - s2)))

    l1 = wakefield_abf(beta1, se1, w=w)
    l2 = wakefield_abf(beta2, se2, w=w)

    # exact log-likelihoods of the five hypotheses (per-SNP priors inside)
    lH0 = 0.0
    lH1 = logsumexp(l1) + np.log(p1) - np.log(n)
    lH2 = logsumexp(l2) + np.log(p2) - np.log(n)
    # H3: distinct causal SNPs -> sum over j != k ABF1j*ABF2k
    sum1 = np.exp(logsumexp(l1))
    sum2 = np.exp(logsumexp(l2))
    sum12 = np.exp(logsumexp(l1 + l2))
    h3sum = sum1 * sum2 - sum12
    lH3 = np.log(max(h3sum, np.finfo(float).tiny)) + np.log(p1) + np.log(p2) - 2.0 * np.log(n)
    # H4: shared causal SNP
    lH4 = logsumexp(l1 + l2) + np.log(p12) - np.log(n)

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
        "n1": n1,
        "n2": n2,
        "s2": s2,
        "p1": p1,
        "p2": p2,
        "p12": p12,
        "w": w,
    }


def load_smr_plot(path):
    """Parse the SMR --plot text format.

    Layout:
      $probe <n> <probeID>
      <n probe summary lines>
      $SNP <nsnp>
      <nsnp lines: SNP chr bp A1 A2>
      $GWAS <ngwas>
      <ngwas lines: SNP b_GWAS se_GWAS>
      $eQTL <nprobe>
      <per probe: 'probe n_eqtl' then n_eqtl lines 'SNP b_eQTL se_eQTL p_eQTL'>

    Returns (snp_info, gwas_df, eqtl_by_probe) where each is a DataFrame.
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rows.append(line.rstrip("\n").split())

    if not rows:
        raise ValueError("empty plot file")
    assert rows[0][0] == "$probe"
    n_probe = int(rows[0][1])

    # find section markers
    idx_snp = next(i for i, r in enumerate(rows) if r and r[0] == "$SNP")
    idx_gwas = next(i for i, r in enumerate(rows) if r and r[0] == "$GWAS")
    idx_eqtl = next(i for i, r in enumerate(rows) if r and r[0] == "$eQTL")

    n_snp = int(rows[idx_snp][1])
    n_gwas = int(rows[idx_gwas][1])
    n_probe_eqtl = int(rows[idx_eqtl][1])
    assert idx_snp + 1 + n_snp == idx_gwas
    assert idx_gwas + 1 + n_gwas == idx_eqtl

    snp_info = pd.DataFrame(
        rows[idx_snp + 1 : idx_snp + 1 + n_snp],
        columns=["SNP", "Chr", "BP", "A1", "A2"],
    )
    gwas = pd.DataFrame(
        rows[idx_gwas + 1 : idx_gwas + 1 + n_gwas],
        columns=["SNP", "b_GWAS", "se_GWAS"],
    )
    gwas["b_GWAS"] = pd.to_numeric(gwas["b_GWAS"], errors="coerce")
    gwas["se_GWAS"] = pd.to_numeric(gwas["se_GWAS"], errors="coerce")

    eqtl_by_probe = {}
    i = idx_eqtl + 1
    for _ in range(n_probe_eqtl):
        probe, n_str = rows[i]
        n = int(n_str)
        block = pd.DataFrame(
            rows[i + 1 : i + 1 + n],
            columns=["SNP", "b_eQTL", "se_eQTL", "p_eQTL"],
        )
        for c in ["b_eQTL", "se_eQTL", "p_eQTL"]:
            block[c] = pd.to_numeric(block[c], errors="coerce")
        eqtl_by_probe[probe] = block
        i += 1 + n
    assert i == len(rows), (i, len(rows))
    return snp_info, gwas, eqtl_by_probe, {"n_probe": n_probe, "n_snp": n_snp,
                                          "n_gwas": n_gwas,
                                          "n_probe_eqtl": n_probe_eqtl}


def run():
    import json

    plot = sys.argv[1] if len(sys.argv) > 1 else (
        r"plot\fdps.ENSG00000160752.txt")
    out_json = sys.argv[2] if len(sys.argv) > 2 else (
        r"outputs\chol_metab_signature\coloc_fdps_eqtlgen_bcac.json")

    snp_info, gwas, eqtl_by_probe, meta = load_smr_plot(plot)
    print("meta:", meta)
    probe = "ENSG00000160752"
    if probe not in eqtl_by_probe:
        raise KeyError(f"{probe} not in eQTL blocks; available: "
                       f"{list(eqtl_by_probe)[:5]}...")
    eqtl = eqtl_by_probe[probe]
    print(f"GWAS SNPs: {len(gwas)}, FDPS eQTL SNPs: {len(eqtl)}")

    # merge on SNP ID (SMR plot output is allele-aligned to A1 for both)
    m = eqtl.merge(gwas, on="SNP", how="inner")
    m = m.merge(snp_info[["SNP", "Chr", "BP", "A1", "A2"]], on="SNP", how="left")
    m["BP"] = pd.to_numeric(m["BP"], errors="coerce")
    m["Chr"] = m["Chr"].astype(str)
    m = m.dropna(subset=["b_eQTL", "se_eQTL", "b_GWAS", "se_GWAS"])
    m = m[(m["se_eQTL"] > 0) & (m["se_GWAS"] > 0)]
    m = m.dropna(subset=["BP"])
    m = m.sort_values("BP").reset_index(drop=True)
    print("merged usable SNPs:", len(m))
    print("region:", m["Chr"].iloc[0], int(m["BP"].min()), "-",
          int(m["BP"].max()))

    z1 = m["b_eQTL"] / m["se_eQTL"]
    z2 = m["b_GWAS"] / m["se_GWAS"]
    m["z_eQTL"] = z1
    m["z_GWAS"] = z2
    top_e = m.loc[m["z_eQTL"].abs().idxmax()]
    top_g = m.loc[m["z_GWAS"].abs().idxmax()]
    print("top eQTL SNP:", top_e["SNP"], f"z={top_e['z_eQTL']:.2f}",
          f"b={top_e['b_eQTL']:.4g} se={top_e['se_eQTL']:.4g}")
    print("top GWAS SNP:", top_g["SNP"], f"z={top_g['z_GWAS']:.2f}",
          f"b={top_g['b_GWAS']:.4g} se={top_g['se_GWAS']:.4g}")
    print(f"corr(z_eQTL, z_GWAS) = {np.corrcoef(z1, z2)[0, 1]:.3f}")

    def run_one(frac, label, out_path):
        if frac is None:
            sub = m
        else:
            probe_bp = 155284498
            sub = m[m["BP"].between(probe_bp - frac / 2.0,
                                    probe_bp + frac / 2.0)]
        if len(sub) < 2:
            print("skip", label, "fewer than 2 SNPs")
            return None
        res = coloc_abf(
            rsid=sub["SNP"].values,
            beta1=sub["b_eQTL"].values,
            se1=sub["se_eQTL"].values,
            beta2=sub["b_GWAS"].values,
            se2=sub["se_GWAS"].values,
            n1=31684,            # eQTLGen
            n2=247173,           # BCAC overall
            type1="quant",
            type2="cc",
            s2=0.2929,
            sdY1=1.0,
            corrections=True,
            p1=1e-4,
            p2=1e-4,
            p12=1e-5,
            w=0.2,
        )
        res["window"] = label
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
        print("\n=== coloc result:", label, "===")
        for k in ["n_snps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                  "PP.H3.abf", "PP.H4.abf", "lead_snp", "lead_snp_pp_h4"]:
            print(f"{k:16s} {res[k]:.6g}" if isinstance(res[k], float)
                  else f"{k:16s} {res[k]}")
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=2, ensure_ascii=False)
        print("saved", out_path)
        return res

    all_res = []
    # full cis window (+/- 2 Mb)
    all_res.append(run_one(None, "full_4Mb",
                           r"outputs\chol_metab_signature\coloc_fdps_full.json"))
    # sensitivity windows
    for kb in [1000, 500, 250]:
        all_res.append(run_one(
            kb * 1000, f"pm{kb}kb",
            rf"outputs\chol_metab_signature\coloc_fdps_{kb}kb.json"))

    summary = pd.DataFrame([r for r in all_res if r])
    if len(summary):
        cols = ["window", "n_snps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                "PP.H3.abf", "PP.H4.abf", "lead_snp", "lead_snp_pp_h4"]
        out_csv = r"outputs\chol_metab_signature\coloc_fdps_summary.csv"
        summary[cols].to_csv(out_csv, index=False, encoding="utf-8-sig")
        print("\nsaved", out_csv)

    # keep the old full-run output path for backward compatibility
    res = coloc_abf(
        rsid=m["SNP"].values,
        beta1=m["b_eQTL"].values,
        se1=m["se_eQTL"].values,
        beta2=m["b_GWAS"].values,
        se2=m["se_GWAS"].values,
        n1=31684,            # eQTLGen sample size (N = 31,684)
        n2=247173,           # BCAC overall cases+controls
        type1="quant",
        type2="cc",
        s2=0.2929,           # ~72,474 cases / 247,173
        sdY1=1.0,
        corrections=True,
        p1=1e-4,
        p2=1e-4,
        p12=1e-5,
        w=0.2,
    )
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print("saved", out_json)


if __name__ == "__main__":
    run()
