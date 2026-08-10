import sys
sys.path.insert(0, r"work")
from coloc_fdps import load_smr_plot, coloc_abf
import numpy as np, pandas as pd

snp_info, gwas, eqtl, meta = load_smr_plot(r"plot\fdps.ENSG00000160752.txt")
m = eqtl["ENSG00000160752"].merge(gwas, on="SNP", how="inner").merge(
    snp_info[["SNP", "Chr", "BP", "A1", "A2"]], on="SNP", how="left")
m["BP"] = pd.to_numeric(m["BP"])
m = m.dropna(subset=["b_eQTL", "se_eQTL", "b_GWAS", "se_GWAS", "BP"])
m = m[(m["se_eQTL"] > 0) & (m["se_GWAS"] > 0)].sort_values("BP")
print("n", len(m))

rows = []
for p12 in [1e-8, 1e-5, 1e-4]:
    for corr in [True, False]:
        r = coloc_abf(
            m["SNP"].values, m["b_eQTL"].values, m["se_eQTL"].values,
            m["b_GWAS"].values, m["se_GWAS"].values,
            n1=31684, n2=247173, s2=0.2929, corrections=corr,
            p1=1e-4, p2=1e-4, p12=p12)
        rows.append({
            "p12": p12, "variance_corrections": corr,
            "PP.H4.abf": r["PP.H4.abf"], "PP.H3.abf": r["PP.H3.abf"],
            "PP.H1.abf": r["PP.H1.abf"], "lead_snp": r["lead_snp"],
        })
        print("p12={:.0e} corr={}: PP.H4={:.6f} PP.H3={:.6f} lead={}".format(
            p12, corr, r["PP.H4.abf"], r["PP.H3.abf"], r["lead_snp"]))

df = pd.DataFrame(rows)
df.to_csv(r"outputs\chol_metab_signature\coloc_fdps_sensitivity.csv",
          index=False, encoding="utf-8-sig")
print("saved")
