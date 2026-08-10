#!/usr/bin/env bash
# Step 5 (template): SMR/HEIDI for hub genes -> breast cancer risk.
# Run on a Linux/macOS machine WITH internet. All paths/URLs below are the
# official sources as of 2026-08-07 (verify before bulk download).

set -euo pipefail
OUT=outputs/chol_metab_signature/smr
mkdir -p "$OUT"
cd "$OUT"

# 1) SMR software (Linux x86_64)
# https://yanglab.westlake.edu.cn/software/smr/#Download
# wget https://yanglab.westlake.edu.cn/software/smr/download/SMR-1.3.1-linux-x86_64.tar.gz
# tar xzf SMR-1.3.1-linux-x86_64.tar.gz

# 2) 1000 Genomes reference (plink bfile, GRCh37) - required for LD
# e.g. from SMR site: https://yanglab.westlake.edu.cn/software/smr/#DataResource
# wget .../gwas.1000G.hs37d5.v1.tar.gz  (or use 1000G Phase3 plink files)

# 3) eQTL summary in SMR BESD format
# eQTLGen (blood, 31,684 individuals): https://www.eqtlgen.org/cis-eqtls.html
#   -> "Summary statistics in SMR format" (e.g. eqtlgen_smr.zip)
# GTEx v8 breast-mammary eQTLs (SMRdb): https://yanglab.westlake.edu.cn/data/SMRdb/
#   -> GTEx_V8_eQTL_SMR/breast_mammary_tissue.tar.bz2

# 4) BCAC breast cancer GWAS summary statistics
# Zhang H et al. Nat Genet 2020 (OncoArray + iCOGS + 11 GWAS):
#   https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-icogs-summary-results/
#   (or https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations)
# The provided file already has columns: SNP A1 A2 freq b se p n.
# GWAS Catalog alternative: GCST90132337 (breast cancer, overall, Zhang 2020).

BFILE=1000G.phase3.v5a/1000G.phase3.v5a  # adjust
EQTL=eqtlgen_smr/eQTLGen_cis_eQTLs_SMR  # BESD basename (eQTLGen) or GTEx breast BESD
GWAS=BCAC_meta/BCAC_Onco_iCOGS_meta_BC_overall_2020.txt
SMR=./smr-1.3.1-linux-x86_64/smr-1.3.1

HUB="ABCG1 DHCR24 DHCR7 FDXR G6PD HMGCS2 HSD17B7 LIMA1 NSDHL PRKAA1 VLDLR"

# 5) Run SMR + HEIDI for each hub gene (cis window 2000kb, HEIDI >0.05 = no LD confound)
for G in $HUB; do
  "$SMR" \
    --bfile "$BFILE" \
    --gwas-summary "$GWAS" \
    --beqtl-summary "$EQTL" \
    --gene-list "$G" \
    --out "$OUT/smr_${G}" \
    --thread-num 4 \
    --maf 0.01 \
    --peqtl 5e-8 \
    --heidi-mth 2 \
    --heidi-min-n 3
done

# 6) Results: smr_*.txt with columns: probe, topSNP, b_GWAS, se_GWAS, p_GWAS,
#    b_eQTL, se_eQTL, p_eQTL, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI.
# Interpretation: p_SMR < 0.05/n_genes (Bonferroni) + p_HEIDI > 0.05 -> causal gene.
# Optional sensitivity: also run with GTEx breast eQTL and with
# ER-negative / ER-positive BCAC subsets (available on the BCAC summary page).
echo "SMR template DONE - inspect $OUT/smr_*.txt"