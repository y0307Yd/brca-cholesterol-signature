"""Download the public data used by the pipeline.

Usage:
  python download_data.py --dir /path/to/data

Large files (eQTLGen BESD, 1000 Genomes, OncoArray) are multi-GB; run each
section separately if needed. Verify checksums where provided.
"""
import argparse
import subprocess
from pathlib import Path

URLS = {
    # NCBI GEO
    "GSE21653_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21653/matrix/GSE21653_series_matrix.txt.gz",
    # ArrayExpress mirror of GSE7390 (GPL570 probe annotation is embedded in
    # the GSE21653 series-matrix file; GSE7390 probes are mapped with
    # hgu133a.db, installed from Bioconductor)
    "GSE7390_processed": "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-7390/E-GEOD-7390.processed.1.zip",
    "GSE7390_sdrf": "https://www.ebi.ac.uk/biostudies/files/E-GEOD-7390/E-GEOD-7390.sdrf.txt",
    # cBioPortal
    "cbioportal_methylation": "https://www.cbioportal.org/api/molecular-profiles/brca_tcga_methylation_hm450/molecular-data?molecularProfileId=brca_tcga_methylation_hm450&sampleListId=brca_tcga_all&entrezGeneId=2224",
    # eQTLGen SMR-format BESD (~1 GB)
    "eqtlgen_besd": "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz",
    # GTEx v8 SMR-format BESD (Yang Lab SMR database)
    "gtex_blood": "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary/Whole_Blood.zip",
    "gtex_breast": "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary/Breast_Mammary_Tissue.zip",
    # BCAC GWAS Catalog harmonised summary statistics
    "bcac_overall": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010098/harmonised/GCST010098.h.tsv.gz",
    "bcac_tnbc": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010100/harmonised/GCST010100.h.tsv.gz",
    "bcac_oncoarray": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
    # 1000 Genomes EUR (PLINK bfiles) - see SMR data resource page for the
    # canonical URL; the files used here are the g1000_eur reference.
    "g1000_eur": "https://yanglab.westlake.edu.cn/software/smr/download/1000G.phase3.v5a.tar.gz",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="data")
    ap.add_argument("--only", default=None,
                    help="comma-separated keys from URLS")
    args = ap.parse_args()
    out = Path(args.dir)
    out.mkdir(parents=True, exist_ok=True)
    keys = URLS if not args.only else [k.strip() for k in args.only.split(",")]
    for k in keys:
        url = URLS[k]
        dest = out / (k + (".gz" if url.endswith(".gz") else
                           (".zip" if ".zip" in url else ".txt.gz")))
        print(f"[{k}] {url}")
        subprocess.run(["curl", "-L", "-C", "-", "-o", str(dest), url], check=False)
        print("  ->", dest)


if __name__ == "__main__":
    main()
