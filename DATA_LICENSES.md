# Data sources, accessions and terms of use

| Resource | Accession / source | Notes |
|---|---|---|
| TCGA-BRCA | NCI GDC portal | Public; cite TCGA Network 2012 and GDC |
| TCGA-CDR | Liu et al. 2018 | Public |
| METABRIC | cBioPortal DataHub | Public; cite Curtis et al. 2012 |
| PAM50 calls | UCSC Xena GDC hub / cBioPortal | Public |
| FDPS HM450 methylation | cBioPortal brca_tcga_methylation_hm450 | Public |
| GSE21653 | NCBI GEO (GPL570) | Public |
| GSE7390 | ArrayExpress E-GEOD-7390 (GPL96) | Public |
| GSE176078 | CELLxGENE (Wu et al. 2021) | Public |
| GSE161529 | Mendeley Data mirror of Pal et al. 2021 | Public |
| eQTLGen cis-eQTLs | eQTLGen consortium | Academic use; cite Vosa et al. 2021 |
| GTEx v8 | Yang Lab SMR database / GTEx Portal | GTEx data-use agreement; cite GTEx Consortium |
| BCAC overall / TNBC | GWAS Catalog GCST010098 / GCST010100 | GWAS Catalog terms; cite Zhang et al. 2020 |
| BCAC OncoArray ER+/ER- | GWAS Catalog GCST004988 | GWAS Catalog terms; cite Michailidou et al. 2017. The original BCAC release prohibits reposting to third-party sites; use the download script and cite the source |
| 1000 Genomes Phase 3 EUR | 1000 Genomes project | Public |
| LM22 signature matrix | Newman et al. 2015 (Nat Methods); mirror in `data/lm22.txt` from the public CIBERSORTx repository | Used for CIBERSORT deconvolution; cite Newman et al. 2015 |
| HPA FDPS entry | Human Protein Atlas (proteinatlas.org/ENSG00000160752) | Public; cite Uhlen et al. 2015 |
| cBioPortal pan-cancer FDPS RNA-seq | cBioPortal datahub | Public; cite Cerami et al. 2012 |

Raw data are not redistributed in this repository. Run `data/download_data.py`
to fetch them, and comply with each source's terms of use.
