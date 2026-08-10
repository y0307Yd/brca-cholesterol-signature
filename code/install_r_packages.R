# Install the R/Bioconductor packages used by the pipeline (R >= 4.6).
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c(
  "org.Hs.eg.db",   # cholesterol gene set annotation
  "hgu133a.db"      # GSE7390 (U133A) probe mapping
))
install.packages(c("susieR", "coloc"), repos = "https://cloud.r-project.org")

# Versions used in the manuscript analyses:
#   susieR 0.14.2, coloc 5.2.3, org.Hs.eg.db 3.23.1,
#   hgu133a.db 3.13.0, BiocManager 1.30.27
cat("Installation complete.\n")
