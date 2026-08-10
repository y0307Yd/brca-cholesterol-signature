args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else "outputs/chol_metab_signature/susie_fdps"

snp <- read.csv(file.path(dir, "snp_table.csv"), stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
names(snp)[1] <- "SNP"

read_bin <- function(fn) {
  con <- file(fn, "rb")
  on.exit(close(con))
  dims <- readBin(con, integer(), n = 2, size = 4)
  m <- matrix(readBin(con, double(), n = prod(dims), size = 8),
              nrow = dims[1], byrow = TRUE)
  m
}
R <- read_bin(file.path(dir, "ld_matrix.bin"))
Z <- read_bin(file.path(dir, "z_matrix.bin"))
stopifnot(nrow(Z) == nrow(snp), nrow(R) == nrow(snp), ncol(R) == nrow(snp))
colnames(R) <- rownames(R) <- snp$SNP
rownames(Z) <- snp$SNP
colnames(Z) <- c("eQTL", "GWAS")
cat("SNPs:", nrow(snp), "\n")

library(susieR)
library(coloc)

set.seed(20260810)
fit1 <- susie_rss(z = Z[, 1], R = R, L = 10, n = 31684,
                  estimate_residual_variance = TRUE,
                  coverage = 0.95)
fit2 <- susie_rss(z = Z[, 2], R = R, L = 10, n = 228951,
                  estimate_residual_variance = TRUE,
                  coverage = 0.95)

cat("eQTL credible sets:", if (is.null(fit1$sets)) 0 else length(fit1$sets$cs), "\n")
cat("GWAS credible sets:", if (is.null(fit2$sets)) 0 else length(fit2$sets$cs), "\n")

cs1 <- fit1$sets
cs2 <- fit2$sets
if (!is.null(cs1$cs) && !is.null(cs2$cs) && length(cs1$cs) > 0 && length(cs2$cs) > 0) {
  idx1 <- cs1$cs_index
  idx2 <- cs2$cs_index
  bf1 <- fit1$lbf_variable[idx1, , drop = FALSE]
  bf2 <- fit2$lbf_variable[idx2, , drop = FALSE]
  colnames(bf1) <- snp$SNP
  colnames(bf2) <- snp$SNP
  ret <- coloc:::coloc.bf_bf(bf1, bf2)
  saveRDS(list(fit1 = fit1, fit2 = fit2, bf1 = bf1, bf2 = bf2, ret = ret),
          file.path(dir, "debug_susie.rds"))
  cat("bf1 dim:", dim(bf1), "bf2 dim:", dim(bf2), "\n")
  cat("bf colnames overlap:", length(intersect(colnames(bf1), colnames(bf2))), "\n")
  str(ret)
  summ <- as.data.frame(ret$summary)
  summ$idx1 <- cs1$cs_index[idx1]
  summ$idx2 <- cs2$cs_index[idx2]
  write.csv(summ, file.path(dir, "coloc_susie_summary.csv"), row.names = FALSE)
  cat("coloc.susie PP.H4:", summ$PP.H4, "\n")
  print(summ)
} else {
  cat("No credible sets in one or both fits; coloc.susie skipped.\n")
}

out <- data.frame(SNP = snp$SNP,
                  PIP_eQTL = fit1$pip,
                  PIP_GWAS = fit2$pip,
                  stringsAsFactors = FALSE)
write.csv(out, file.path(dir, "susie_pip.csv"), row.names = FALSE)

cat("DONE\n")
