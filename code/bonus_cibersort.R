# CIBERSORT-style deconvolution (nu-SVR, linear kernel) with LM22.
# Core algorithm faithful to Newman et al. 2015 / CIBERSORT.R CoreAlg.
suppressMessages(library(e1071))

work <- "C:/Users/Y/Documents/Codex/2026-08-06/new-chat/work"

core_alg <- function(X, y) {
  set.seed(1)
  best <- NULL
  best_rss <- Inf
  for (nu in c(0.25, 0.5, 0.75)) {
    svr <- svm(X, y, scale = TRUE, type = "nu-regression",
               kernel = "linear", nu = nu, cost = 1)
    rss <- sum(svr$residuals^2)
    if (rss < best_rss) {
      best_rss <- rss
      best <- svr
    }
  }
  w <- t(best$coefs) %*% best$SV        # 1 x n_genes weight vector
  w <- as.numeric(w)
  w[w < 0] <- 0
  if (sum(w) > 0) w <- w / sum(w)
  list(frac = w, rss = best_rss)
}

lm22 <- read.delim("C:/Users/Y/Documents/Codex/2026-08-06/new-chat/data/lm22.txt",
                   row.names = 1, check.names = FALSE)
# plain expression matrix: gene symbols in rows, samples in columns
mixture <- read.delim(file.path(work, "estimate_input.txt"),
                      row.names = 1, check.names = FALSE)
genes <- intersect(rownames(lm22), rownames(mixture))
cat("common LM22 genes:", length(genes), "\n")
X <- as.matrix(lm22[genes, ])
Mi <- as.matrix(mixture[genes, ])

celltypes <- colnames(X)
frac <- matrix(NA_real_, ncol(Mi), length(celltypes),
               dimnames = list(colnames(Mi), celltypes))
rss_vec <- numeric(ncol(Mi))
for (i in seq_len(ncol(Mi))) {
  out <- core_alg(X, Mi[, i])
  frac[i, ] <- out$frac
  rss_vec[i] <- out$rss
  if (i %% 200 == 0) cat("sample", i, "done\n")
}

res <- data.frame(sample_id = rownames(frac), frac, stringsAsFactors = FALSE)
res$RMSE <- sqrt(rss_vec / nrow(X))
write.csv(res, file.path(work, "bonus_cibersort_fractions.csv"),
          row.names = FALSE)
cat("saved fractions; mean RMSE:", mean(res$RMSE), "\n")
