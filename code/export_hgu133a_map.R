suppressMessages(library(hgu133a.db))
suppressMessages(library(AnnotationDbi))
m <- select(hgu133a.db, keys = keys(hgu133a.db, keytype = "PROBEID"),
            columns = "SYMBOL", keytype = "PROBEID")
m <- m[!is.na(m$SYMBOL), ]
write.csv(m, "outputs/chol_metab_signature/hgu133a_probe_map.csv",
          row.names = FALSE)
cat("hgu133a map rows:", nrow(m), "\n")
