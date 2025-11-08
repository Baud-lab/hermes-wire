#!/usr/bin/env Rscript
# collect_fullmed_genes_fold.R
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--med_tsv",   type="character", help="med__GENUS__foldXX.tsv (must have columns 'conclusion' and 'gene_name')"),
  make_option("--med_rdata", type="character", help="mediation_raw_data__GENUS__foldXX.RData (must contain object DF)"),
  make_option("--out_genes", type="character", help="Output: fullmed_genes__GENUS__foldXX.txt"),
  make_option("--out_df",    type="character", help="Output: DF_fullmed__GENUS__foldXX.tsv or .tsv.gz"),
  make_option("--verbose",   action="store_true", default=TRUE)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$med_tsv), !is.null(opt$med_rdata), !is.null(opt$out_genes), !is.null(opt$out_df))

msg <- function(...) if (isTRUE(opt$verbose)) cat(sprintf(...), "\n")

# helper: write tsv (supports .gz)
write_tsv <- function(DT, path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open="wb")
    on.exit(close(con), add=TRUE)
    write.table(DT, file=con, sep="\t", quote=FALSE, row.names=FALSE)
  } else {
    fwrite(DT, file=path, sep="\t")
  }
}

# 1) Read mediation TSV → gene list with Full mediation
msg("[fullmed] Reading: %s", opt$med_tsv)
M <- fread(opt$med_tsv, sep="\t", na.strings=c("NA",""))
stopifnot("conclusion" %in% names(M), "gene_name" %in% names(M))
G <- unique(M[conclusion == "Full mediation", gene_name])

if (length(G) == 0L) {
  msg("[fullmed] No 'Full mediation' rows; writing empty outputs")
  writeLines(character(), con=opt$out_genes)
  write_tsv(data.table(), opt$out_df)
  quit(save="no", status=0)
}

writeLines(G, con=opt$out_genes)
msg("[fullmed] Genes with Full mediation: %d", length(G))

# 2) Load RData (expects object DF)
msg("[fullmed] Loading: %s", opt$med_rdata)
load(opt$med_rdata)  # must create object DF
stopifnot(exists("DF"))
stopifnot("id" %in% colnames(DF))

# gene columns present in DF
g_in_df <- intersect(G, colnames(DF))
if (length(g_in_df) == 0L) {
  msg("[fullmed] None of the Full mediation genes are columns in DF; writing empty DF")
  write_tsv(data.table(), opt$out_df)
  quit(save="no")
}

k <- ncol(DF)
stopifnot(k >= 7)  # id + >=1 gene + last 5
last5_idx <- (k-4):k
keep <- c("id", g_in_df, colnames(DF)[last5_idx])

DF2 <- as.data.table(DF[, keep, drop=FALSE])
write_tsv(DF2, opt$out_df)
msg("[fullmed] Wrote: %s (rows=%d, cols=%d)", opt$out_df, nrow(DF2), ncol(DF2))
