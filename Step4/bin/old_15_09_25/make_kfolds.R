#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opt_list <- list(
  make_option("--scan_ids", type="character", help="genesis_final_scan_ids.txt"),
  make_option("--k",        type="integer",  default=5, help="Number of folds"),
  make_option("--seed",     type="integer",  default=17),
  make_option("--outdir",   type="character", default="kfold_ids")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$scan_ids))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

ids <- fread(opt$scan_ids, header = FALSE)$V1
ids <- unique(ids)
set.seed(opt$seed)
ids <- sample(ids, length(ids), replace = FALSE)

# split into ~equal K folds
fold_id <- ceiling(seq_along(ids) / (length(ids)/opt$k))
fold_id[fold_id > opt$k] <- opt$k

for (f in seq_len(opt$k)) {
  test_ids  <- ids[ fold_id == f ]
  train_ids <- ids[ fold_id != f ]
  d <- file.path(opt$outdir, sprintf("fold%02d", f))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  fwrite(data.table(id = train_ids), file.path(d, "scan_ids_train.txt"), col.names = FALSE)
  fwrite(data.table(id = test_ids ), file.path(d, "scan_ids_test.txt" ), col.names = FALSE)
}
cat(sprintf("[kfold] Wrote %d folds under %s\n", opt$k, opt$outdir))

