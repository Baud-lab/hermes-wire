#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
})

# -------- CLI --------
option_list <- list(
  make_option("--metadata",  type="character",
              help="TSV with sample IDs (RFID or sample), dam, cage, etc."),
  make_option("--outdir",    type="character",
              help="Output directory", default="."),
  make_option("--keep_ids",  type="character",
              help="Optional file with authoritative sample IDs (one per line) to subset+order SRMs.", default=NULL)
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -------- I/O helpers --------
read_keep_ids <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  x <- readLines(path, warn = FALSE)
  x <- trimws(sub("\r$", "", x))
  x <- x[nzchar(x)]
  unique(x)
}

write_grm_txt <- function(M, ids, path) {
  # GCTA .grm.gz: i j N_ij value  (we use N_ij=1 for design/SRM covariances)
  con <- gzfile(path, "wt"); on.exit(close(con))
  n <- length(ids)
  for (i in seq_len(n)) {
    for (j in i:n) {
      val <- M[i, j]
      # Always write every upper-triangular entry (including zeros), as expected by GCTA
      writeLines(paste(i, j, 1, sprintf("%.6f", val)), con)
    }
  }
  idf <- data.frame(FID = ids, IID = ids, check.names = FALSE)
  write.table(idf, sub(".grm.gz$", ".grm.id", path),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# -------- Build SRM (symmetric sparse) --------
# 1 on diagonal; 1 for pairs sharing the same grouping label; 0 otherwise
make_srm_sparse <- function(ids, group) {
  stopifnot(length(ids) == length(group))
  n <- length(ids)
  # Start with identity (diagonal = 1)
  M <- Diagonal(n)  # dgCMatrix identity
  
  # Ensure factor to get stable levels, drop unused
  group <- droplevels(factor(group))
  
  # For each group, set all within-group pairs to 1
  levs <- levels(group)
  for (lv in levs) {
    idx <- which(group == lv)
    if (length(idx) > 1) {
      # Add a full block of 1s for this group
      # (diagonal entries may become 2; we'll cap to 1 below)
      M <- M + sparseMatrix(
        i = rep(idx, each = length(idx)),
        j = rep(idx, times = length(idx)),
        x = 1, dims = c(n, n)
      )
    }
    # if length(idx) == 1, the Diagonal already set M[idx, idx] = 1
  }
  
  # Cap any entries > 1 to 1, enforce symmetry, and set dimnames
  M@x[M@x > 1] <- 1
  M <- forceSymmetric(M, uplo = "U")
  M <- as(M, "dsCMatrix")  # symmetric sparse numeric
  rownames(M) <- colnames(M) <- ids
  M
}

# -------- Main --------
cat("[build_srm] Reading metadata:", opt$metadata, "\n")
meta <- fread(opt$metadata)

# Standardize ID column name to 'sample'
if ("RFID" %in% names(meta) && !"sample" %in% names(meta)) {
  setnames(meta, "RFID", "sample")
}

req <- c("sample", "dam", "cage")
if (!all(req %in% names(meta))) {
  stop("[build_srm] Metadata must contain columns: ", paste(req, collapse = ", "))
}

# Basic cleaning: drop rows with missing keys; coerce dam/cage to factor
meta <- meta[!is.na(sample) & !is.na(dam) & !is.na(cage)]
meta[, dam  := factor(dam)]
meta[, cage := factor(cage)]

# De-duplicate by sample (keep first), warn if duplicates were present
dups <- meta[duplicated(sample), unique(sample)]
if (length(dups)) {
  warning("[build_srm] Found duplicated sample IDs; keeping first occurrence for: ",
          paste(head(dups, 10), collapse = ", "),
          if (length(dups) > 10) sprintf(" ... (+%d more)", length(dups) - 10) else "")
  meta <- meta[!duplicated(sample)]
}

# Optional authoritative ordering/subsetting (e.g., SeqArray sample.id)
keep_ids <- read_keep_ids(opt$keep_ids)

if (!is.null(keep_ids)) {
  # Report overlap and losses
  in_both <- intersect(keep_ids, meta$sample)
  if (!length(in_both)) {
    stop("[build_srm] No overlap between keep_ids and metadata sample IDs.")
  }
  if (length(in_both) < length(keep_ids)) {
    missing <- setdiff(keep_ids, in_both)
    warning("[build_srm] ", length(missing), " keep_ids not found in metadata; e.g.: ",
            paste(head(missing, 10), collapse = ", "),
            if (length(missing) > 10) sprintf(" ... (+%d more)", length(missing) - 10) else "")
  }
  if (length(in_both) < nrow(meta)) {
    dropped <- setdiff(meta$sample, in_both)
    if (length(dropped)) {
      message("[build_srm] Dropping ", length(dropped), " metadata rows not in keep_ids.")
    }
  }
  # Subset and order metadata to match keep_ids
  meta <- meta[match(keep_ids, meta$sample), , nomatch = 0]
}

# Final IDs vector and quick summary
ids <- meta$sample
cat("[build_srm] N samples:", length(ids),
    " unique dams:", length(unique(meta$dam)),
    " unique cages:", length(unique(meta$cage)), "\n")

# Build SRMs
cat("[build_srm] Building SRM (dam)...\n")
srm_dam  <- make_srm_sparse(ids, meta$dam)

cat("[build_srm] Building SRM (cage)...\n")
srm_cage <- make_srm_sparse(ids, meta$cage)

# Save RDS (symmetric sparse) and the ID list
saveRDS(srm_dam,  file = file.path(opt$outdir, "srm_dam.rds"))
saveRDS(srm_cage, file = file.path(opt$outdir, "srm_cage.rds"))
writeLines(ids, file.path(opt$outdir, "srm_ids.txt"))

cat("[build_srm] SRM matrices saved. Writing GREML .grm[.gz]/.grm.id ...\n")
write_grm_txt(srm_dam,  ids, file.path(opt$outdir, "srm_dam.grm.gz"))
write_grm_txt(srm_cage, ids, file.path(opt$outdir, "srm_cage.grm.gz"))

cat("[build_srm] Done.\n")
