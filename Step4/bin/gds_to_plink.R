#!/usr/bin/env Rscript
# Purpose: Convert GDS genotypes to PLINK .bed/.bim/.fam for downstream tools,
#          but first restrict to the cohort (scan_ids) and drop SNPs that do
#          NOT show all three dosages (0,1,2) among those samples.
#
# Notes:
# - Uses SNPRelate (fast GDS blocks) to count per-variant dosage presence.
# - Missing-coded genotypes (e.g., 3) are set to NA. Any non-integers are rounded.
# - Keeps a SNP iff it has at least 'min_dosage_levels' distinct dosages among
#   the cohort (default=3, i.e., must have 0, 1, and 2 represented).
# - Finally calls snpgdsGDS2BED with both sample.id and snp.id subsets.
#
# References:
# - SNPRelate genotype coding (0/1/2/3, where 3 = missing) and usage appear in the
#   official vignette and reference manual. See Bioconductor SNPRelate docs. [1][2]

suppressPackageStartupMessages({
  library(optparse)
  library(SNPRelate)
  library(gdsfmt)
})

# ---------------- CLI ----------------
option_list <- list(
  make_option("--gds",  type = "character", help = "Input GDS file (candidate SNPs)"),
  make_option("--out",  type = "character", help = "Output PLINK base path (no extension)"),
  make_option("--scan_ids", type = "character", default = NULL,
              help = "Optional: path to genesis_final_scan_ids.txt to restrict samples"),
  make_option("--min_dosage_levels", type = "integer", default = 3,
              help = "Distinct dosage levels (among 0,1,2) required to KEEP a SNP [default: %default]"),
  make_option("--block_size", type = "integer", default = 5000,
              help = "Variants per block when scanning genotypes [default: %default]"),
  make_option("--write_sample_map", action = "store_true", default = FALSE,
              help = "Also write <out>.sample_map.tsv with sample.id and optional sample.label")
)
opt <- parse_args(OptionParser(option_list = option_list))

gds_path <- opt$gds
out_base <- opt$out
if (is.null(gds_path) || is.null(out_base))
  stop("[gds_to_plink] You must provide --gds and --out")
if (!file.exists(gds_path))
  stop("[gds_to_plink] GDS not found: ", gds_path)

cat("[gds_to_plink] Opening GDS (read-only):", gds_path, "\n")
g <- snpgdsOpen(gds_path, readonly = TRUE)
on.exit({ try(snpgdsClose(g), silent = TRUE) }, add = TRUE)

# -------------- helpers --------------
has_node <- function(g, node) {
  out <- try(index.gdsn(g, node), silent = TRUE)
  !inherits(out, "try-error")
}
require_nodes <- function(g, nodes) {
  missing <- nodes[!vapply(nodes, function(n) has_node(g, n), logical(1))]
  if (length(missing)) {
    stop("[gds_to_plink] Missing required GDS node(s): ", paste(missing, collapse = ", "))
  }
}

# -------------- validate structure --------------
required <- c("genotype", "sample.id", "snp.id", "snp.chromosome", "snp.position")
require_nodes(g, required)

sample_id <- read.gdsn(index.gdsn(g, "sample.id"))
snp_id    <- read.gdsn(index.gdsn(g, "snp.id"))
n_samp    <- length(sample_id)
n_snp     <- length(snp_id)
if (n_samp == 0L || n_snp == 0L)
  stop("[gds_to_plink] Empty sample or SNP set in GDS (samples=", n_samp, ", snps=", n_snp, ")")

has_allele  <- has_node(g, "snp.allele")
has_label   <- has_node(g, "sample.label")

geno_node <- index.gdsn(g, "genotype")
geno_dim  <- objdesp.gdsn(geno_node)$dim
if (length(geno_dim) != 2L)
  stop("[gds_to_plink] 'genotype' node is not a 2D array; got dim=", paste(geno_dim, collapse = "x"))

snpfirstdim <- NA
if (geno_dim[1] == n_snp && geno_dim[2] == n_samp) {
  snpfirstdim <- TRUE
} else if (geno_dim[1] == n_samp && geno_dim[2] == n_snp) {
  snpfirstdim <- FALSE
} else {
  stop(sprintf("[gds_to_plink] Genotype dims (%d x %d) inconsistent with n_snp=%d, n_samp=%d",
               geno_dim[1], geno_dim[2], n_snp, n_samp))
}
cat(sprintf("[gds_to_plink] Detected genotype dims: %d x %d  => snpfirstdim=%s\n",
            geno_dim[1], geno_dim[2], ifelse(snpfirstdim, "TRUE", "FALSE")))
cat(sprintf("[gds_to_plink] Samples: %d  SNPs: %d\n", n_samp, n_snp))

# -------------- Sample restriction (cohort) --------------
sample_keep <- sample_id
if (!is.null(opt$scan_ids)) {
  if (!file.exists(opt$scan_ids))
    stop("[gds_to_plink] --scan_ids file not found: ", opt$scan_ids)
  keep_scan <- scan(opt$scan_ids, what = character(), quiet = TRUE)
  # Coerce to the same type as sample_id for matching
  if (is.numeric(sample_id)) {
    keep_scan <- as.numeric(keep_scan)
  } else {
    keep_scan <- as.character(keep_scan)
  }
  in_gds <- keep_scan[keep_scan %in% sample_id]
  if (length(in_gds) == 0L)
    stop("[gds_to_plink] None of the provided scan_ids exist in GDS sample.id")
  sample_keep <- in_gds
  cat("[gds_to_plink] Restricting to cohort samples: ", length(sample_keep), " of ", n_samp, "\n", sep = "")
}

# -------------- NEW: Cohort-aware variant filtering before write --------------
cat("[gds_to_plink] Scanning variants to keep only those with at least ",
    opt$min_dosage_levels, " distinct dosages (among 0/1/2) in the cohort...\n", sep = "")

bsz <- max(1000L, as.integer(opt$block_size))
ns  <- length(snp_id)
blocks <- split(seq_len(ns), ceiling(seq_len(ns) / bsz))

keep_snp <- rep(FALSE, ns)
for (bi in seq_along(blocks)) {
  idx <- blocks[[bi]]
  geno <- snpgdsGetGeno(g,
                        sample.id   = sample_keep,
                        snp.id      = snp_id[idx],
                        snpfirstdim = snpfirstdim,
                        with.id     = FALSE,
                        verbose     = FALSE)
  # snpgdsGetGeno returns matrix [variant x sample] if snpfirstdim==TRUE, else [sample x variant]
  if (!snpfirstdim) geno <- t(geno)
  # Recode missing (3) -> NA; round any non-integers
  geno[geno > 2] <- NA
  mode(geno) <- "numeric"
  geno <- round(geno)
  # Count presence of 0/1/2 per row
  has0 <- rowSums(geno == 0, na.rm = TRUE) > 0
  has1 <- rowSums(geno == 1, na.rm = TRUE) > 0
  has2 <- rowSums(geno == 2, na.rm = TRUE) > 0
  keep_vec <- (has0 + has1 + has2) >= opt$min_dosage_levels
  keep_snp[idx] <- keep_vec
  if ((bi %% 20) == 0 || bi == length(blocks)) {
    cat(sprintf("[gds_to_plink]   block %d/%d ... kept so far: %d\n",
                bi, length(blocks), sum(keep_snp)))
  }
}
snp_keep_ids <- snp_id[keep_snp]
cat("[gds_to_plink] Variants kept after dosage filter: ", length(snp_keep_ids), " / ", ns, "\n", sep = "")
if (length(snp_keep_ids) == 0L) {
  stop("[gds_to_plink] No variants passed the dosage diversity filter among the cohort.")
}

# -------------- conversion --------------
cat("[gds_to_plink] Writing PLINK trio to base:", out_base, "\n")
snpgdsGDS2BED(
  gdsobj      = g,
  bed.fn      = out_base,
  sample.id   = sample_keep,      # cohort
  snp.id      = snp_keep_ids,     # filtered variants
  snpfirstdim = snpfirstdim,
  verbose     = TRUE
)

# -------------- post-write sanity --------------
bed <- paste0(out_base, ".bed")
bim <- paste0(out_base, ".bim")
fam <- paste0(out_base, ".fam")
out_files <- c(bed, bim, fam)
if (!all(file.exists(out_files))) {
  stop("[gds_to_plink] Expected PLINK files not found: ",
       paste(out_files[!file.exists(out_files)], collapse = ", "))
}
fi <- file.info(out_files)
cat("[gds_to_plink] File sizes (bytes):\n")
print(fi[, "size", drop = FALSE])

count_lines <- function(path) {
  con <- file(path, open = "r"); on.exit(close(con), add = TRUE)
  n <- 0L
  repeat { ln <- readLines(con, n = 100000L, warn = FALSE); if (!length(ln)) break; n <- n + length(ln) }
  n
}
bim_n <- count_lines(bim)
fam_n <- count_lines(fam)
cat(sprintf("[gds_to_plink] Row counts:  .bim = %d  .fam = %d  (after filtering)\n", bim_n, fam_n))

# -------------- optional: sample map --------------
if (isTRUE(opt$write_sample_map)) {
  out_map <- paste0(out_base, ".sample_map.tsv")
  if (has_label) {
    sample_label <- read.gdsn(index.gdsn(g, "sample.label"))
    df <- data.frame(sample.id = sample_id,
                     sample.label = if (length(sample_label) == length(sample_id)) sample_label else NA_character_,
                     stringsAsFactors = FALSE)
  } else {
    df <- data.frame(sample.id = sample_id, stringsAsFactors = FALSE)
  }
  write.table(df, file = out_map, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[gds_to_plink] Wrote sample map:", out_map, "\n")
}

cat("[gds_to_plink] Done. Created files:",
    paste(out_files, collapse = " "), "\n")
