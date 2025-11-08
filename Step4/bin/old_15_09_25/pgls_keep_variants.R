#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--corr", type="character", help="PGLS correlations_by_gene__{genus}.tsv"),
  make_option("--q",    type="double",   help="q-value threshold", default=0.10),
  make_option("--out",  type="character",help="Output TSV",        default="pgls_keep.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$corr) || !file.exists(opt$corr)) stop("Missing --corr or file not found")

DT <- fread(opt$corr, sep = "\t")
setDT(DT)
setnames(DT, tolower(names(DT)))

# Identify columns
q_col <- c("qvalue","q","q_used","qval")
p_col <- c("p_used","p","pvalue")
q_name <- q_col[q_col %in% names(DT)][1]
p_name <- p_col[p_col %in% names(DT)][1]

if (is.na(q_name)) {
  if (is.na(p_name)) stop("Neither q- nor p-value column found in PGLS file")
  # compute BH over all rows as a fallback
  DT[, qtmp := p.adjust(get(p_name), method = "BH")]
  q_name <- "qtmp"
}

# SNP identifier: prefer explicit SNP/rsid, else sentinel_snp, else chr:pos
snp_col <- c("snp","sentinel_snp","variant")
snp_name <- snp_col[snp_col %in% names(DT)][1]

# normalize chr/pos if needed
if (!("chr" %in% names(DT)) && "chromosome" %in% names(DT)) setnames(DT, "chromosome", "chr")
if (!("pos" %in% names(DT)) && "position"   %in% names(DT)) setnames(DT, "position",   "pos")

if (is.na(snp_name)) {
  if (!all(c("chr","pos") %in% names(DT))) {
    stop("Cannot build SNP ID: missing 'snp'/'sentinel_snp'/'variant' and also no 'chr'+'pos'.")
  }
  DT[, chr := sub("^chr", "", as.character(chr))]
  DT[, snp := paste0(chr, ":", as.integer(pos))]
} else {
  setnames(DT, snp_name, "snp", skip_absent = TRUE)
}

# Filter by q
KEEP <- DT[ is.finite(get(q_name)) & get(q_name) <= opt$q,
            .(snp,
              chr = if ("chr" %in% names(DT)) chr else NA_character_,
              pos = if ("pos" %in% names(DT)) pos else NA_real_,
              p   = if (!is.na(p_name)) get(p_name) else NA_real_,
              q   = get(q_name)) ]

# Deduplicate
KEEP <- unique(KEEP, by = "snp")

if (nrow(KEEP) == 0L) {
  warning("No SNPs pass PGLS q ≤ ", opt$q, " in: ", opt$corr)
}

fwrite(KEEP, file = opt$out, sep = "\t")
cat("[pgls-keep] wrote", opt$out, "rows=", nrow(KEEP), "\n")
