#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

opt_list <- list(
  make_option("--sentinels", type="character", help="Sentinels TSV"),
  make_option("--assoc_species_dir", type="character", help="assoc_species_full dir (has by_trait)"),
  make_option("--bim", type="character", help="PLINK .bim (unused, kept for compatibility)"),
  make_option("--genera", type="character", help="Comma-separated genera list (optional)"),
  make_option("--taxonomy", type="character", help="Taxonomy file (unused here)"),
  make_option("--out", type="character", default="beta_cf_by_species__from_genus.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$sentinels), file.exists(opt$sentinels))
stopifnot(!is.null(opt$assoc_species_dir), dir.exists(opt$assoc_species_dir))
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

logf <- function(...) message(sprintf("[beta-cf] %s", sprintf(...)))

tolower_names <- function(dt) setnames(dt, names(dt), tolower(names(dt)))
pick_col <- function(nms, candidates) {
  x <- intersect(tolower(candidates), nms)
  if (length(x)) x[1] else NA_character_
}

## ---------- read sentinels ----------
sent <- fread(opt$sentinels, sep="\t", header=TRUE)
stopifnot(nrow(sent) > 0)
tolower_names(sent)

# normalize expected columns; prefer 'variant_id' if present
vid_col <- pick_col(names(sent), c("variant_id","variant.id","vid","index","row"))
snp_col <- pick_col(names(sent), c("snp","rsid","marker","id"))
chr_col <- pick_col(names(sent), c("chr","chrom","chromosome"))
bp_col  <- pick_col(names(sent), c("bp","pos","position"))
p_col   <- pick_col(names(sent), c("p","pval","p.value","pvalue"))

if (!is.na(vid_col)) setnames(sent, vid_col, "variant_id")
if (!is.na(snp_col)) setnames(sent, snp_col, "snp")
if (!is.na(chr_col)) setnames(sent, chr_col, "chr")
if (!is.na(bp_col))  setnames(sent,  bp_col, "pos")
if (!is.na(p_col))   setnames(sent,   p_col, "p")

# Derive chr/pos from SNP if needed
if (!"chr" %in% names(sent) || !"pos" %in% names(sent)) {
  if ("snp" %in% names(sent) && all(grepl("^\\w+:\\d+$", sent$snp))) {
    parts <- tstrsplit(sent$snp, ":", fixed=TRUE)
    sent[, chr := parts[[1]]]
    sent[, pos := as.integer(parts[[2]])]
  }
}

# Keep only usable columns
keep <- intersect(c("gene","genus","variant_id","snp","p","chr","pos"), names(sent))
sent <- unique(sent[, ..keep])[]
if (!("genus" %in% names(sent))) {
  stop("Sentinels file must contain a 'genus' column for strict scoping.")
}
sent[, genus := trimws(as.character(genus))]
logf("Sentinels: %d rows; genera present: %s",
     nrow(sent), paste(sort(unique(sent$genus)), collapse = ","))

## ---------- list assoc files ----------
assoc_dir <- file.path(opt$assoc_species_dir, "by_trait")
assoc_files <- list.files(assoc_dir, pattern = "_assoc\\.tsv(\\.gz)?$", full.names = TRUE)
stopifnot(length(assoc_files) > 0)

species_from_file <- function(f) {
  b <- basename(f)
  b <- sub("\\.tsv(\\.gz)?$", "", b, ignore.case=TRUE)
  sub("_assoc$", "", b, ignore.case=TRUE)
}
extract_genus <- function(sp) sub("^s__([^_]+).*", "\\1", sp)

## ---------- iterate ----------
gset <- if (!is.null(opt$genera) && nzchar(opt$genera)) trimws(strsplit(opt$genera, ",", fixed=TRUE)[[1]]) else character(0)
added_vid <- 0L
added_pos <- 0L
skipped_by_genus <- 0L

out_list <- vector("list", length(assoc_files))
for (i in seq_along(assoc_files)) {
  f  <- assoc_files[i]
  sp <- species_from_file(f)
  sp=sp[!grepl("s__Bacteroides_F_",sp)]
  sp_genus <- extract_genus(sp)
  
  # hard guard: only use sentinel rows that match THIS species file’s genus
  sent_g <- sent[genus == sp_genus]
  if (length(gset)) sent_g <- sent_g[genus %in% gset]
  if (!nrow(sent_g)) { skipped_by_genus <- skipped_by_genus + 1L; next }
  
  x <- tryCatch(fread(f, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE),
                error=function(e) NULL)
  if (is.null(x) || !nrow(x)) next
  tolower_names(x)
  
  # Identify standard columns in assoc file
  vid_x  <- pick_col(names(x), c("variant_id","variant.id","index","row"))
  chr_x  <- pick_col(names(x), c("chr","chrom","chromosome"))
  pos_x  <- pick_col(names(x), c("pos","bp","position"))
  beta_x <- pick_col(names(x), c("beta","est","estimate","beta_hat","effect"))
  se_x   <- pick_col(names(x), c("se_beta","se","est.se","beta_se","se_effect"))
  p_x    <- pick_col(names(x), c("p","pval","p.value","pvalue","score.pval"))
  
  mm <- NULL
  
  # (A) Join by variant_id if available
  if (!is.na(vid_x) && ("variant_id" %in% names(sent_g))) {
    setnames(x, vid_x, "variant_id", skip_absent = TRUE)
    if (!is.na(chr_x)) setnames(x, chr_x, "chr",  skip_absent = TRUE)
    if (!is.na(pos_x)) setnames(x, pos_x, "pos",  skip_absent = TRUE)
    mm <- merge(x, sent_g, by="variant_id", all.y=FALSE, all.x=FALSE)
    if (nrow(mm)) added_vid <- added_vid + nrow(mm)
  }
  
  # (B) Fallback: join by CHR:POS restricted to this genus’ sentinels
  if (is.null(mm) || !nrow(mm)) {
    if (!is.na(chr_x) && !is.na(pos_x)) {
      setnames(x, chr_x, "chr", skip_absent = TRUE)
      setnames(x, pos_x, "pos", skip_absent = TRUE)
      mm <- merge(x, sent_g[!is.na(chr) & !is.na(pos)], by=c("chr","pos"),
                  all.y=FALSE, all.x=FALSE)
      if (nrow(mm)) added_pos <- added_pos + nrow(mm)
    } else {
      next
    }
  }
  
  if (!nrow(mm)) next
  
  out <- mm[, .(
    gene    = if ("gene" %in% names(mm)) gene else NA_character_,
    genus   = sp_genus,                                    # enforce genus from file name
    snp     = if ("snp" %in% names(mm)) snp else NA_character_,
    chr     = if ("chr" %in% names(mm)) as.character(chr) else NA_character_,
    pos     = if ("pos" %in% names(mm)) as.integer(pos) else NA_integer_,
    species = sp,
    beta    = if (!is.na(beta_x)) as.numeric(get(beta_x)) else NA_real_,
    se_beta = if (!is.na(se_x))   as.numeric(get(se_x))   else NA_real_,
    p       = if (!is.na(p_x))    as.numeric(get(p_x))    else NA_real_
  )]
  out$chr=sent$chr[match(out$gene,sent$gene)]
  out$pos=sent$pos[match(out$gene,sent$gene)]
  
  out_list[[i]] <- out
}

res <- rbindlist(out_list, use.names=TRUE, fill=TRUE)

# De-dup within (gene,genus,species) – keep first (usually single match)
if (nrow(res)) res <- unique(res, by=c("gene","genus","species"))

fwrite(res, file=opt$out, sep="\t")
logf("Wrote %s (rows=%d, matched via variant_id=%d, via chr:pos=%d, species_files_skipped_by_genus=%d)",
     opt$out, nrow(res), added_vid, added_pos, skipped_by_genus)
