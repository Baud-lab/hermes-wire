#!/usr/bin/env Rscript
# Gene-level association (Burden/SKAT-O/SMMAT/...) via GENESIS with LOCO + SRMs
# Writes: outdir/by_gene_by_trait/<trait>.geneset.tsv.gz  (GENE_ID, chr, start, end, n_var, p)

.required <- c(
  "optparse","data.table","Biobase","SeqArray","SeqVarTools",
  "GENESIS","GenomicRanges","survey","CompQuadForm"
)
.missing <- .required[!vapply(.required, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing)) {
  stop(sprintf("Missing R packages: %s\nInstall them in the container and re-run.",
               paste(.missing, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Biobase)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(GenomicRanges)
})

LOG <- function(...) cat("[set_assoc]", ..., "\n", sep = "")

opt_list <- list(
  make_option("--gds",            type="character"),
  make_option("--scan_lookup",    type="character"),
  make_option("--scan_ids",       type="character"),
  make_option("--label_ids",      type="character"),
  make_option("--resid_rda",      type="character"),
  make_option("--metadata",       type="character"),
  make_option("--grm_rdata",      type="character"),
  make_option("--srm_dam",        type="character"),
  make_option("--srm_cage",       type="character"),
  make_option("--genes_grl_rds",  type="character"),
  make_option("--genera",         type="character", default="Prevotella,Bacteroides"),
  make_option("--traits",         type="character"),
  make_option("--trait_prefix",   type="character", default="g__"),
  make_option("--outdir",         type="character", default="assoc_set"),
  make_option("--maf_min",        type="double",  default=NA_real_),
  make_option("--mac_min",        type="integer", default=NA_integer_),
  make_option("--miss_max",       type="double",  default=NA_real_),
  make_option("--hwe_p_min",      type="double",  default=NA_real_),
  make_option("--hwe_maf_min",    type="double",  default=0.05),
  make_option("--loco_grm_rds",   type="character", default=NA_character_),
  make_option("--method",         type="character", default="SMMAT") # SMMAT, SKAT, fastSKAT, Burden, SKATO, ...
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_bytrait <- file.path(opt$outdir, "by_gene_by_trait")
dir.create(out_bytrait, showWarnings = FALSE, recursive = TRUE)

## ---------- Residuals: choose genus or species matrix ----------
load(opt$resid_rda)
if (!exists("residuals_qned_counts_objs"))
  stop("Expected object 'residuals_qned_counts_objs' in --resid_rda")

resid_mat <- if (opt$trait_prefix == "g__") {
  t(residuals_qned_counts_objs[[length(residuals_qned_counts_objs)-1]])
} else {
  t(residuals_qned_counts_objs[[length(residuals_qned_counts_objs)]])
}
rownames(resid_mat) <- sapply(strsplit(rownames(resid_mat), "_"), `[`, 1)

genera <- trimws(strsplit(opt$genera, ",", fixed = TRUE)[[1]])
if (opt$trait_prefix == "g__") {
  wanted <- paste0("g__", genera)
  traits_all <- colnames(resid_mat)[colnames(resid_mat) %in% wanted]
} else {
  wanted_prefix <- paste0("^(", paste0("s__", genera, "_", collapse="|"), ")")
  traits_all <- grep(wanted_prefix, colnames(resid_mat), value = TRUE)
  traits_all <- traits_all[!grepl("^s__Bacteroides_F_", traits_all)]
}

traits <- traits_all
if (!is.null(opt$traits) && nzchar(opt$traits)) {
  want <- trimws(strsplit(opt$traits, ",", fixed=TRUE)[[1]])
  traits <- intersect(traits_all, want)
}
if (!length(traits)) {
  stop(sprintf("No traits matched --traits='%s' after genus filtering. Example available: %s",
               opt$traits %||% "", paste(head(traits_all,5), collapse=", ")))
}

LOG("Traits selected: ", length(traits))

## ---------- Cohort alignment & random effects ----------
meta <- fread(opt$metadata)
id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
if (is.na(id_col)) stop("Metadata must have an ID column among sample/Sample/RFID/ID/iid")
meta <- meta[!duplicated(meta[[id_col]])]

load(opt$grm_rdata); if (!exists("grm")) stop("No 'grm' in --grm_rdata")
srm_dam  <- readRDS(opt$srm_dam)
srm_cage <- readRDS(opt$srm_cage)

gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
lookup <- fread(opt$scan_lookup)
stopifnot(all(c("scan_id","sample_label") %in% names(lookup)))
label2int <- setNames(lookup$scan_id,  lookup$sample_label)
int2label <- setNames(lookup$sample_label, lookup$scan_id)

sel_int <- scan(opt$scan_ids,  what=integer(),  quiet=TRUE)
lab     <- scan(opt$label_ids, what=character(), quiet=TRUE)
if (length(sel_int) != length(lab)) stop("scan_ids and label_ids length mismatch")

lab_map <- as.character(int2label[as.character(sel_int)])
if (anyNA(lab_map)) stop("Some scan_ids missing in lookup")
if (!all(lab == lab_map)) lab <- lab_map

sel_char <- as.character(sel_int)
resid_mat <- resid_mat[lab, , drop=FALSE]
grm       <- grm[lab, lab];        rownames(grm)      <- colnames(grm)      <- sel_char
srm_dam   <- srm_dam[lab, lab];    rownames(srm_dam)  <- colnames(srm_dam)  <- sel_char
srm_cage  <- srm_cage[lab, lab];   rownames(srm_cage) <- colnames(srm_cage) <- sel_char

annot0  <- AnnotatedDataFrame(data.frame(
  sample.id = seqGetData(gds, "sample.id"),
  row.names = as.character(seqGetData(gds, "sample.id"))
))
svd <- SeqVarData(gds, sampleData = annot0)
all_ids_full <- seqGetData(svd, "sample.id")
cohort_mask  <- all_ids_full %in% sel_int
seqResetFilter(svd, verbose=FALSE)
seqSetFilter(svd, sample.sel = cohort_mask, action="set", verbose=FALSE)
cur <- seqGetFilter(svd)
LOG("After cohort filter  samples=", sum(cur$sample.sel), " variants=", sum(cur$variant.sel))

## ---------- LOCO map ----------
use_loco <- (!is.na(opt$loco_grm_rds) && nzchar(opt$loco_grm_rds))
loco_map <- if (use_loco) readRDS(opt$loco_grm_rds) else NULL

read_grm_prefix <- function(prefix, ids_char) {
  idp <- paste0(prefix, ".grm.id"); bin <- paste0(prefix, ".grm.bin")
  if (!file.exists(idp) || !file.exists(bin)) stop("Missing GRM files for prefix: ", prefix)
  ids <- read.table(idp, stringsAsFactors=FALSE)
  iid <- as.character(ids[[2]]); n <- length(iid)
  con <- file(bin, "rb"); on.exit(close(con), add=TRUE)
  lt  <- readBin(con, what=numeric(), n=n*(n+1)/2, size=4)
  M <- matrix(0, n, n); k <- 1L
  for (i in 1:n) for (j in 1:i) { M[i,j] <- lt[k]; M[j,i] <- lt[k]; k <- k+1L }
  rownames(M) <- iid; colnames(M) <- iid
  if (!all(ids_char %in% iid)) stop("LOCO GRM missing IDs")
  M[ids_char, ids_char, drop=FALSE]
}

## ---------- QC helper ----------
apply_qc <- function(sd, miss_max, maf_min, mac_min, hwe_p_min, hwe_maf_min) {
  nvar <- length(seqGetData(sd, "variant.id"))
  if (nvar == 0) return(rep(FALSE, 0))
  miss  <- missingGenotypeRate(sd, margin="by.variant")
  refAF <- alleleFrequency(sd, n=0); maf <- pmin(refAF, 1 - refAF)
  nsamp <- length(seqGetData(sd, "sample.id"))
  n_called <- round((1 - miss) * nsamp)
  mac <- round(maf * 2L * n_called)
  miss[is.na(miss)] <- 1; maf[is.na(maf)] <- 0; mac[is.na(mac)] <- 0L
  
  hwe_keep <- rep(TRUE, length(maf))
  if (!is.na(hwe_p_min)) {
    hw <- try(hwe(sd, permute=FALSE), silent=TRUE)
    if (!inherits(hw, "try-error") && is.data.frame(hw) && "p" %in% names(hw)) {
      hwe_p <- hw$p; hwe_p[is.na(hwe_p)] <- 0
      maf_cut <- ifelse(is.na(hwe_maf_min), 0.05, hwe_maf_min)
      test_ok <- (maf >= maf_cut)
      hwe_keep[test_ok] <- (hwe_p[test_ok] >= hwe_p_min)
    }
  }
  
  keep <- rep(TRUE, length(miss))
  if (!is.na(miss_max)) keep <- keep & (miss <= miss_max)
  if (!is.na(maf_min))  keep <- keep & (maf  >= maf_min)
  if (!is.na(mac_min))  keep <- keep & (mac  >= mac_min)
  keep & hwe_keep
}

seqResetFilter(svd, verbose=FALSE)
seqSetFilter(svd, sample.sel = cohort_mask, action="set", verbose=FALSE)
all_chr <- seqGetData(svd, "chromosome")
keep_train <- apply_qc(svd, opt$miss_max, opt$maf_min, opt$mac_min, opt$hwe_p_min, opt$hwe_maf_min)
chr_levels <- sort(unique(all_chr[keep_train]), na.last=NA)
keep_idx_by_chr <- lapply(chr_levels, function(ch) which((all_chr == ch) & keep_train))
names(keep_idx_by_chr) <- as.character(chr_levels)
LOG("QC kept variants:", sum(keep_train))

# Load gene GRangesList and keep single-chromosome genes
grl <- readRDS(opt$genes_grl_rds)
if (!inherits(grl, "GRangesList")) stop("genes_grl_rds is not a GRangesList")

gene_chr <- vapply(grl, function(g) {
  ch <- as.character(runValue(seqnames(reduce(g))))
  if (length(ch) == 1) ch else NA_character_
}, character(1))
keep_single <- !is.na(gene_chr)
grl <- grl[keep_single]
gene_chr <- gene_chr[keep_single]

# Sample-level design
cur_ids <- sel_char
pd_full <- data.frame(sample.id = as.integer(cur_ids), row.names = cur_ids)

## ---------- Main loop: per trait ----------
for (tr in traits) {
  y <- as.numeric(resid_mat[as.character(int2label[cur_ids]), tr])
  if (length(y) != length(cur_ids)) stop("y length mismatch after label→scan alignment.")
  
  out_file <- file.path(out_bytrait, paste0(tr, ".geneset.tsv.gz"))
  res_rows <- list(); part <- 1L
  
  for (chr in names(keep_idx_by_chr)) {
    idx <- keep_idx_by_chr[[chr]]
    if (length(idx) == 0) next
    
    seqResetFilter(svd, verbose=FALSE)
    seqSetFilter(svd, sample.sel = cohort_mask, variant.sel = idx, action="set", verbose=FALSE)
    
    cur <- seqGetFilter(svd)
    nvar_chr <- sum(cur$variant.sel)
    LOG("[", tr, "] chr=", chr, " nvar_postQC=", nvar_chr)
    if (nvar_chr == 0) next
    
    g_this <- grl[gene_chr == chr]
    if (length(g_this) == 0) next
    pos <- seqGetData(svd, "position")
    has_hits <- vapply(g_this, function(g) any(pos >= start(g)[1] & pos <= end(g)[length(g)]),
                       logical(1))
    nonempty <- which(has_hits)
    if (!length(nonempty)) next
    
    pd <- pd_full; pd$y <- y
    covlist <- list(
      dam  = srm_dam [cur_ids, cur_ids, drop=FALSE],
      cage = srm_cage[cur_ids, cur_ids, drop=FALSE],
      grm  = if (use_loco) read_grm_prefix(loco_map[[as.character(chr)]], cur_ids)
      else          grm[cur_ids, cur_ids, drop=FALSE]
    )
    nm <- try(fitNullModel(pd, outcome="y", covars=NULL,
                           cov.mat=covlist, family="gaussian", verbose=FALSE),
              silent=TRUE)
    if (inherits(nm, "try-error")) {
      LOG("[", tr, "] fitNullModel failed on chr ", chr, "; skipping chr")
      next
    }
    
    agg <- SeqVarListIterator(svd, variantRanges = g_this[nonempty], verbose = FALSE)
    ag  <- try(assocTestAggregate(agg, null.model = nm, test = opt$method, verbose = FALSE),
               silent = TRUE)
    if (inherits(ag, "try-error") || is.null(ag)) { rm(nm, agg); gc(); next }
    res <- try(ag$results, silent = TRUE)
    if (inherits(res, "try-error") || is.null(res) || !nrow(res)) { rm(nm, agg, ag); gc(); next }
    res_df <- as.data.frame(res)
    
    p_col <- switch(opt$method,
                    "SMMAT"      = "pval_SMMAT",
                    "SKAT"       = if ("pval" %in% names(res_df)) "pval" else "p.value",
                    "fastSKAT"   = if ("pval" %in% names(res_df)) "pval" else "p.value",
                    "Burden"     = if ("Score.pval" %in% names(res_df)) "Score.pval" else "p.value",
                    "SKATO"      = if ("pval_SKATO" %in% names(res_df)) "pval_SKATO" else "min.pval",
                    "p.value")
    if (!p_col %in% names(res_df)) {
      p_candidates <- intersect(c("p","p.value","pval","p_value","pval_SMMAT","pval_SKATO","Score.pval"),
                                names(res_df))
      if (!length(p_candidates)) { rm(nm, agg, ag); gc(); next }
      p_col <- p_candidates[1]
    }
    
    n_col <- if ("n.site" %in% names(res_df)) "n.site" else NULL
    
    gr_one <- unlist(range(g_this[nonempty]))
    starts <- start(gr_one); ends <- end(gr_one)
    
    df <- data.frame(
      GENE_ID = names(g_this)[nonempty],
      chr     = as.character(chr),
      start   = starts,
      end     = ends,
      n_var   = if (length(n_col)) res_df[[n_col]] else NA_integer_,
      p       = as.numeric(res_df[[p_col]]),
      stringsAsFactors = FALSE
    )
    res_rows[[part]] <- df; part <- part + 1L
    
    rm(nm, agg, ag, res_df); gc()
  } # chr
  
  if (length(res_rows)) {
    OUT <- rbindlist(res_rows, use.names = TRUE, fill = TRUE)
    setorder(OUT, chr, start, p)
    fwrite(OUT, file = out_file, sep = "\t", compress = "gzip")
  } else {
    con <- gzfile(out_file, "wt"); writeLines("", con); close(con)
  }
  LOG("Wrote: ", out_file)
}

fwrite(data.table(ok=TRUE, n_traits=length(traits)),
       file.path(opt$outdir, "set_done.tsv"), sep="\t")
LOG("DONE")

`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || !nzchar(a)) b else a
