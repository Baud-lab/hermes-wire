#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(Biobase)
  library(SeqArray); library(SeqVarTools); library(GENESIS)
})

# Per-block variant count (override with GENESIS_VARBLOCK env var)
VAR_BLOCK <- as.integer(Sys.getenv("GENESIS_VARBLOCK", "500"))

# ---------------------- CLI ----------------------
opt <- parse_args(OptionParser(option_list=list(
  make_option("--gds",           type="character"),
  make_option("--scan_lookup",   type="character"),
  make_option("--scan_ids",      type="character"),
  make_option("--label_ids",     type="character"),
  make_option("--resid_rda",     type="character"),
  make_option("--metadata",      type="character"),
  make_option("--grm_rdata",     type="character"),
  make_option("--srm_dam",       type="character"),
  make_option("--srm_cage",      type="character"),
  make_option("--genera",        type="character", default="Prevotella,Bacteroides"),
  make_option("--trait_prefix",  type="character", default="s__"),
  make_option("--outdir",        type="character", default="assoc"),
  # QC thresholds (applied on TRAIN only)
  make_option("--maf_min",       type="double",  default=NA),
  make_option("--mac_min",       type="integer", default=NA),
  make_option("--miss_max",      type="double",  default=NA),
  make_option("--hwe_p_min",     type="double",  default=NA),
  make_option("--hwe_maf_min",   type="double",  default=NA),
  # LOCO
  make_option("--loco_grm_rds",  type="character", default=NA),
  # K-folds
  make_option("--k_folds",       type="integer", default=5),
  make_option("--seed",          type="integer", default=17)
)))

log <- function(...) cat("[species_assoc]", ..., "\n")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ------------------ Load residuals ------------------
load(opt$resid_rda)
if (!exists("residuals_qned_counts_objs"))
  stop("Expected 'residuals_qned_counts_objs' in --resid_rda")
residuals <- t(residuals_qned_counts_objs[[length(residuals_qned_counts_objs)]])
rownames(residuals) <- sapply(strsplit(rownames(residuals), "_"), `[`, 1)

genera <- strsplit(opt$genera, ",", fixed=TRUE)[[1]]
wanted <- unlist(lapply(genera, function(g) paste0(opt$trait_prefix, g, "_")))
traits <- grep(paste0("^(", paste0(wanted, collapse="|"), ")"), colnames(residuals), value=TRUE)
traits <- traits[!grepl("s__Bacteroides_F_", traits)]
stopifnot(length(traits) > 0)
resid_s <- residuals[, traits, drop=FALSE]

# --------------- Metadata / GRM / SRMs ---------------
meta <- fread(opt$metadata)
id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
if (is.na(id_col)) stop("Metadata must have an ID column among sample/Sample/RFID/ID/iid")
meta <- meta[!duplicated(meta[[id_col]])]

load(opt$grm_rdata); if (!exists("grm")) stop("No 'grm' in --grm_rdata")
srm_dam  <- readRDS(opt$srm_dam)
srm_cage <- readRDS(opt$srm_cage)

# ------------------- GDS + mapping -------------------
gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)

lookup <- fread(opt$scan_lookup)
stopifnot(all(c("scan_id","sample_label") %in% names(lookup)))
label2int <- setNames(lookup$scan_id,  lookup$sample_label)
int2label <- setNames(lookup$sample_label, lookup$scan_id)

# ---------------- Fixed (848) cohort ----------------
sel_int <- scan(opt$scan_ids,  what=integer(),  quiet=TRUE)
lab     <- scan(opt$label_ids, what=character(), quiet=TRUE)
if (length(sel_int) != length(lab))
  stop("scan_ids and label_ids lengths differ (must be paired)")

lab_map <- as.character(int2label[as.character(sel_int)])
if (anyNA(lab_map)) stop("Some scan_ids have no label in scan_lookup")
if (!all(lab == lab_map)) {
  warning("label_ids do not exactly match lookup; using mapping derived from scan_ids")
  lab <- lab_map
}

sel_char <- as.character(sel_int)
resid_s  <- resid_s[lab,,drop=FALSE]
grm      <- grm[lab, lab];        rownames(grm)      <- colnames(grm)      <- sel_char
srm_dam  <- srm_dam[lab, lab];    rownames(srm_dam)  <- colnames(srm_dam)  <- sel_char
srm_cage <- srm_cage[lab, lab];   rownames(srm_cage) <- colnames(srm_cage) <- sel_char

annot0 <- AnnotatedDataFrame(data.frame(
  sample.id = seqGetData(gds, "sample.id"),
  row.names = as.character(seqGetData(gds, "sample.id"))
))
seqData <- SeqVarData(gds, sampleData=annot0)

all_ids_full <- seqGetData(seqData, "sample.id")
cohort_mask  <- all_ids_full %in% sel_int

seqResetFilter(seqData, verbose=FALSE)
seqSetFilter(seqData, sample.sel = cohort_mask, action="set", verbose=FALSE)
cur <- seqGetFilter(seqData)
cat("[species_assoc] after cohort filter  samples=",
    sum(cur$sample.sel), " variants=", sum(cur$variant.sel), "\n", sep="")
seqResetFilter(seqData, verbose=FALSE)

# ------------------- Helpers -------------------
apply_qc_f <- function(sd, keep_var_base, miss_max, maf_min, mac_min, hwe_p_min, hwe_maf_min) {
  nvar <- length(seqGetData(sd, "variant.id"))
  if (nvar == 0) return(rep(FALSE, 0))
  miss  <- missingGenotypeRate(sd, margin="by.variant")
  refAF <- alleleFrequency(sd, n=0); maf <- pmin(refAF, 1-refAF)
  nsamp <- length(seqGetData(sd,"sample.id"))
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
  keep & hwe_keep & keep_var_base
}

use_loco <- !is.na(opt$loco_grm_rds) && nzchar(opt$loco_grm_rds)
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

# --------------- Grouped K-folds ---------------
set.seed(opt$seed)
group_col <- if ("dam" %in% names(meta)) "dam" else if ("cage" %in% names(meta)) "cage" else NA
groups <- if (!is.na(group_col)) meta[match(lab, meta[[id_col]]), ..group_col][[1]] else lab
groups <- as.character(groups)
uniq_g <- unique(groups); k <- min(opt$k_folds, length(uniq_g))
perm   <- sample(uniq_g, length(uniq_g))
fold_g <- split(perm, rep(1:k, length.out=length(perm)))
folds  <- lapply(fold_g, function(gs) lab[groups %in% gs])  # test labels per fold

log(sprintf("Cohort fixed from MAKE_GENESIS_COHORT_IDS: n=%d", length(lab)))

sum_rows <- list()

# ================= MAIN LOOP =================
for (fi in seq_along(folds)) {
  test_lab  <- folds[[fi]]
  train_lab <- setdiff(lab, test_lab)
  
  fold_dir <- file.path(opt$outdir, sprintf("fold%02d", fi))
  dir.create(fold_dir, showWarnings=FALSE, recursive=TRUE)
  
  # -------- TRAIN QC ONCE PER FOLD --------
  train_mask_full <- cohort_mask & (all_ids_full %in% as.integer(label2int[train_lab]))
  
  seqResetFilter(seqData, verbose=FALSE)
  seqSetFilter(seqData, sample.sel = train_mask_full, action="set", verbose=FALSE)
  
  all_vid <- seqGetData(seqData,"variant.id")
  all_chr <- seqGetData(seqData,"chromosome")
  keep_train <- apply_qc_f(seqData, rep(TRUE, length(all_vid)),
                           opt$miss_max, opt$maf_min, opt$mac_min,
                           opt$hwe_p_min, opt$hwe_maf_min)
  chr_levels <- sort(unique(all_chr[keep_train]), na.last=NA)
  keep_idx_by_chr <- lapply(chr_levels, function(ch) which((all_chr==ch) & keep_train))
  names(keep_idx_by_chr) <- as.character(chr_levels)
  
  # closure to run train/test subsets using TRAIN-QC variant indices
  run_subset <- function(cur_lab, subset_name) {
    subset_mask_full <- cohort_mask & (all_ids_full %in% as.integer(label2int[cur_lab]))
    
    seqResetFilter(seqData, verbose=FALSE)
    seqSetFilter(seqData, sample.sel = subset_mask_full, action="set", verbose=FALSE)
    samp_ids_now <- seqGetData(seqData, "sample.id")
    pd_all <- data.frame(sample.id = samp_ids_now,
                         row.names = as.character(samp_ids_now),
                         stringsAsFactors = FALSE)
    
    nonempty_traits <- 0L
    for (tr in colnames(resid_s)) {
      y_all <- rep(NA_real_, nrow(pd_all))
      cur_chr_ids <- as.character(as.integer(label2int[cur_lab]))
      y_all[match(cur_chr_ids, rownames(pd_all))] <- as.numeric(resid_s[cur_lab, tr])
      pd_all$y <- y_all
      pd_sel <- pd_all[cur_chr_ids, , drop=FALSE]
      
      res_parts <- list(); part <- 1L
      for (chr in names(keep_idx_by_chr)) {
        idx <- keep_idx_by_chr[[chr]]
        if (length(idx) == 0) next
        
        seqResetFilter(seqData, verbose=FALSE)
        seqSetFilter(seqData,
                     sample.sel  = subset_mask_full,
                     variant.sel = idx,
                     action="set", verbose=FALSE)
        
        covlist <- list(
          dam  = srm_dam[cur_chr_ids, cur_chr_ids, drop=FALSE],
          cage = srm_cage[cur_chr_ids, cur_chr_ids, drop=FALSE],
          grm  = if (use_loco) read_grm_prefix(loco_map[[as.character(chr)]], cur_chr_ids) else grm[cur_chr_ids, cur_chr_ids, drop=FALSE]
        )
        
        nm <- try(fitNullModel(pd_sel, outcome="y", covars=NULL,
                               cov.mat=covlist, family="gaussian", verbose=FALSE), silent=TRUE)
        if (inherits(nm, "try-error")) {
          log(sprintf("[fold%02d][%s][%s][chr=%s] fitNullModel failed; skipping", fi, subset_name, tr, chr))
          next
        }
        
        it <- SeqVarBlockIterator(seqData, variantBlock=VAR_BLOCK, verbose=FALSE)
        a  <- try(assocTestSingle(it, nm, test="Score", verbose=FALSE), silent=TRUE)
        if (!inherits(a, "try-error") && !is.null(a) && nrow(a) > 0) {
          # Keep only essential columns to reduce memory/IO
          keep_cols <- intersect(c("variant.id","chr","pos","Score.pval","Score","Score.SE","Est","Est.SE"), names(a))
          a <- a[, keep_cols, drop=FALSE]
          res_parts[[part]] <- as.data.frame(a); part <- part + 1L
        }
        rm(nm, it, a); gc()
      } # chr
      
      out_dir <- file.path(fold_dir, subset_name); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
      fn <- file.path(out_dir, paste0(tr, "_assoc.tsv.gz"))
      if (length(res_parts)) {
        res_all <- rbindlist(res_parts, use.names=TRUE, fill=TRUE)
        fwrite(res_all, fn, sep="\t", compress="gzip")
        nonempty_traits <- nonempty_traits + 1L
      } else {
        con <- gzfile(fn, "wt"); writeLines("", con); close(con)
      }
    } # trait
    nonempty_traits
  } # run_subset
  
  n_train <- run_subset(train_lab, "train")
  n_test  <- run_subset(test_lab,  "test")
  sum_rows[[length(sum_rows)+1L]] <- data.table(fold=sprintf("fold%02d",fi), train_nonempty=n_train, test_nonempty=n_test)
  log(sprintf("[fold%02d] train_nonempty=%d test_nonempty=%d", fi, n_train, n_test))
} # folds

fwrite(rbindlist(sum_rows), file.path(opt$outdir,"fold_summary.tsv"), sep="\t")
log("DONE")
