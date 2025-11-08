#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Biobase)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
})

LOG <- function(...) cat("[genus_assoc]", ..., "\n", sep = "")

main <- function() {
  option_list <- list(
    make_option("--gds",           type="character"),
    make_option("--scan_lookup",   type="character"),
    make_option("--scan_ids",      type="character"),
    make_option("--label_ids",     type="character"),
    make_option("--trait",         type="character"),
    make_option("--resid_rda",     type="character"),
    make_option("--metadata",      type="character"),
    make_option("--grm_rdata",     type="character"),
    make_option("--srm_dam",       type="character"),
    make_option("--srm_cage",      type="character"),
    make_option("--genera",        type="character", default="Prevotella,Bacteroides"),
    make_option("--trait_prefix",  type="character", default="g__"),
    make_option("--outdir",        type="character", default="assoc_genus"),
    make_option("--maf_min",       type="double",  default=NA_real_),
    make_option("--mac_min",       type="integer", default=NA_integer_),
    make_option("--miss_max",      type="double",  default=NA_real_),
    make_option("--hwe_p_min",     type="double",  default=NA_real_),
    make_option("--hwe_maf_min",   type="double",  default=0.05),
    make_option("--loco_grm_rds",  type="character", default=NA_character_)
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  VAR_BLOCK <- as.integer(Sys.getenv("GENESIS_VARBLOCK", "500"))
  LOG("VAR_BLOCK=", VAR_BLOCK)
  
  ## ---- Residuals (GENUS) ----
  load(opt$resid_rda)
  if (opt$trait=="Taxa"){
    if (!exists("residuals_qned_counts_objs"))
      stop("Expected 'residuals_qned_counts_objs' in --resid_rda")
    genus_layer_idx <- length(residuals_qned_counts_objs) - 1L
    residuals <- t(residuals_qned_counts_objs[[genus_layer_idx]])
    
  } else if (opt$trait=="Guild") {
    residuals <- t(residuals_qned_counts_clusters_objs[[1]])
    
  } else if (opt$trait=="Beta") {
    residuals <- t(residuals_qned_counts_beta_objs[[1]])
    
  } else if (opt$trait=="Alpha") {
    residuals <- t(residuals_qned_counts_alpha_objs[[1]])
  }
  rownames(residuals) <- sapply(strsplit(rownames(residuals), "_"), `[`, 1)
  
  genera <- strsplit(opt$genera, ",", fixed=TRUE)[[1]]
  # robust match: g__Prevotella or g__Prevotella_*
  patt <- paste0("^", opt$trait_prefix, "(", paste(genera, collapse="|"), ")", "(_|$)")
  traits <- grep(patt, colnames(residuals), value=TRUE)
  traits <- traits[traits!="g__Bacteroides_F"]
  stopifnot(length(traits) > 0)
  resid_g <- residuals[, traits, drop=FALSE]
  LOG("Traits selected: ", length(traits))
  
  ## ---- Metadata / GRM / SRMs ----
  meta <- fread(opt$metadata)
  id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
  if (is.na(id_col)) stop("Metadata must have an ID column among sample/Sample/RFID/ID/iid")
  meta <- meta[!duplicated(meta[[id_col]])]
  
  load(opt$grm_rdata); if (!exists("grm")) stop("No 'grm' in --grm_rdata")
  srm_dam  <- readRDS(opt$srm_dam)
  srm_cage <- readRDS(opt$srm_cage)
  
  ## ---- GDS + lookup ----
  gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
  lookup <- fread(opt$scan_lookup)
  stopifnot(all(c("scan_id","sample_label") %in% names(lookup)))
  int2label <- setNames(lookup$sample_label, lookup$scan_id)
  
  ## ---- Cohort ----
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
  
  ## Align
  resid_g  <- resid_g[lab,,drop=FALSE]
  grm      <- grm[lab, lab];        rownames(grm)      <- colnames(grm)      <- sel_char
  srm_dam  <- srm_dam[lab, lab];    rownames(srm_dam)  <- colnames(srm_dam)  <- sel_char
  srm_cage <- srm_cage[lab, lab];   rownames(srm_cage) <- colnames(srm_cage) <- sel_char
  
  ## ---- SeqArray subset ----
  seqData <- SeqVarData(gds, sampleData = AnnotatedDataFrame(
    data.frame(sample.id = seqGetData(gds, "sample.id"),
               row.names   = as.character(seqGetData(gds, "sample.id")))
  ))
  cohort_mask  <- seqGetData(seqData, "sample.id") %in% sel_int
  seqResetFilter(seqData, verbose=FALSE)
  seqSetFilter(seqData, sample.sel = cohort_mask, action="set", verbose=FALSE)
  cur <- seqGetFilter(seqData)
  LOG("After cohort filter  samples=", sum(cur$sample.sel), " variants=", sum(cur$variant.sel))
  
  ## ---- LOCO (optional) ----
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
  
  ## ---- QC on full cohort once ----
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
  
  seqResetFilter(seqData, verbose=FALSE)
  seqSetFilter(seqData, sample.sel = cohort_mask, action="set", verbose=FALSE)
  all_chr <- seqGetData(seqData, "chromosome")
  keep_train <- apply_qc(seqData, opt$miss_max, opt$maf_min, opt$mac_min,
                         opt$hwe_p_min, opt$hwe_maf_min)
  chr_levels <- sort(unique(all_chr[keep_train]), na.last=NA)
  keep_idx_by_chr <- lapply(chr_levels, function(ch) which((all_chr == ch) & keep_train))
  names(keep_idx_by_chr) <- as.character(chr_levels)
  LOG("QC kept variants: ", sum(keep_train))
  
  ## ---- Association ----
  cur_chr_ids <- sel_char
  cur_lab     <- as.character(int2label[cur_chr_ids])
  if (anyNA(cur_lab)) stop("int2label mapping produced NA for some scan IDs; check lookup.")
  pd_full <- data.frame(sample.id = as.integer(cur_chr_ids),
                        row.names = cur_chr_ids, stringsAsFactors = FALSE)
  
  for (tr in colnames(resid_g)) {
    y <- as.numeric(resid_g[cur_lab, tr])
    if (length(y) != length(cur_chr_ids)) stop("y length mismatch after label→scan alignment.")
    pd <- pd_full; pd$y <- y
    
    outdir_trait <- file.path(opt$outdir, "by_trait")
    dir.create(outdir_trait, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(outdir_trait, paste0(tr, ".tsv.gz"))
    
    res_parts <- list(); part <- 1L
    for (chr in names(keep_idx_by_chr)) {
      idx <- keep_idx_by_chr[[chr]]; if (length(idx) == 0) next
      seqResetFilter(seqData, verbose=FALSE)
      seqSetFilter(seqData, sample.sel = cohort_mask, variant.sel = idx, action="set", verbose=FALSE)
      
      covlist <- list(
        dam  = srm_dam [cur_chr_ids, cur_chr_ids, drop=FALSE],
        cage = srm_cage[cur_chr_ids, cur_chr_ids, drop=FALSE],
        grm  = if (use_loco) read_grm_prefix(loco_map[[as.character(chr)]], cur_chr_ids)
        else          grm[cur_chr_ids, cur_chr_ids, drop=FALSE]
      )
      
      nm <- try(fitNullModel(pd, outcome="y", covars=NULL,
                             cov.mat=covlist, family="gaussian", verbose=FALSE),
                silent=TRUE)
      if (inherits(nm, "try-error")) {
        LOG("[", tr, "] fitNullModel failed on chr ", chr, "; skipping chr")
        next
      }
      
      it <- SeqVarBlockIterator(seqData, variantBlock = VAR_BLOCK, verbose = FALSE)
      a  <- try(assocTestSingle(it, nm, test="Score", verbose=FALSE), silent=TRUE)
      if (!inherits(a, "try-error") && !is.null(a) && nrow(a) > 0) {
        keep_cols <- intersect(c("variant.id","chr","pos","Score.pval","Score","Score.SE","Est","Est.SE"),
                               names(a))
        res_parts[[part]] <- as.data.frame(a)[, keep_cols, drop=FALSE]
        part <- part + 1L
      }
      rm(nm, it, a); gc()
    }
    
    if (length(res_parts)) {
      res_all <- rbindlist(res_parts, use.names=TRUE, fill=TRUE)
      res_all[, logP := -log10(Score.pval)]
      res_all[, bonf_p := p.adjust(Score.pval, method = "bonferroni")]
      res_all[, qvalue := p.adjust(Score.pval, method = "fdr")]
      res_all[, Direction := fifelse(Est > 0, "Positive", fifelse(Score < 0, "Negative", "Zero"))]
      res_all[, Significance := fifelse(qvalue < 0.10 & bonf_p < 0.05, "Significant (Bonferroni)",
                                    fifelse(qvalue < 0.10, "Significant (FDR)", "Non-significant"))]
      res_all = res_all[order(res_all$logP, decreasing = T),]
      fwrite(res_all, out_file, sep="\t", compress="gzip")
    } else {
      con <- gzfile(out_file, "wt"); writeLines("", con); close(con)
    }
    LOG("Wrote: ", out_file)
  }
  
  fwrite(data.table(ok=TRUE, n_traits=length(traits)),
         file.path(opt$outdir, "genus_done.tsv"), sep="\t")
  LOG("DONE")
}

tryCatch(main(), error = function(e) { message("FATAL: ", conditionMessage(e)); quit(save="no", status=1) })
