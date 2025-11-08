#!/usr/bin/env Rscript

## ======================================================================
## SPECIES_BETAS_FOR_SENTINELS_GENUS_FOLD — robust + verbose
## - Basic syntax only
## - Hard guards before using rownames()/colnames()
## - Mirrors run_genesis_genus.R alignment:
##     * Use LABEL domain to extract phenotypes (y)
##     * Subset GRM/SRMs by LABELS, then rename to SCAN IDs
## - Prints before/after counts and brief LMM summaries
## ======================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Biobase)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
})

LOG <- function(...) cat("[species_oof]", paste(..., collapse = " "), "\n", sep = "")

## ---------- helpers (simple, defensive) ----------
stopf <- function(...) { msg <- paste(..., collapse=" "); cat("FATAL:", msg, "\n"); quit(save="no", status=1) }

# Check that an object has rownames/colnames before using them.
check_has_rownames <- function(x, what, required=TRUE, preview_n=5L) {
  rn <- rownames(x)
  has <- !is.null(rn)
  LOG(sprintf("%s rownames present : %s", what, if (has) "YES" else "NO"))
  if (has) {
    nprev <- min(preview_n, length(rn))
    LOG(sprintf("%s example rownames : %s", what, paste(utils::head(rn, nprev), collapse=", ")))
  }
  if (required && !has) stopf(what, " is missing rownames; cannot match IDs safely.")
  invisible(has)
}
check_has_colnames <- function(x, what, required=TRUE, preview_n=5L) {
  cn <- colnames(x)
  has <- !is.null(cn)
  LOG(sprintf("%s colnames present : %s", what, if (has) "YES" else "NO"))
  if (has) {
    nprev <- min(preview_n, length(cn))
    LOG(sprintf("%s example colnames : %s", what, paste(utils::head(cn, nprev), collapse=", ")))
  }
  if (required && !has) stopf(what, " is missing colnames; cannot match traits safely.")
  invisible(has)
}

trim_after_underscore <- function(x) {
  x <- as.character(x)
  sub("^([^_]+).*", "\\1", x, perl=TRUE)
}

norm_chr_string <- function(x) {
  # Normalize chromosome values from sentinels: drop "chr", trim spaces, lowercase→character
  y <- as.character(x)
  y <- trimws(y)
  y <- tolower(y)
  y <- sub("^chr", "", y)
  y
}

## ---------- CLI ----------
main <- function() {
  option_list <- list(
    make_option("--gds",              type="character"),
    make_option("--scan_lookup",      type="character"),
    make_option("--scan_ids_test",    type="character"),
    make_option("--label_ids_test",   type="character"),
    make_option("--trait",            type="character"),
    make_option("--resid_rda",        type="character"),
    make_option("--metadata",         type="character"),
    make_option("--grm_rdata",        type="character"),
    make_option("--srm_dam",          type="character"),
    make_option("--srm_cage",         type="character"),
    make_option("--loco_grm_rds",     type="character", default=NA_character_),
    make_option("--genus",            type="character"),
    make_option("--species_prefix",   type="character", default="s__"),
    make_option("--sentinels_tsv",    type="character"),
    make_option("--out",              type="character", default="beta_cf_by_species.tsv"),
    make_option("--var_block",        type="integer",   default=500),
    make_option("--maf_min",          type="double",    default=NA_real_),
    make_option("--mac_min",          type="integer",   default=NA_integer_),
    make_option("--miss_max",         type="double",    default=NA_real_),
    make_option("--hwe_p_min",        type="double",    default=NA_real_),
    make_option("--hwe_maf_min",      type="double",    default=0.05)
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  
  ## ---------- sanity on files ----------
  if (!file.exists(opt$gds))            stopf("GDS not found:", opt$gds)
  if (!file.exists(opt$scan_lookup))    stopf("scan_lookup TSV not found:", opt$scan_lookup)
  if (!file.exists(opt$scan_ids_test))  stopf("scan_ids_test not found:", opt$scan_ids_test)
  if (!file.exists(opt$label_ids_test)) stopf("label_ids_test not found:", opt$label_ids_test)
  if (!file.exists(opt$resid_rda))      stopf("resid_rda not found:", opt$resid_rda)
  if (!file.exists(opt$metadata))       stopf("metadata not found:", opt$metadata)
  if (!file.exists(opt$grm_rdata))      stopf("grm_rdata not found:", opt$grm_rdata)
  if (!file.exists(opt$srm_dam))        stopf("srm_dam not found:", opt$srm_dam)
  if (!file.exists(opt$srm_cage))       stopf("srm_cage not found:", opt$srm_cage)
  if (!file.exists(opt$sentinels_tsv))  stopf("sentinels_tsv not found:", opt$sentinels_tsv)
  
  ## ---------- residuals: load + pick species level ----------
  load(opt$resid_rda)
  if (opt$trait=="Taxa"){
    if (!exists("residuals_qned_counts_objs"))
      stop("Expected 'residuals_qned_counts_objs' in --resid_rda")
    genus_layer_idx <- length(residuals_qned_counts_objs)
    species_mat <- t(residuals_qned_counts_objs[[genus_layer_idx]])

  } else if (opt$trait=="Guild") {
    species_mat <- t(residuals_qned_counts_clusters_objs[[1]])

  } else if (opt$trait=="Beta") {
    species_mat <- t(residuals_qned_counts_beta_objs[[1]])
    
  } else if (opt$trait=="Alpha") {
    species_mat <- t(residuals_qned_counts_alpha_objs[[1]])
  }
  rownames(species_mat) <- sapply(strsplit(rownames(species_mat), "_"), `[`, 1)

  # Guards before ANY matching:
  check_has_rownames(species_mat, "Residuals(species)", required=TRUE)
  check_has_colnames(species_mat, "Residuals(species)", required=TRUE)
  
  ## ---------- select species traits for the requested genus ----------
  all_traits <- colnames(species_mat)
  if (opt$trait=="Taxa"){
    patt <- paste0("^", opt$species_prefix, opt$genus, "(_|$)")
    traits <- grep(patt, all_traits, value=TRUE)
    traits=traits[!grepl("s__Bacteroides_F_",traits)]
  } else {
    patt <- paste0(opt$species_prefix, opt$genus)
    traits <- all_traits[all_traits==patt]
  }
  LOG("Species traits for genus=", opt$genus, " :", length(traits))
  if (length(traits) == 0) {
    LOG("Residuals(species) trait prefix:", opt$species_prefix,
        " | example traits:", paste(utils::head(all_traits, 5), collapse=", "))
    stopf("No species traits found for genus=", opt$genus)
  }
  species_mat <- species_mat[, traits, drop=FALSE]
  
  ## ---------- metadata (just for reporting ID column) ----------
  meta <- fread(opt$metadata)
  id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
  if (is.na(id_col)) stopf("Metadata must contain an ID column among sample/Sample/RFID/ID/iid")
  meta <- meta[!duplicated(meta[[id_col]])]
  LOG("Metadata rows       :", nrow(meta), "| ID column detected:", id_col)
  
  ## ---------- GRM / SRMs (guards on dimnames) ----------
  load(opt$grm_rdata); if (!exists("grm")) stopf("No object 'grm' in --grm_rdata")
  srm_dam  <- readRDS(opt$srm_dam)
  srm_cage <- readRDS(opt$srm_cage)
  LOG("GRM dims            :", nrow(grm), "\u00D7", ncol(grm))
  LOG("SRM(dam) dims       :", nrow(srm_dam), "\u00D7", ncol(srm_dam))
  LOG("SRM(cage) dims      :", nrow(srm_cage), "\u00D7", ncol(srm_cage))
  
  check_has_rownames(grm,     "GRM",      required=TRUE)
  check_has_colnames(grm,     "GRM",      required=TRUE)
  check_has_rownames(srm_dam, "SRM(dam)", required=TRUE)
  check_has_colnames(srm_dam, "SRM(dam)", required=TRUE)
  check_has_rownames(srm_cage,"SRM(cage)",required=TRUE)
  check_has_colnames(srm_cage,"SRM(cage)",required=TRUE)
  
  ## ---------- GDS + lookup ----------
  gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
  lookup <- fread(opt$scan_lookup)
  if (!all(c("scan_id","sample_label") %in% names(lookup)))
    stopf("scan_lookup must have columns: scan_id, sample_label")
  LOG("Lookup rows         :", nrow(lookup))
  int2label <- setNames(as.character(lookup$sample_label), as.character(lookup$scan_id))
  
  ## ---------- test cohort ----------
  sel_int <- scan(opt$scan_ids_test,  what=integer(),  quiet=TRUE)
  lab_raw <- scan(opt$label_ids_test, what=character(), quiet=TRUE)
  LOG("TEST scan_ids (raw) :", length(sel_int))
  LOG("TEST label_ids (raw):", length(lab_raw))
  if (length(sel_int) != length(lab_raw))
    stopf("scan_ids_test and label_ids_test lengths differ (must be paired)")
  
  # labels derived from scan_ids via lookup (authoritative)
  lab_map <- as.character(int2label[as.character(sel_int)])
  if (any(is.na(lab_map))) stopf("Some scan_ids have no label in scan_lookup")
  n_disagree <- sum(lab_raw != lab_map, na.rm=TRUE)
  LOG("Label disagreement  :", n_disagree, "(between file labels and lookup mapping)")
  lab <- lab_map
  lab_norm <- trim_after_underscore(lab)    # label domain used by residuals
  
  ## ---------- residuals coverage BEFORE indexing ----------
  res_ids <- rownames(species_mat)  # guaranteed non-NULL by our guard
  cover_label <- sum(lab_norm %in% res_ids)
  cover_scan  <- sum(as.character(sel_int) %in% res_ids)
  LOG(sprintf("Residuals coverage  : label=%d/%d | scan_id=%d/%d",
              cover_label, length(lab_norm), cover_scan, length(sel_int)))
  domain_used <- if (cover_label == length(lab_norm)) "label"
  else if (cover_scan == length(sel_int)) "scan_id"
  else "mixed"
  LOG(sprintf("Residuals rownames  : domain=%s | coverage(label)=%d/%d | coverage(scan_id)=%d/%d",
              domain_used, cover_label, length(lab_norm), cover_scan, length(sel_int)))
  if (!(cover_label == length(lab_norm) || cover_scan == length(sel_int))) {
    stopf("Residuals rownames do not match label or scan_id entirely. Please standardize residual IDs.")
  }
  
  # Align residuals to test cohort using LABEL domain (mirrors genus script)
  y_mat <- species_mat[lab_norm, , drop=FALSE]
  LOG("TEST after resid    :", nrow(y_mat), " samples ( ", length(unique(rownames(y_mat))), " ids unique)")
  
  ## ---------- SeqArray: subset to cohort ----------
  seqData <- SeqVarData(gds, sampleData = AnnotatedDataFrame(
    data.frame(sample.id = seqGetData(gds, "sample.id"),
               row.names   = as.character(seqGetData(gds, "sample.id")))
  ))
  cohort_mask <- seqGetData(seqData, "sample.id") %in% sel_int
  seqResetFilter(seqData, verbose=FALSE)
  seqSetFilter(seqData, sample.sel = cohort_mask, action="set", verbose=FALSE)
  cur <- seqGetFilter(seqData)
  nsamp_filt <- sum(cur$sample.sel)
  nvar_all   <- sum(cur$variant.sel)
  LOG("GDS sample.id count :", length(seqGetData(seqData, "sample.id")))
  LOG(sprintf("GDS TEST matches    : %d/%d (TEST in GDS)", nsamp_filt, length(sel_int)))
  if (nsamp_filt != length(sel_int)) stopf("Not all TEST scan_ids found in GDS sample.id")
  
  ## ---------- QC (once) ----------
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
  LOG("After cohort filter  samples=", sum(cur$sample.sel), " variants=", nvar_all)
  LOG("QC kept variants    :", sum(keep_train))
  
  ## ---------- Align GRM/SRMs to labels, then rename to scan_ids ----------
  sel_char <- as.character(sel_int)  # scan IDs as character
  grm_cut     <- grm[lab_norm, lab_norm, drop=FALSE]
  srm_dam_cut <- srm_dam[lab_norm, lab_norm, drop=FALSE]
  srm_cag_cut <- srm_cage[lab_norm, lab_norm, drop=FALSE]
  # Rename to scan IDs so covariances align with sample.id in GENESIS
  rownames(grm_cut)     <- colnames(grm_cut)     <- sel_char
  rownames(srm_dam_cut) <- colnames(srm_dam_cut) <- sel_char
  rownames(srm_cag_cut) <- colnames(srm_cag_cut) <- sel_char
  
  ## ---------- Sentinels: read + validate ----------
  sent <- fread(opt$sentinels_tsv)
  # Accept typical headers: gene, genus, variant_id, SNP, P, CHR, BP (case-insensitive)
  cn <- tolower(names(sent))
  names(sent) <- cn
  if (!("chr" %in% names(sent) || "bp" %in% names(sent))) {
    LOG("Sentinels header:", paste(names(sent), collapse=","))
    stopf("Sentinels TSV must include 'CHR' and 'BP' (chromosome and position).")
  }
  # Normalize chr + keep required cols
  if ("chr" %in% names(sent))  sent$chr <- norm_chr_string(sent$chr)
  if ("bp"  %in% names(sent))  sent$bp  <- as.integer(sent$bp)
  sent <- sent[!is.na(chr) & !is.na(bp)]
  LOG("Sentinels header:", paste(names(sent), collapse=","))
  LOG("Sentinels rows      :", nrow(sent))
  
  # Pre-compute how many targets per chr and whether we can match them post-QC
  target_by_chr <- split(sent$bp, sent$chr)
  for (chrk in names(target_by_chr)) {
    want <- length(target_by_chr[[chrk]])
    idx_pool <- keep_idx_by_chr[[chrk]]
    matched <- 0L
    if (!is.null(idx_pool) && length(idx_pool) > 0) {
      seqResetFilter(seqData, verbose=FALSE)
      seqSetFilter(seqData, sample.sel = cohort_mask, variant.sel = idx_pool, action="set", verbose=FALSE)
      cur_pos <- seqGetData(seqData, "position")
      matched <- sum(cur_pos %in% target_by_chr[[chrk]])
    }
    LOG("[targets chr ", chrk, "] targets=", want, " | matched=", matched)
  }
  
  ## ---------- Per species loop (iterate over residual columns) ----------
  out_rows <- list()
  out_i <- 1L
  
  # Phenotype dataframe for GENESIS (sample.id = scan_id; rownames=scan_id)
  pd_full <- data.frame(sample.id = as.integer(sel_char),
                        row.names = sel_char, stringsAsFactors = FALSE)
  
  for (tr in traits) {
    species_name <- tr
    
    # Build per-chromosome index set from sentinels after QC (exact position match)
    sp_rows <- sent  # sentinels are genus-level; reuse same set for each species
    sp_by_chr <- split(sp_rows, sp_rows$chr)
    matched_idx_by_chr <- list()
    
    for (chrk in names(sp_by_chr)) {
      sr <- sp_by_chr[[chrk]]
      tar_pos <- as.integer(sr$bp)
      idx_pool <- keep_idx_by_chr[[chrk]]
      if (is.null(idx_pool) || length(idx_pool) == 0) {
        LOG("[ ", species_name, " ] chr ", chrk, " : targets=", length(tar_pos),
            " | matched=0 (no QC-kept variants)")
        next
      }
      seqResetFilter(seqData, verbose=FALSE)
      seqSetFilter(seqData, sample.sel = cohort_mask, variant.sel = idx_pool, action="set", verbose=FALSE)
      cur_pos <- seqGetData(seqData, "position")
      idx_match <- idx_pool[cur_pos %in% tar_pos]
      LOG("[ ", species_name, " ] chr ", chrk, " : targets=", length(tar_pos),
          " | matched=", length(idx_match), " (by_pos=", length(idx_match), ")")
      if (length(idx_match) > 0) matched_idx_by_chr[[chrk]] <- idx_match
    }
    
    if (length(matched_idx_by_chr) == 0) {
      next
    }
    
    ## ================== IMPORTANT FIX ==================
    ## Phenotype vector MUST be indexed by LABELS (lab_norm),
    ## not by scan IDs, because y_mat rownames are labels.
    ## This matches the genus script and avoids "subscript out of bounds".
    y <- as.numeric(y_mat[lab_norm, tr])
    if (length(y) != length(sel_char)) {
      stopf("y length mismatch after label→scan alignment (expected ", length(sel_char),
            ", got ", length(y), ")")
    }
    ## ===================================================
    
    # Fit null and test per chr block
    res_parts <- list(); part <- 1L
    for (chrk in names(matched_idx_by_chr)) {
      idx <- matched_idx_by_chr[[chrk]]
      if (length(idx) == 0) next
      
      seqResetFilter(seqData, verbose=FALSE)
      seqSetFilter(seqData, sample.sel = cohort_mask, variant.sel = idx, action="set", verbose=FALSE)
      
      covlist <- list(
        dam  = srm_dam_cut[sel_char, sel_char, drop=FALSE],
        cage = srm_cag_cut[sel_char, sel_char, drop=FALSE],
        grm  = grm_cut    [sel_char, sel_char, drop=FALSE]
      )
      # NOTE: If using LOCO later, ensure GRM for that chr is read and subset to sel_char.
      # GENESIS requires cov.mat dimnames to match sample.id. (See docs.)
      # (https://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_single.html)
      
      # Build a trait-specific pd quickly (keeps order identical to sel_char)
      pd <- pd_full
      pd$y <- y
      
      nm <- try(fitNullModel(pd, outcome="y", covars=NULL,
                             cov.mat=covlist, family="gaussian", verbose=FALSE),
                silent=TRUE)
      if (inherits(nm, "try-error")) {
        LOG("[ ", species_name, " ] fitNullModel failed on chr ", chrk, "; skipping chr")
        next
      }
      
      it <- SeqVarBlockIterator(seqData, variantBlock = opt$var_block, verbose = FALSE)
      a  <- try(assocTestSingle(it, nm, test="Score", verbose=FALSE), silent=TRUE)
      if (!inherits(a, "try-error") && !is.null(a) && nrow(a) > 0) {
        aa <- as.data.frame(a)
        keep_cols <- intersect(c("variant.id","chr","pos","Score.pval","Score","Score.SE","Est","Est.SE"),
                               names(aa))
        aa <- aa[, keep_cols, drop=FALSE]
        aa$trait <- species_name
        res_parts[[part]] <- aa
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
      # brief per-trait summary (top 3)
      ord <- order(res_all$Score.pval, decreasing=FALSE)
      top_n <- min(3L, nrow(res_all))
      top_msg <- paste(utils::head(sprintf("chr%s:%s p=%.3g Est=%.3g",
                                           res_all$chr[ord], res_all$pos[ord],
                                           res_all$Score.pval[ord], res_all$Est[ord]), top_n),
                       collapse=" | ")
      LOG("[ ", species_name, " ] LMM top hits   :", top_msg)
      
      out_rows[[out_i]] <- res_all
      out_i <- out_i + 1L
    } else {
      LOG("[ ", species_name, " ] No tested variants produced results after filtering.")
    }
  }
  
  ## ---------- write output ----------
  if (length(out_rows)) {
    OUT <- rbindlist(out_rows, use.names=TRUE, fill=TRUE)
    data.table::fwrite(OUT, opt$out, sep="\t")
    LOG("Wrote:", opt$out, " (rows=", nrow(OUT), ")")
  } else {
    fwrite(data.table(), opt$out, sep="\t")
    LOG("Wrote:", opt$out, " (rows=0)")
  }
  
  LOG("DONE")
}

tryCatch(main(), error=function(e) { cat("FATAL:", conditionMessage(e), "\n"); quit(save="no", status=1) })
