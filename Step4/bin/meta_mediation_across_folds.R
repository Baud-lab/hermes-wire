#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--inputs-list", type="character", help="One per-fold TSV per line (now may contain both mediation and mvmr rows)"),
  make_option("--sig_metric",  type="character", default="q", help="p or q for meta-level calls (mediation branch)"),
  make_option("--fdr_level",   type="double",    default=0.10),
  make_option("--p_level",     type="double",    default=0.05),
  # IV meta options (mediation branch)
  make_option("--f_threshold",       type="double", default=10),
  make_option("--steiger_min_prop",  type="double", default=0.60),
  make_option("--out",               type="character", default="mediation_meta.tsv")
)))

stopifnot(file.exists(opt$`inputs-list`))
files <- unique(trimws(readLines(opt$`inputs-list`)))
files <- files[nzchar(files) & file.exists(files)]
if (!length(files)) stop("No readable inputs")

DT <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  setDT(x); setnames(x, tolower(names(x)))
  x
}), fill=TRUE)

if (!"analysis" %in% names(DT)) DT[, analysis := "med"]  # back-compat

## ========= helpers =========
meta_one <- function(est, se) {
  w <- 1/(se^2); w[!is.finite(w)] <- NA_real_
  if (all(is.na(w))) {
    estm <- if (all(is.na(est))) NA_real_ else mean(est, na.rm=TRUE)
    return(list(est=estm, se=NA_real_, z=NA_real_, p=NA_real_, Q=NA_real_, Qdf=NA_real_, Qp=NA_real_, I2=NA_real_))
  }
  w[is.na(w)] <- 0
  sw <- sum(w)
  if (sw <= 0) {
    estm <- if (all(is.na(est))) NA_real_ else mean(est, na.rm=TRUE)
    return(list(est=estm, se=NA_real_, z=NA_real_, p=NA_real_, Q=NA_real_, Qdf=NA_real_, Qp=NA_real_, I2=NA_real_))
  }
  estm <- sum(w*est, na.rm=TRUE)/sw
  sem  <- sqrt(1/sw)
  z    <- estm/sem
  p    <- 2*pnorm(abs(z), lower.tail = FALSE)
  Q    <- sum(w*(est - estm)^2, na.rm=TRUE)
  k    <- sum(is.finite(w) & w>0)
  Qdf  <- max(k-1, 0)
  Qp   <- if (Qdf > 0) 1 - pchisq(Q, df=Qdf) else NA_real_
  I2   <- if (Qdf > 0 && is.finite(Q) && Q>0) max((Q - Qdf)/Q, 0) * 100 else 0
  list(est=estm, se=sem, z=z, p=p, Q=Q, Qdf=Qdf, Qp=Qp, I2=I2)
}
col_or_na <- function(D, nm) if (nm %in% names(D)) D[[nm]] else rep(NA_real_, nrow(D))
safe_sign <- function(x) ifelse(x > 0, 1, ifelse(x < 0, -1, 0))
mode_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


use_q <- tolower(opt$sig_metric) == "q"
pthr  <- opt$p_level
qthr  <- opt$fdr_level
F_thr <- opt$f_threshold
steiger_thr <- opt$steiger_min_prop

## ========= branch A: mediation rows (your original meta) =========
MED <- DT[analysis %in% c("med","mediation","")]
META_MED <- data.table()

if (nrow(MED)) {
  by_keys <- c("gene_name","snp","chr","pos","variant_id","outcome","mediator")
  miss_id <- setdiff(by_keys, names(MED))
  if (length(miss_id)) stop("Mediation rows missing required columns: ", paste(miss_id, collapse=", "))
  
  META_MED <- MED[, {
    ia <- meta_one(col_or_na(.SD,"indirect"), col_or_na(.SD,"indirect_se"))
    id <- meta_one(col_or_na(.SD,"direct"),   col_or_na(.SD,"direct_se"))
    ba <- meta_one(col_or_na(.SD,"alpha"),    col_or_na(.SD,"alpha_se"))
    bb <- meta_one(col_or_na(.SD,"beta"),     col_or_na(.SD,"beta_se"))
    ct <- meta_one(col_or_na(.SD,"c_total"),  col_or_na(.SD,"c_total_se"))
    iv <- meta_one(col_or_na(.SD,"iv_theta"), col_or_na(.SD,"iv_theta_se"))
    
    total_vec     <- if ("total"    %in% names(.SD)) col_or_na(.SD,"total")    else col_or_na(.SD,"direct")+col_or_na(.SD,"indirect")
    prop_vec      <- if ("prop_med" %in% names(.SD)) col_or_na(.SD,"prop_med") else 100 * col_or_na(.SD,"indirect") / (col_or_na(.SD,"direct")+col_or_na(.SD,"indirect"))
    sign_coh_vec  <- safe_sign(col_or_na(.SD,"indirect")) == safe_sign(col_or_na(.SD,"c_total"))
    sign_coh_prop <- mean(sign_coh_vec, na.rm=TRUE)
    
    F_alpha_vec   <- col_or_na(.SD,"F_alpha")
    R2_alpha_vec  <- col_or_na(.SD,"R2_alpha")
    F_total_vec   <- col_or_na(.SD,"F_total")
    R2_total_vec  <- col_or_na(.SD,"R2_total")
    steiger_vec   <- if ("steiger_ok" %in% names(.SD)) .SD[["steiger_ok"]] else rep(NA, .N)
    iv_ok_vec     <- if ("iv_ok"      %in% names(.SD)) .SD[["iv_ok"]]      else rep(NA, .N)
    
    n_folds          <- .N
    n_alpha_total    <- sum(col_or_na(.SD,"n_alpha"),   na.rm=TRUE)
    n_outcome_total  <- sum(col_or_na(.SD,"n_outcome"), na.rm=TRUE)
    
    F_alpha_meta <- if (is.finite(ba$se) && ba$se > 0) (ba$est / ba$se)^2 else NA_real_
    
    total_meta_ab   <- id$est + ia$est
    prop_med_meta   <- if (is.finite(total_meta_ab) && total_meta_ab != 0) 100 * ia$est / total_meta_ab else NA_real_
    delta_c_vs_total_meta <- ct$est - total_meta_ab
    
    list(
      analysis = "med",
      n_folds = n_folds,
      n_alpha_sum = n_alpha_total,
      n_outcome_sum = n_outcome_total,
      alpha_meta   = ba$est, alpha_se_meta   = ba$se, alpha_z_meta   = ba$z, alpha_p_meta   = ba$p,
      alpha_Q      = ba$Q,   alpha_Qdf       = ba$Qdf, alpha_Qp      = ba$Qp, alpha_I2      = ba$I2,
      beta_meta    = bb$est, beta_se_meta    = bb$se, beta_z_meta    = bb$z, beta_p_meta    = bb$p,
      beta_Q       = bb$Q,   beta_Qdf        = bb$Qdf, beta_Qp       = bb$Qp, beta_I2       = bb$I2,
      direct_meta  = id$est, direct_se_meta  = id$se, direct_z_meta  = id$z, direct_p_meta  = id$p,
      direct_Q     = id$Q,   direct_Qdf      = id$Qdf, direct_Qp     = id$Qp, direct_I2     = id$I2,
      indirect_meta= ia$est, indirect_se_meta= ia$se, indirect_z_meta= ia$z, indirect_p_meta= ia$p,
      indirect_Q   = ia$Q,   indirect_Qdf    = ia$Qdf, indirect_Qp   = ia$Qp, indirect_I2   = ia$I2,
      c_total_meta = ct$est, c_total_se_meta = ct$se, c_total_z_meta = ct$z, c_total_p_meta = ct$p,
      c_total_Q    = ct$Q,   c_total_Qdf     = ct$Qdf, c_total_Qp    = ct$Qp, c_total_I2    = ct$I2,
      iv_theta_meta= iv$est, iv_theta_se_meta= iv$se, iv_theta_z_meta= iv$z, iv_theta_p_meta= iv$p,
      iv_theta_Q   = iv$Q,   iv_theta_Qdf    = iv$Qdf, iv_theta_Qp   = iv$Qp, iv_theta_I2   = iv$I2,
      total_meta_ab         = total_meta_ab,
      prop_med_meta         = prop_med_meta,
      delta_c_vs_total_meta = delta_c_vs_total_meta,
      total_mean      = mean(total_vec, na.rm=TRUE),
      total_median    = median(total_vec, na.rm=TRUE),
      prop_med_mean   = mean(prop_vec, na.rm=TRUE),
      prop_med_median = median(prop_vec, na.rm=TRUE),
      prop_med_min    = suppressWarnings(min(prop_vec, na.rm=TRUE)),
      prop_med_max    = suppressWarnings(max(prop_vec, na.rm=TRUE)),
      sign_coherence_prop = sign_coh_prop,
      F_alpha_meta = F_alpha_meta,
      F_alpha_mean = mean(F_alpha_vec, na.rm=TRUE),
      F_alpha_median = median(F_alpha_vec, na.rm=TRUE),
      F_alpha_min  = suppressWarnings(min(F_alpha_vec, na.rm=TRUE)),
      R2_alpha_mean= mean(R2_alpha_vec, na.rm=TRUE),
      R2_alpha_median= median(R2_alpha_vec, na.rm=TRUE),
      R2_alpha_min = suppressWarnings(min(R2_alpha_vec, na.rm=TRUE)),
      F_total_mean = mean(F_total_vec, na.rm=TRUE),
      F_total_median= median(F_total_vec, na.rm=TRUE),
      F_total_min  = suppressWarnings(min(F_total_vec, na.rm=TRUE)),
      R2_total_mean= mean(R2_total_vec, na.rm=TRUE),
      R2_total_median= median(R2_total_vec, na.rm=TRUE),
      R2_total_min = suppressWarnings(min(R2_total_vec, na.rm=TRUE)),
      steiger_prop = mean(steiger_vec == TRUE, na.rm=TRUE),
      steiger_all  = all(steiger_vec == TRUE, na.rm=TRUE),
      iv_ok_prop   = mean(iv_ok_vec == TRUE, na.rm=TRUE)
    )
  }, by = by_keys]
  
  # BH calls (unchanged)
  META_MED[, q_alpha_meta    := p.adjust(alpha_p_meta,    method="BH"), by=.(outcome)]
  META_MED[, q_beta_meta     := p.adjust(beta_p_meta,     method="BH"), by=.(outcome)]
  META_MED[, q_direct_meta   := p.adjust(direct_p_meta,   method="BH"), by=.(outcome)]
  META_MED[, q_indirect_meta := p.adjust(indirect_p_meta, method="BH"), by=.(outcome)]
  META_MED[, q_c_total_meta  := p.adjust(c_total_p_meta,  method="BH"), by=.(outcome)]
  META_MED[, q_iv_theta_meta := p.adjust(iv_theta_p_meta, method="BH"), by=.(outcome)]
  
  META_MED[, `:=`(
    sig_alpha_meta    = if (use_q) q_alpha_meta    < qthr else alpha_p_meta    < pthr,
    sig_beta_meta     = if (use_q) q_beta_meta     < qthr else beta_p_meta     < pthr,
    sig_direct_meta   = if (use_q) q_direct_meta   < qthr else direct_p_meta   < pthr,
    sig_indirect_meta = if (use_q) q_indirect_meta < qthr else indirect_p_meta < pthr,
    sig_c_total_meta  = if (use_q) q_c_total_meta  < qthr else c_total_p_meta  < pthr,
    sig_iv_theta_meta = if (use_q) q_iv_theta_meta < qthr else iv_theta_p_meta < pthr
  )]
  
  classify <- function(has_med, has_dir, has_assoc, ind_est, dir_est) {
    if (isTRUE(has_med) && !isTRUE(has_dir)) return("Full Mediation")
    if (isTRUE(has_med) &&  isTRUE(has_dir))
      return(ifelse((ind_est < 0 & dir_est > 0) | (ind_est > 0 & dir_est < 0),
                    "Partial Mediation - Reverse", "Partial Mediation"))
    if (!isTRUE(has_med) && isTRUE(has_dir) && isTRUE(has_assoc))
      return(ifelse(dir_est > 0, "Pleiotropy - Positive", "Pleiotropy - Negative"))
    if (!isTRUE(has_med) && isTRUE(has_dir) && !isTRUE(has_assoc))
      return(ifelse(dir_est > 0, "Only Direct Effect - Positive", "Only Direct Effect - Negative"))
    "No Effect"
  }
  META_MED[, conclusion_meta :=
             mapply(classify, sig_indirect_meta, sig_direct_meta, sig_alpha_meta,
                    indirect_meta, direct_meta)]
  META_MED[, sig_rule_meta := if (use_q) sprintf("FDR<%.2f", qthr) else sprintf("p<%.3f", pthr)]
  
  META_MED[, strong_iv_meta    := is.finite(alpha_se_meta) & is.finite(alpha_meta) & ((alpha_meta/alpha_se_meta)^2 >= F_thr)]
  META_MED[, steiger_ok_loose  := is.finite(steiger_prop) & (steiger_prop >= steiger_thr)]
  META_MED[, steiger_ok_strict := isTRUE(steiger_all)]
  META_MED[, iv_ok_meta_loose  := strong_iv_meta & sig_alpha_meta & !sig_direct_meta & steiger_ok_loose]
  META_MED[, iv_ok_meta_strict := strong_iv_meta & sig_alpha_meta & !sig_direct_meta & steiger_ok_strict]
}

## ========= branch B: MVMR rows =========
MVMR <- DT[analysis == "mvmr"]
META_MVMR <- data.table()

if (nrow(MVMR)) {
  keys <- c("anchor","outcome","exp1","exp2")
  miss <- setdiff(keys, names(MVMR))
  if (length(miss)) stop("MVMR rows missing required columns: ", paste(miss, collapse=", "))
  
  META_MVMR <- MVMR[, {
    ## meta of direct effects for exp1 and exp2 separately
    m1 <- meta_one(col_or_na(.SD,"beta_exp1_dir"), col_or_na(.SD,"se_exp1_dir"))
    m2 <- meta_one(col_or_na(.SD,"beta_exp2_dir"), col_or_na(.SD,"se_exp2_dir"))
    
    ## summarize conditional F across folds
    # Robust aggregators (avoid Inf/-Inf when all NA)
    safe_mean   <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm=TRUE)
    safe_median <- function(v) if (all(is.na(v))) NA_real_ else median(v, na.rm=TRUE)
    safe_min    <- function(v) if (all(is.na(v))) NA_real_ else suppressWarnings(min(v, na.rm=TRUE))
    safe_max    <- function(v) if (all(is.na(v))) NA_real_ else suppressWarnings(max(v, na.rm=TRUE))
    
    # NOTE: all column names have already been lower-cased above
    F1v <- col_or_na(.SD, "f_cond_exp1")
    F2v <- col_or_na(.SD, "f_cond_exp2")
    Qv  <- col_or_na(.SD, "q_mvmr")
    Qpv <- col_or_na(.SD, "qp_mvmr")
    NW  <- col_or_na(.SD, "n_snps")
    WV  <- .SD[["weak_iv_gate"]]  # already lower-cased; may be absent -> NULL
    
    list(
      analysis       = "mvmr",
      n_folds        = .N,
      n_snps_mean    = safe_mean(NW),
      n_snps_min     = safe_min(NW),
      n_snps_max     = safe_max(NW),
      F1_mean        = safe_mean(F1v),
      F1_median      = safe_median(F1v),
      F1_min         = safe_min(F1v),
      F1_max         = safe_max(F1v),
      F1_prop_ge10   = if (all(is.na(F1v))) NA_real_ else mean(F1v >= 10, na.rm=TRUE),
      F2_mean        = safe_mean(F2v),
      F2_median      = safe_median(F2v),
      F2_min         = safe_min(F2v),
      F2_max         = safe_max(F2v),
      F2_prop_ge10   = if (all(is.na(F2v))) NA_real_ else mean(F2v >= 10, na.rm=TRUE),
      
      # meta for direct effects (unchanged)
      beta1_meta     = m1$est, se1_meta = m1$se, z1_meta = m1$z, p1_meta = m1$p, Q1 = m1$Q, Q1df = m1$Qdf, Q1p = m1$Qp, I2_1 = m1$I2,
      beta2_meta     = m2$est, se2_meta = m2$se, z2_meta = m2$z, p2_meta = m2$p, Q2 = m2$Q, Q2df = m2$Qdf, Q2p = m2$Qp, I2_2 = m2$I2,
      
      # pleiotropy summaries
      Q_mvmr_mean    = safe_mean(Qv),
      Qp_mvmr_mean   = safe_mean(Qpv),
      Qp_mvmr_prop5  = if (all(is.na(Qpv))) NA_real_ else mean(Qpv < 0.05, na.rm=TRUE),
      
      # weak IV bookkeeping across folds
      weak_iv_prop   = if (is.null(WV)) NA_real_ else mean(WV == TRUE, na.rm=TRUE),
      weak_iv_any    = if (is.null(WV)) NA else any(WV == TRUE, na.rm=TRUE)
    )
    
    nsv <- col_or_na(.SD,"n_snps")
    
    ## weak-IV gating summaries (if present)
    wgv <- if ("weak_iv_gate" %in% names(.SD)) .SD[["weak_iv_gate"]] else rep(NA, .N)
    weak_prop <- mean(wgv == TRUE, na.rm = TRUE)
    weak_any  <- any(wgv == TRUE,  na.rm = TRUE)
    weak_all  <- all(wgv == TRUE,  na.rm = TRUE)
    
    ## fold-level conclusion mode (if present)
    cons_vec  <- if ("conclusion_simple" %in% names(.SD)) as.character(.SD[["conclusion_simple"]]) else rep(NA_character_, .N)
    cons_mode <- mode_chr(cons_vec)
    cons_prop <- mean(cons_vec == cons_mode, na.rm = TRUE)
    
    list(
      analysis       = "mvmr",
      n_folds        = .N,
      
      ## ---- fold-like columns: provide aggregates under familiar names ----
      n_snps_median  = suppressWarnings(median(nsv, na.rm=TRUE)),
      n_snps_mean    = mean(nsv, na.rm=TRUE),
      n_snps_min     = suppressWarnings(min(nsv, na.rm=TRUE)),
      n_snps_max     = suppressWarnings(max(nsv, na.rm=TRUE)),
      
      F1_mean        = mean(F1v, na.rm=TRUE),
      F1_median      = median(F1v, na.rm=TRUE),
      F1_min         = suppressWarnings(min(F1v, na.rm=TRUE)),
      F1_max         = suppressWarnings(max(F1v, na.rm=TRUE)),
      F1_prop_ge10   = mean(F1v >= 10, na.rm=TRUE),
      
      F2_mean        = mean(F2v, na.rm=TRUE),
      F2_median      = median(F2v, na.rm=TRUE),
      F2_min         = suppressWarnings(min(F2v, na.rm=TRUE)),
      F2_max         = suppressWarnings(max(F2v, na.rm=TRUE)),
      F2_prop_ge10   = mean(F2v >= 10, na.rm=TRUE),
      
      ## meta of per-exposure effects
      beta1_meta     = m1$est, se1_meta = m1$se, z1_meta = m1$z, p1_meta = m1$p, Q1 = m1$Q, Q1df = m1$Qdf, Q1p = m1$Qp, I2_1 = m1$I2,
      beta2_meta     = m2$est, se2_meta = m2$se, z2_meta = m2$z, p2_meta = m2$p, Q2 = m2$Q, Q2df = m2$Qdf, Q2p = m2$Qp, I2_2 = m2$I2,
      
      ## heterogeneity (Q) summary across folds
      Qp_mvmr_mean   = mean(Qpv, na.rm=TRUE),
      Qp_mvmr_min    = suppressWarnings(min(Qpv, na.rm=TRUE)),
      Qp_mvmr_prop5  = mean(Qpv < 0.05, na.rm=TRUE),
      
      ## weak-IV gate summaries
      weak_iv_prop   = weak_prop,
      weak_iv_any    = weak_any,
      weak_iv_all    = weak_all,
      
      ## fold conclusion mode (diagnostic)
      conclusion_simple_fold_mode      = cons_mode,
      conclusion_simple_fold_mode_prop = cons_prop
    )
  }, by = keys]
  
  
  # BH by outcome for each exposure's direct effect
  META_MVMR[, q1_meta := p.adjust(p1_meta, method="BH"), by=.(outcome)]
  META_MVMR[, q2_meta := p.adjust(p2_meta, method="BH"), by=.(outcome)]
  
  META_MVMR[, `:=`(
    sig1_meta = if (use_q) q1_meta < qthr else p1_meta < pthr,
    sig2_meta = if (use_q) q2_meta < qthr else p2_meta < pthr
  )]
  
  ## ----- simple meta-level conclusion (same 4 buckets you already use) -----
  META_MVMR[, conclusion_simple_meta := {
    s1 <- is.finite(F1_median) && F1_median >= 10
    s2 <- is.finite(F2_median) && F2_median >= 10
    pleio <- is.finite(Qp_mvmr_mean) && Qp_mvmr_mean < 0.05
    weak_any <- isTRUE(weak_iv_any)
    
    if (s1 && s2 && sig1_meta && sig2_meta) {
      "Bidirectional"
    } else if (s1 && sig1_meta && !sig2_meta) {
      "True causal"
    } else if (pleio || weak_any || (sig1_meta && !s1) || (sig2_meta && !s2)) {
      "Correlated pleiotropy"
    } else {
      "Not causal"
    }
  }, by = .(anchor, outcome)]
  
  
  ## ----- add fold-like aliases so downstream code finds expected names -----
  META_MVMR[, `:=`(
    # keep your meta summaries…
    # …and ALSO provide the fold-lookalike columns as sensible aggregates:
    conclusion_simple = conclusion_simple_meta,
    n_snps            = n_snps_median,
    F_cond_exp1       = F1_median,
    F_cond_exp2       = F2_median,
    beta_exp1_dir     = beta1_meta,
    se_exp1_dir       = se1_meta,
    p_exp1_dir        = p1_meta,
    beta_exp2_dir     = beta2_meta,
    se_exp2_dir       = se2_meta,
    p_exp2_dir        = p2_meta,
    Qp_mvmr           = Qp_mvmr_mean,
    weak_iv_gate      = weak_iv_any
  )]
  
}

## ----- Simple meta-level conclusions (same 4 buckets) -----

# Mediation/univariable (per mediator/outcome) -- reduce vectors to scalars per group
if (nrow(META_MED)) {
  META_MED[, conclusion_simple_meta := {
    # Use the strict IV validity flag you already computed upstream
    iv_ok_grp  <- any(iv_ok_meta_strict == TRUE, na.rm = TRUE)
    wald_grp   <- any(sig_iv_theta_meta == TRUE, na.rm = TRUE)
    assoc_grp  <- any(sig_c_total_meta  == TRUE, na.rm = TRUE)
    direct_grp <- any(sig_direct_meta   == TRUE, na.rm = TRUE)
    
    if (iv_ok_grp && wald_grp) {
      "True causal"
    } else if (assoc_grp && (!iv_ok_grp || direct_grp)) {
      "Correlated pleiotropy"
    } else {
      "Not causal"
    }
  }, by = .(mediator, outcome)]
}

## ========= combine & write =========
OUT <- rbindlist(list(META_MED, META_MVMR), fill=TRUE)

# Light column ordering: keep keys & useful fields in front
front <- c("analysis","gene_name","snp","chr","pos","variant_id","mediator","anchor","exp1","exp2","outcome")
setcolorder(OUT, c(intersect(front, names(OUT)), setdiff(names(OUT), front)))

fwrite(OUT, file = opt$out, sep = "\t")
cat("[meta_mediation] WROTE", opt$out, "rows=", nrow(OUT), "\n")
