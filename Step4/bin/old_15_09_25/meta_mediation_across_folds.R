#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--inputs-list", type="character", help="One per-fold mediation TSV per line"),
  make_option("--sig_metric",  type="character", default="q", help="p or q for meta-level calls"),
  make_option("--fdr_level",   type="double",    default=0.10),
  make_option("--p_level",     type="double",    default=0.05),
  # IV meta options
  make_option("--f_threshold",       type="double", default=10,    help="F-stat threshold for 'strong' IV"),
  make_option("--steiger_min_prop",  type="double", default=0.60,  help="Min fraction of folds with Steiger OK"),
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

# Required identifiers
by_keys <- c("gene_name","snp","chr","pos","variant_id","outcome","mediator")
miss_id <- setdiff(by_keys, names(DT))
if (length(miss_id)) stop("Missing required columns in per-fold inputs: ", paste(miss_id, collapse=", "))

# Helpers
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

use_q <- tolower(opt$sig_metric) == "q"
pthr  <- opt$p_level
qthr  <- opt$fdr_level
F_thr <- opt$f_threshold
steiger_thr <- opt$steiger_min_prop

META <- DT[, {
  # --- meta core effects ---
  ia <- meta_one(col_or_na(.SD,"indirect"), col_or_na(.SD,"indirect_se"))
  id <- meta_one(col_or_na(.SD,"direct"),   col_or_na(.SD,"direct_se"))
  ba <- meta_one(col_or_na(.SD,"alpha"),    col_or_na(.SD,"alpha_se"))
  bb <- meta_one(col_or_na(.SD,"beta"),     col_or_na(.SD,"beta_se"))
  ct <- meta_one(col_or_na(.SD,"c_total"),  col_or_na(.SD,"c_total_se"))
  iv <- meta_one(col_or_na(.SD,"iv_theta"), col_or_na(.SD,"iv_theta_se"))
  
  # --- mediation-style fold summaries we want to KEEP in meta ---
  total_vec     <- if ("total"    %in% names(.SD)) col_or_na(.SD,"total")    else col_or_na(.SD,"direct")+col_or_na(.SD,"indirect")
  prop_vec      <- if ("prop_med" %in% names(.SD)) col_or_na(.SD,"prop_med") else 100 * col_or_na(.SD,"indirect") / (col_or_na(.SD,"direct")+col_or_na(.SD,"indirect"))
  # sign coherence between indirect and c_total per fold
  sign_coh_vec  <- safe_sign(col_or_na(.SD,"indirect")) == safe_sign(col_or_na(.SD,"c_total"))
  sign_coh_prop <- mean(sign_coh_vec, na.rm=TRUE)
  
  # IV diagnostics (fold-level)
  F_alpha_vec   <- col_or_na(.SD,"F_alpha")
  R2_alpha_vec  <- col_or_na(.SD,"R2_alpha")
  F_total_vec   <- col_or_na(.SD,"F_total")
  R2_total_vec  <- col_or_na(.SD,"R2_total")
  steiger_vec   <- if ("steiger_ok" %in% names(.SD)) .SD[["steiger_ok"]] else rep(NA, .N)
  iv_ok_vec     <- if ("iv_ok"      %in% names(.SD)) .SD[["iv_ok"]]      else rep(NA, .N)
  
  n_folds          <- .N
  n_alpha_total    <- sum(col_or_na(.SD,"n_alpha"),   na.rm=TRUE)
  n_outcome_total  <- sum(col_or_na(.SD,"n_outcome"), na.rm=TRUE)
  
  # meta F for alpha
  F_alpha_meta <- if (is.finite(ba$se) && ba$se > 0) (ba$est / ba$se)^2 else NA_real_
  
  # mediation-style META aggregates
  total_meta_ab   <- id$est + ia$est
  prop_med_meta   <- if (is.finite(total_meta_ab) && total_meta_ab != 0) 100 * ia$est / total_meta_ab else NA_real_
  delta_c_vs_total_meta <- ct$est - total_meta_ab
  
  list(
    n_folds = n_folds,
    n_alpha_sum = n_alpha_total,
    n_outcome_sum = n_outcome_total,
    
    # meta estimates (paths)
    alpha_meta   = ba$est, alpha_se_meta   = ba$se, alpha_z_meta   = ba$z, alpha_p_meta   = ba$p,
    alpha_Q      = ba$Q,   alpha_Qdf       = ba$Qdf, alpha_Qp      = ba$Qp, alpha_I2      = ba$I2,
    
    beta_meta    = bb$est, beta_se_meta    = bb$se, beta_z_meta    = bb$z, beta_p_meta    = bb$p,
    beta_Q       = bb$Q,   beta_Qdf        = bb$Qdf, beta_Qp       = bb$Qp, beta_I2       = bb$I2,
    
    direct_meta  = id$est, direct_se_meta  = id$se, direct_z_meta  = id$z, direct_p_meta  = id$p,
    direct_Q     = id$Q,   direct_Qdf      = id$Qdf, direct_Qp     = id$Qp, direct_I2     = id$I2,
    
    indirect_meta= ia$est, indirect_se_meta= ia$se, indirect_z_meta= ia$z, indirect_p_meta= ia$p,
    indirect_Q   = ia$Q,   indirect_Qdf    = ia$Qdf, indirect_Qp   = ia$Qp, indirect_I2   = ia$I2,
    
    # total effect from Y~G (no M)
    c_total_meta = ct$est, c_total_se_meta = ct$se, c_total_z_meta = ct$z, c_total_p_meta = ct$p,
    c_total_Q    = ct$Q,   c_total_Qdf     = ct$Qdf, c_total_Qp    = ct$Qp, c_total_I2    = ct$I2,
    
    # IV Wald ratio
    iv_theta_meta= iv$est, iv_theta_se_meta= iv$se, iv_theta_z_meta= iv$z, iv_theta_p_meta= iv$p,
    iv_theta_Q   = iv$Q,   iv_theta_Qdf    = iv$Qdf, iv_theta_Qp   = iv$Qp, iv_theta_I2   = iv$I2,
    
    # --- mediation-style meta add-backs ---
    total_meta_ab         = total_meta_ab,          # = direct_meta + indirect_meta
    prop_med_meta         = prop_med_meta,          # %
    delta_c_vs_total_meta = delta_c_vs_total_meta,  # should be ~0 if models align
    
    # fold summaries to keep in meta output
    total_mean      = mean(total_vec, na.rm=TRUE),
    total_median    = median(total_vec, na.rm=TRUE),
    prop_med_mean   = mean(prop_vec, na.rm=TRUE),
    prop_med_median = median(prop_vec, na.rm=TRUE),
    prop_med_min    = suppressWarnings(min(prop_vec, na.rm=TRUE)),
    prop_med_max    = suppressWarnings(max(prop_vec, na.rm=TRUE)),
    sign_coherence_prop = sign_coh_prop,
    
    # IV summaries
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
    iv_ok_prop   = mean(iv_ok_vec == TRUE, na.rm=TRUE),
    
    F_alpha_meta = F_alpha_meta
  )
}, by = by_keys]

# clean up -Inf from min()
for (nm in c("F_alpha_min","R2_alpha_min","F_total_min","R2_total_min","prop_med_min","prop_med_max")) {
  if (nm %in% names(META)) set(META, i=which(!is.finite(META[[nm]])), j=nm, value=NA_real_)
}

# Multiple testing by outcome
META[, q_alpha_meta    := p.adjust(alpha_p_meta,    method="BH"), by=.(outcome)]
META[, q_beta_meta     := p.adjust(beta_p_meta,     method="BH"), by=.(outcome)]
META[, q_direct_meta   := p.adjust(direct_p_meta,   method="BH"), by=.(outcome)]
META[, q_indirect_meta := p.adjust(indirect_p_meta, method="BH"), by=.(outcome)]
META[, q_c_total_meta  := p.adjust(c_total_p_meta,  method="BH"), by=.(outcome)]
META[, q_iv_theta_meta := p.adjust(iv_theta_p_meta, method="BH"), by=.(outcome)]

META[, `:=`(
  sig_alpha_meta    = if (use_q) q_alpha_meta    < qthr else alpha_p_meta    < pthr,
  sig_beta_meta     = if (use_q) q_beta_meta     < qthr else beta_p_meta     < pthr,
  sig_direct_meta   = if (use_q) q_direct_meta   < qthr else direct_p_meta   < pthr,
  sig_indirect_meta = if (use_q) q_indirect_meta < qthr else indirect_p_meta < pthr,
  sig_c_total_meta  = if (use_q) q_c_total_meta  < qthr else c_total_p_meta  < pthr,
  sig_iv_theta_meta = if (use_q) q_iv_theta_meta < qthr else iv_theta_p_meta < pthr
)]

# Narrative conclusion (meta)
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
META[, conclusion_meta :=
       mapply(classify, sig_indirect_meta, sig_direct_meta, sig_alpha_meta,
              indirect_meta, direct_meta)]
META[, sig_rule_meta := if (use_q) sprintf("FDR<%.2f", qthr) else sprintf("p<%.3f", pthr)]

# Meta IV flags
META[, strong_iv_meta    := is.finite(F_alpha_meta) & (F_alpha_meta >= F_thr)]
META[, strong_iv_any_fold:= is.finite(F_alpha_min)  & (F_alpha_min  >= F_thr)]
META[, steiger_ok_loose  := is.finite(steiger_prop) & (steiger_prop >= steiger_thr)]
META[, steiger_ok_strict := isTRUE(steiger_all)]
META[, iv_ok_meta_loose  := strong_iv_meta & sig_alpha_meta & !sig_direct_meta & steiger_ok_loose]
META[, iv_ok_meta_strict := strong_iv_meta & sig_alpha_meta & !sig_direct_meta & steiger_ok_strict]

# Column order
col_order <- c(by_keys,
               "n_folds","n_alpha_sum","n_outcome_sum",
               # path metas
               "alpha_meta","alpha_se_meta","alpha_z_meta","alpha_p_meta","q_alpha_meta","alpha_Q","alpha_Qdf","alpha_Qp","alpha_I2",
               "beta_meta","beta_se_meta","beta_z_meta","beta_p_meta","q_beta_meta","beta_Q","beta_Qdf","beta_Qp","beta_I2",
               "direct_meta","direct_se_meta","direct_z_meta","direct_p_meta","q_direct_meta","direct_Q","direct_Qdf","direct_Qp","direct_I2",
               "indirect_meta","indirect_se_meta","indirect_z_meta","indirect_p_meta","q_indirect_meta","indirect_Q","indirect_Qdf","indirect_Qp","indirect_I2",
               "c_total_meta","c_total_se_meta","c_total_z_meta","c_total_p_meta","q_c_total_meta","c_total_Q","c_total_Qdf","c_total_Qp","c_total_I2",
               "iv_theta_meta","iv_theta_se_meta","iv_theta_z_meta","iv_theta_p_meta","q_iv_theta_meta","iv_theta_Q","iv_theta_Qdf","iv_theta_Qp","iv_theta_I2",
               # mediation add-backs
               "total_meta_ab","prop_med_meta","delta_c_vs_total_meta",
               "total_mean","total_median",
               "prop_med_mean","prop_med_median","prop_med_min","prop_med_max",
               "sign_coherence_prop",
               # IV summaries / gates
               "F_alpha_meta","F_alpha_mean","F_alpha_median","F_alpha_min",
               "R2_alpha_mean","R2_alpha_median","R2_alpha_min",
               "F_total_mean","F_total_median","F_total_min",
               "R2_total_mean","R2_total_median","R2_total_min",
               "steiger_prop","steiger_all","iv_ok_prop",
               "strong_iv_meta","strong_iv_any_fold","steiger_ok_loose","steiger_ok_strict",
               "iv_ok_meta_loose","iv_ok_meta_strict",
               # calls
               "sig_alpha_meta","sig_beta_meta","sig_direct_meta","sig_indirect_meta","sig_c_total_meta","sig_iv_theta_meta",
               "conclusion_meta","sig_rule_meta"
)
setcolorder(META, intersect(col_order, names(META)))

fwrite(META, file = opt$out, sep = "\t")
cat("[meta_mediation] WROTE", opt$out, "rows=", nrow(META), "\n")
