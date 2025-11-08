#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------------- CLI ----------------------
opt <- parse_args(OptionParser(option_list=list(
  make_option("--meta",              type="character", help="mediation_meta__<Genus>.tsv"),
  make_option("--out",               type="character", default="mediation_sig_iv"),
  make_option("--q",                 type="double",    default=0.10, help="FDR cutoff for significance"),
  make_option("--Fmin",              type="double",    default=10,   help="Min meta F for strong instrument"),
  make_option("--steiger_min_prop",  type="double",    default=0.60, help="Min fraction of folds Steiger OK"),
  make_option("--sign_coh_min",      type="double",    default=0.60, help="Min fraction of folds with sign(indirect)==sign(c_total)"),
  make_option("--delta_frac_tol",    type="double",    default=0.05, help="|c_total - (direct+indirect)| <= tol * max(|c_total|,|direct+indirect|)"),
  make_option("--require_ivtheta",   action="store_true", default=FALSE,
              help="If set, also require IV Wald ratio meta to be significant"),
  make_option("--hetero_max_I2",     type="double",    default=50,   help="Optional I2 ceiling for Wald ratio (set high to disable)"),
  make_option("--sort_by",           type="character", default="q_iv_theta_meta",
              help="Primary sort column: q_iv_theta_meta|q_indirect_meta|F_alpha_meta")
)))
stopifnot(file.exists(opt$meta))

DT <- fread(opt$meta)
setDT(DT); setnames(DT, tolower(names(DT)))

# ---------------------- Helpers ----------------------
nznum <- function(x) ifelse(is.finite(x), x, NA_real_)
coh_ok <- function(delta, total_ab, c_total, frac_tol=0.05) {
  # |c_total - total_ab| <= frac_tol * max(|c_total|, |total_ab|)
  M <- pmax(abs(total_ab), abs(c_total), 1e-12)
  ok <- abs(delta) <= frac_tol * M
  ifelse(is.finite(ok), ok, FALSE)
}

# ---------------------- Core filters ----------------------
q     <- opt$q
Fmin  <- opt$Fmin
propS <- opt$steiger_min_prop
propC <- opt$sign_coh_min

# 1) Instrument validity (meta-level)
DT[, pass_alpha     := sig_alpha_meta == TRUE]
DT[, pass_strength  := is.finite(f_alpha_meta) & f_alpha_meta >= Fmin]
DT[, pass_exclusion := sig_direct_meta != TRUE]        # want NOT significant
DT[, pass_steiger   := (is.finite(steiger_prop) & steiger_prop >= propS) | (steiger_all == TRUE)]

# 2) Full mediation structure (meta-level)
DT[, pass_mediation := sig_indirect_meta == TRUE & sig_direct_meta != TRUE]

# 3) Coherence checks (recommended)
DT[, delta_ok := coh_ok(delta_c_vs_total_meta, total_meta_ab, c_total_meta, frac_tol = opt$delta_frac_tol)]
DT[, sign_coh_ok := is.finite(sign_coherence_prop) & sign_coherence_prop >= propC]

# 4) Optional Wald ratio significance & reasonable heterogeneity
DT[, pass_wald := sig_iv_theta_meta == TRUE]
DT[, pass_het  := is.finite(iv_theta_i2) & iv_theta_i2 <= opt$hetero_max_I2]

# Composite instrument gate
DT[, iv_valid_loose  := pass_alpha & pass_strength & pass_exclusion & pass_steiger]
DT[, iv_valid_strict := iv_valid_loose & sign_coh_ok & delta_ok]

# Final selection
if (opt$require_ivtheta) {
  SEL <- DT[ iv_valid_strict & pass_mediation & pass_wald & pass_het ]
} else {
  SEL <- DT[ iv_valid_strict & pass_mediation ]
}

# ---------------------- Pretty summary table ----------------------
# Tiers for quick reading
SEL[, tier := fifelse(sig_iv_theta_meta==TRUE & (iv_theta_i2 <= 25 | iv_theta_qp >= 0.10), "gold",
                      fifelse(sig_iv_theta_meta==TRUE, "silver", "bronze"))]

# Pick the most informative columns
keep <- c(
  "gene_name","snp","chr","pos","mediator","outcome",
  # instrument path
  "alpha_meta","alpha_se_meta","q_alpha_meta","F_alpha_meta","steiger_prop","steiger_all",
  # mediator->outcome and mediation
  "beta_meta","beta_se_meta","q_beta_meta",
  "indirect_meta","indirect_se_meta","q_indirect_meta",
  "direct_meta","direct_se_meta","q_direct_meta",
  # totals (two views)
  "c_total_meta","c_total_se_meta","q_c_total_meta",
  "total_meta_ab","prop_med_meta","delta_c_vs_total_meta","sign_coherence_prop",
  # IV estimate + heterogeneity
  "iv_theta_meta","iv_theta_se_meta","q_iv_theta_meta","iv_theta_I2","iv_theta_Qp",
  # fold counts
  "n_folds","n_alpha_sum","n_outcome_sum",
  # gates + tier
  "iv_valid_loose","iv_valid_strict","pass_mediation","tier"
)
keep <- intersect(keep, names(SEL))
OUT <- SEL[, ..keep]

# Sorting preference
#ord <- switch(opt$sort_by,
#              "q_indirect_meta" = order(OUT$q_indirect_meta, OUT$q_alpha_meta, OUT$f_alpha_meta, na.last=TRUE),
#              "f_alpha_meta"    = order(OUT$f_alpha_meta, OUT$q_alpha_meta, OUT$q_indirect_meta, na.last=TRUE),
#              "q_iv_theta_meta" = order(OUT$q_iv_theta_meta, OUT$q_indirect_meta, OUT$f_alpha_meta, na.last=TRUE),
#              order(OUT$q_iv_theta_meta, OUT$q_indirect_meta, OUT$f_alpha_meta, na.last=TRUE)
#)
#OUT <- OUT[ord]

fwrite(OUT, file=opt$out, sep="\t")
cat(sprintf("[filter_full_mediation] wrote %s (rows=%d)\n", opt$out, nrow(OUT)))
