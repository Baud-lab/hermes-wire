#!/usr/bin/env Rscript

## ================================================================
## Feed-forward SEM with guarded preprocessing to avoid singularities
## Default SEM: Genes → Bact → Prev → Beta → Glucose
## Guild: parsed/residualized and, if --include_Guild TRUE, included as
##         Prev -> Guild -> Beta (plus diagnostics).
## - Locks cohort order
## - Reads 'gold' instruments for Bact and Prev; splits Bact-only / Prev-only / Overlap
## - Drops rows with NA Glucose BEFORE anything else (so every test uses identical samples)
## - Pulls dosages from one GDS (for the filtered cohort)
## - Residualizes targets via sommer (same RE structure), IN-PLACE (replaces raw columns)
## - Pre-SEM filters: drop constants, duplicates, high correlations (|r|>=thr), QR-drop linearly dependent GENE columns only
## - Adds lm() + VIF diagnostics on (residualized) mediators/outcome (with and without Guild)
## - Builds SEM (Guild optional); trims by BH-FDR (10%) with feed-forward cascade
## - Mermaid diagram; blavaan validates the final (trimmed) lavaan model
## - Logs: sommer heads (before/after), per-iteration models, drop logs, VIFs, fits
## ================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(readr); library(readxl)
  library(stringr); library(purrr); library(tidyr); library(tibble)
  library(gdsfmt); library(SNPRelate)
  library(lavaan); library(semTools); library(blavaan)
  library(sommer); library(car)
  library(igraph)   # used for path enumeration when computing indirect/total effects
})

## -------- logging --------
options(warn=1)
msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%F %T"), paste0(...)))
  try(flush(stdout()), silent=TRUE); try(flush(stderr()), silent=TRUE)
}

## -------- helpers --------

## ---- persistent labels (global) ----
label_map <- new.env(parent = emptyenv())
label_ctr <- new.env(parent = emptyenv())

label_for <- function(lhs, rhs, pref) {
  key <- paste(lhs, rhs, sep="~")
  if (exists(key, envir = label_map, inherits = FALSE)) {
    return(get(key, envir = label_map, inherits = FALSE))
  }
  n <- get0(pref, envir = label_ctr, ifnotfound = 0L) + 1L
  assign(pref, n, envir = label_ctr)
  lab <- sprintf("%s%d", pref, n)
  assign(key, lab, envir = label_map)
  lab
}

norm_id <- function(x) base::sub("_.*$", "", as.character(x))

uniq    <- function(x) unique(x[!is.na(x)])

pick_col <- function(M, key) {
  cn <- colnames(M)
  exact <- which(cn == key); if (length(exact)) return(cn[exact[1]])
  st <- which(startsWith(cn, key)); if (length(st)) return(cn[st[1]])
  NA_character_
}

stop_if_missing <- function(who, have, need) {
  missing <- setdiff(need, have)
  if (length(missing)) stop(sprintf("%s is missing %d cohort IDs. Examples: %s",
                                    who, length(missing), paste(head(missing, 10), collapse=", ")))
}

.align_by_key <- function(df, row_key, value_cols, id_cols) {
  vc <- intersect(value_cols, names(df)); if (!length(vc)) return(NULL)
  ic <- intersect(id_cols,    names(df)); if (!length(ic)) return(NULL)
  val <- df[[ vc[1] ]]
  ids <- norm_id(df[[ ic[1] ]])
  ord <- match(norm_id(row_key), ids)
  if (anyNA(ord)) return(NULL)
  as.numeric(val[ord])
}

## robust residual vector from sommer::mmer (works with 4.x; aligns to row_key)
## ---- align a numeric vector to row_key if it has names ----
.align_num_to_key <- function(v, row_key) {
  if (is.null(v)) return(NULL)
  v <- as.numeric(v)
  if (!is.null(names(v))) {
    ord <- match(row_key, as.character(names(v)))
    return(v[ord])
  }
  if (length(v) == length(row_key)) return(v)
  NULL
}

## ---- robust fitted vector from sommer::mmer (prefer object slots) ----
fitted_vec_from_mmer <- function(fit, row_key, trait = NULL) {
  for (nm in c("fitted", "fitted.y", "fitted.values")) {
    fv <- tryCatch(fit[[nm]], error = function(e) NULL)
    fv <- .align_num_to_key(fv, row_key)
    if (is.numeric(fv)) return(as.numeric(fv))
  }
  # generic fallbacks
  fv <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  if (is.numeric(fv) && length(fv)) {
    fv2 <- .align_num_to_key(fv, row_key); if (is.numeric(fv2)) return(fv2)
    if (length(fv) == length(row_key)) return(as.numeric(fv))
  }
  # older sommer: fitted.mmer() -> $pvals
  ff <- tryCatch(sommer:::fitted.mmer(fit), error = function(e) NULL)
  if (is.list(ff) && is.data.frame(ff$pvals)) {
    pv <- ff$pvals
    if (!is.null(trait) && "trait" %in% names(pv)) pv <- pv[pv$trait == trait, , drop = FALSE]
    return(.align_by_key(pv, row_key,
                         value_cols = c("predicted.value","fitted","fitted.value"),
                         id_cols    = c("Row","row","sample","Sample","id","ID")))
  } else if (is.data.frame(ff)) {
    return(.align_by_key(ff, row_key,
                         value_cols = c("predicted.value","fitted","fitted.value"),
                         id_cols    = c("Row","row","sample","Sample","id","ID")))
  }
  NULL
}

log_lavaan <- function(fit, iter, stage, tag = "") {
  hdr <- sprintf("[SEM] === lavaan summary (iter %d | stage: %s%s) ===",
                 iter, stage, if (nzchar(tag)) paste0(" | ", tag) else "")
  msg(hdr)
  out <- utils::capture.output(
    summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE),
    modindices(fit, sort. = TRUE, minimum.value = 3.84)
  )
  # print to console
  writeLines(out)
  # and save to disk
  fn <- file.path(opt$outdir,
                  sprintf("lavaan_iter_%02d_%s%s_summary.txt",
                          iter, stage, if (nzchar(tag)) paste0("_", tag) else ""))
  writeLines(c(hdr, out), con = fn)
}


## ---- robust residual vector (prefer object slot), aligned to row_key ----
resid_vec_from_mmer <- function(fit, y, row_key, trait = NULL) {
  # 0) object slot (works in sommer 4.x)
  v <- tryCatch(fit$residuals, error = function(e) NULL)
  v <- .align_num_to_key(v, row_key)
  if (is.numeric(v) && length(v) == length(y) && all(is.finite(v))) return(as.numeric(v))
  
  # 1) stats::residuals fallback
  v <- tryCatch(stats::residuals(fit), error = function(e) NULL)
  if (is.numeric(v) && length(v) == length(y) && all(is.finite(v))) return(as.numeric(v))
  
  # 2) sommer:::residuals.mmer() (data.frame shapes)
  rr <- tryCatch(sommer:::residuals.mmer(fit), error = function(e) NULL)
  if (is.list(rr) && !is.null(rr$residuals)) {
    v <- .align_by_key(rr$residuals, row_key,
                       value_cols = c("residual","residuals","eBLUP","e"),
                       id_cols    = c("Row","row","sample","Sample","id","ID"))
    if (is.numeric(v) && length(v) == length(y) && all(is.finite(v))) return(as.numeric(v))
  } else if (is.data.frame(rr)) {
    v <- .align_by_key(rr, row_key,
                       value_cols = c("residual","residuals","eBLUP","e"),
                       id_cols    = c("Row","row","sample","Sample","id","ID"))
    if (is.numeric(v) && length(v) == length(y) && all(is.finite(v))) return(as.numeric(v))
  }
  
  # 3) last resort: compute e = y - fitted
  fv <- fitted_vec_from_mmer(fit, row_key, trait)
  if (is.numeric(fv) && length(fv) == length(y)) {
    v <- as.numeric(y) - as.numeric(fv)
    if (all(is.finite(v))) return(v)
  }
  
  # 4) absolute fallback
  warning("Residual extraction failed; falling back to centered y for this trait.")
  as.numeric(y) - mean(as.numeric(y), na.rm = TRUE)
}


## Fixed part Xβ as a numeric vector (works even if beta.hat shape varies)
fixed_part_vec_from_mmer <- function(fit, formula_str, data) {
  X <- stats::model.matrix(as.formula(formula_str), data = data)
  
  # try several places beta estimates can live across sommer versions
  b <- NULL
  # direct slot (common in many versions)
  b <- tryCatch(fit$beta.hat, error = function(e) NULL)
  if (is.null(b)) b <- tryCatch(fit$beta, error = function(e) NULL)
  # summary() table fallback
  if (is.null(b)) {
    sb <- tryCatch(summary(fit), error = function(e) NULL)
    if (is.list(sb)) {
      cand <- NULL
      # common names seen in printed summaries
      if ("beta" %in% names(sb)) cand <- sb$beta
      if (is.null(cand) && "Beta" %in% names(sb)) cand <- sb$Beta
      if (is.data.frame(cand) && "Estimate" %in% names(cand)) b <- cand$Estimate
    }
  }
  # coerce to numeric vector
  if (is.null(b)) {
    # ultra-safe fallback: with ~1, use intercept-only mean of response column if present
    if (ncol(X) == 1L) {
      return(rep(mean(data[[all.vars(as.formula(formula_str))[1]]], na.rm = TRUE), nrow(X)))
    } else {
      stop("Could not locate fixed-effects estimates in 'fit' to build Xβ.")
    }
  }
  b <- as.numeric(b)
  # handle length mismatch (e.g., intercept-only)
  if (length(b) == 1L && ncol(X) >= 1L) b <- rep(b, ncol(X))
  as.numeric(drop(X %*% b))
}

spd_check <- function(df, cols) {
  cc <- stats::complete.cases(df[, cols, drop=FALSE])
  S <- tryCatch(stats::cov(df[cc, cols, drop=FALSE]), error=function(e) NULL)
  if (is.null(S) || anyNA(S)) return(list(min_eig=NA_real_, n_neg=NA_integer_))
  ev <- tryCatch(eigen(S, symmetric=TRUE, only.values=TRUE)$values, error=function(e) NULL)
  if (is.null(ev)) return(list(min_eig=NA_real_, n_neg=NA_integer_))
  list(min_eig=min(ev, na.rm=TRUE), n_neg=sum(ev < 0, na.rm=TRUE))
}

qr_drop_dependent_genes <- function(df, model_cols, mediator_cols, gene_cols) {
  cc <- stats::complete.cases(df[, model_cols, drop=FALSE])
  X  <- as.matrix(df[cc, model_cols, drop=FALSE])
  if (nrow(X) < 2 || ncol(X) < 2) return(character(0))
  Q  <- tryCatch(qr(X), error=function(e) NULL)
  if (is.null(Q)) return(character(0))
  keep_idx  <- sort(Q$pivot[seq_len(Q$rank)])
  keep_cols <- colnames(X)[keep_idx]
  dep_cols  <- setdiff(colnames(X), keep_cols)
  intersect(dep_cols, gene_cols)
}

numify <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.numeric(x))
  if (inherits(x, "Date") || inherits(x, "POSIXt")) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[nchar(x) == 0 | toupper(x) == "NA"] <- NA_character_
  x <- gsub(",", ".", x, fixed=TRUE)
  suppressWarnings(as.numeric(x))
}

## Align a named vector to the filtered cohort ids (strict, informative)
vec_from <- function(v_full, keys = cohort_ids, name = deparse(substitute(v_full))) {
  if (is.null(v_full)) stop(sprintf("%s is NULL.", name))
  if (is.null(names(v_full))) {
    stop(sprintf("%s has no names; cannot align to cohort IDs.", name))
  }
  # normalize both sides like you did elsewhere
  v_names <- norm_id(names(v_full))
  k_norm  <- norm_id(keys)
  
  ord <- match(k_norm, v_names)
  if (anyNA(ord)) {
    missing <- keys[is.na(ord)]
    stop(sprintf("%s is missing %d cohort IDs. Examples: %s",
                 name, length(missing), paste(head(missing, 10), collapse = ", ")))
  }
  setNames(as.numeric(v_full[ord]), keys)
}


## ======= preflight to fail early with clear diagnostics (NO auto-orthogonalization) =======
preflight_sem <- function(df, mediators, genes, outcome, tag = "") {
  message(sprintf("[preflight%s] Running SEM stability checks ...", ifelse(tag=="","",paste0(":",tag))))
  model_cols <- unique(c(mediators, genes, outcome))
  sds <- sapply(df[, model_cols, drop=FALSE], function(z) sd(numify(z), na.rm=TRUE))
  if (any(!is.finite(sds) | sds == 0)) {
    bad <- names(sds)[!is.finite(sds) | sds == 0]
    stop("Constant/invalid variables: ", paste(bad, collapse=", "))
  }
  if (length(mediators) > 1) {
    Cmed <- suppressWarnings(cor(df[, mediators, drop=FALSE], use="pairwise.complete.obs"))
    diag(Cmed) <- 0
    dup_idx <- which(abs(Cmed) >= 0.999999, arr.ind = TRUE)
    if (nrow(dup_idx) > 0) {
      pairs <- apply(dup_idx, 1, function(r) paste(rownames(Cmed)[r[1]], colnames(Cmed)[r[2]], sep=" ~ "))
      stop("Near-duplicate mediators detected (|r|>=0.999999): ",
           paste(unique(pairs), collapse=" ; "),
           "\nThis should be impossible by design. Check earlier steps and sample alignment.")
    }
  }
  cc <- stats::complete.cases(df[, model_cols, drop=FALSE])
  Xcc <- df[cc, model_cols, drop=FALSE]
  if (nrow(Xcc) < 2) stop("Not enough complete cases after listwise deletion.")
  al <- alias(lm(reformulate(model_cols[model_cols != outcome], outcome), data=Xcc))
  if (!is.null(al$Complete)) {
    stop("Exact linear dependencies (alias) detected among predictors:\n",
         paste(capture.output(print(al$Complete)), collapse="\n"))
  }
  qr_rank <- qr(as.matrix(stats::model.matrix(~ . - 1, data=Xcc)))$rank
  if (qr_rank < ncol(Xcc)) {
    stop(sprintf("Rank deficiency: rank=%d < p=%d", qr_rank, ncol(Xcc)))
  }
  S <- cov(Xcc, use="pairwise.complete.obs")
  ev <- eigen(S, symmetric=TRUE, only.values=TRUE)$values
  message(sprintf("[preflight%s] min eigen(pairwise cov) = %.3e",
                  ifelse(tag=="","",paste0(":",tag)), min(ev)))
  if (min(ev) < -1e-10) {
    stop("Sample covariance is not positive-definite. Resolve collinearity/duplicates first.")
  }
  invisible(TRUE)
}

## ======= Residual-covariance discovery (pre-lavaan) =======
# Get ancestors from an adjacency list (character vectors)

# Detect strong residual covariances among candidate variables.
# - skips pairs that are already linked by a directed parent->child edge in `adj`
#   (e.g., if Prev ~ Bact exists, do NOT add Bact ~~ Prev).
# - keeps pairs with |r| >= `thr` and p < `alpha`.
detect_resid_covariances <- function(Z, adj, candidates,
                                     alpha = 0.05, thr = 0.10,
                                     forbid_parent_child = TRUE) {
  candidates <- intersect(unique(candidates), colnames(Z))
  if (length(candidates) < 2) return(character(0))
  
  # Parent-child undirected keys (to block exogenous<->endogenous residual cov)
  parent_child <- character(0)
  if (forbid_parent_child && length(adj)) {
    for (lhs in names(adj)) {
      if (!length(adj[[lhs]])) next
      for (rhs in adj[[lhs]]) {
        parent_child <- c(parent_child, paste(sort(c(lhs, rhs)), collapse = "||"))
      }
    }
    parent_child <- unique(parent_child)
  }
  
  # All unordered candidate pairs
  prs <- combn(candidates, 2, simplify = FALSE)
  keep <- character(0)
  for (p in prs) {
    key <- paste(sort(p), collapse = "||")
    if (forbid_parent_child && key %in% parent_child) next
    
    x <- Z[[p[1]]]; y <- Z[[p[2]]]
    ok <- stats::complete.cases(x, y)
    if (sum(ok) < 10) next
    
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok]))
    if (is.finite(ct$estimate) && abs(ct$estimate) >= thr &&
        is.finite(ct$p.value)   && ct$p.value < alpha) {
      keep <- c(keep, sprintf("%s ~~ %s", p[1], p[2]))
    }
  }
  unique(keep)
}



## -------- CLI --------
opt_list <- list(
  make_option("--cohort", type="character", help="TXT with one sample ID per line (final cohort)"),
  make_option("--med_bac_tsv",  type="character", help="mediation_sig_iv.tsv (Bact pipeline)"),
  make_option("--med_prev_tsv", type="character", help="mediation_sig_iv.tsv (Prev pipeline)"),
  make_option("--putative_genes_xlsx", type="character", help="putative_genes_final.xlsx (optional)"),
  make_option("--selected_genes_rdata", type="character", help="genotypes_candidates.gds"),
  make_option("--resid_rda",    type="character"),
  make_option("--beta_rda",     type="character"),
  make_option("--alpha_rda",    type="character"),
  make_option("--clusters_rda", type="character"),
  make_option("--phen_rda",     type="character"),
  make_option("--grm_rda",      type="character"),
  make_option("--dam_rds",      type="character", default=NA),
  make_option("--cage_rds",     type="character", default=NA),
  make_option("--metadata",     type="character", default=NA),
  make_option("--meta_rfid_col",type="character", default="RFID"),
  make_option("--meta_dam_col", type="character", default="dam"),
  make_option("--meta_cage_col",type="character", default="cage"),
  make_option("--bac_name",    type="character", default="g__Bacteroides"),
  make_option("--prev_name",   type="character", default="g__Prevotella"),
  make_option("--guild_name", type="character", default="cluster__3"),
  make_option("--beta_trait",  type="character", default="beta__PD_PC1"),
  make_option("--alpha_trait",  type="character", default="alpha__PD_q2"),
  make_option("--glucose_col",    type="character", default="Glucose"),
  make_option("--bmi_col", type="character", default="BMI",
              help="Optional second outcome column in phen_rda (e.g., 'BMI')."),
  make_option("--fdr_level",      type="double",  default=0.10),
  make_option("--mi_threshold",   type="double",  default=10.0),
  make_option("--max_iter",       type="integer", default=15),
  make_option("--iv_gold_tier",       type="logical", default=TRUE),
  make_option("--missing_mode",   type="character", default="fiml", help="fiml|listwise"),
  make_option("--drop_corr_threshold", type="double", default=0.99,
              help="Drop one of each gene pair with |r| >= this (pairwise.complete.obs)"),
  make_option("--qr_drop",        type="logical", default=TRUE,
              help="After corr-pruning, QR-drop linearly dependent gene columns (complete cases only)"),
  make_option("--jitter_sd",      type="double", default=0.0,
              help="If still non-SPD, add tiny N(0,sd^2) jitter to gene columns (0 = off)"),
  make_option("--print_mmer_summary", type="logical", default=TRUE,
              help="Print sommer summary() per target"),
  make_option("--drop_genes", type="logical", default=TRUE,
              help="If TRUE, exclude all host-gene instruments and skip genotype reads; Bact is the only exogenous variable."),
  make_option("--components", type="character", default="", help="..."),
  make_option("--model_file", type="character", default=NA, help="..."),
  make_option("--model_text", type="character", default="", help="..."),
  make_option("--outdir", type="character", default="/users/abaud/fmorillo/paper_figures/sem_reports/")
)
opt <- parse_args(OptionParser(option_list = opt_list))
DROP_GENES <- isTRUE(opt$drop_genes)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## ================================================================
## 0) Cohort
## ================================================================
stopifnot(!is.null(opt$cohort), file.exists(opt$cohort))
cohort_raw <- read_lines(opt$cohort) |> (\(z) z[nzchar(z)])()
cohort_ids_all <- uniq(norm_id(cohort_raw))
msg("Cohort loaded: ", length(cohort_ids_all), " unique IDs (after normalization)")

## ================================================================
## 1) Load residual layers / GRM; subset to cohort; DROP NA Glucose EARLY
## ================================================================
msg("Loading residual matrices, phenotypes, and covariance matrices...")
stopifnot(file.exists(opt$resid_rda), file.exists(opt$beta_rda),
          file.exists(opt$clusters_rda), file.exists(opt$phen_rda),
          file.exists(opt$grm_rda))
if (!is.na(opt$alpha_rda)) stopifnot(file.exists(opt$alpha_rda))

load(opt$resid_rda)
load(opt$beta_rda)
if (!is.na(opt$alpha_rda)) load(opt$alpha_rda)
load(opt$clusters_rda); load(opt$phen_rda); load(opt$grm_rda)

use_DamCage_mats <- !is.na(opt$dam_rds) && file.exists(opt$dam_rds) &&
  !is.na(opt$cage_rds) && file.exists(opt$cage_rds)
Dam <- Cage <- NULL
if (use_DamCage_mats) { Dam <- as.matrix(readRDS(opt$dam_rds)); Cage <- as.matrix(readRDS(opt$cage_rds)) }

## genus layer
genus_layer_idx <- length(residuals_qned_counts_objs) - 1L
R_genus <- t(residuals_qned_counts_objs[[genus_layer_idx]])

## DROP an optional female-specific Bacteroides row if present (avoid dimname mismatch)
if ("g__Bacteroides_F" %in% rownames(R_genus)) {
  R_genus <- R_genus[rownames(R_genus) != "g__Bacteroides_F", , drop = FALSE]
}

## THEN normalize row IDs
rownames(R_genus) <- norm_id(rownames(R_genus))

bac_col  <- pick_col(R_genus, opt$bac_name)
prev_col <- pick_col(R_genus, opt$prev_name)
stopifnot(!is.na(bac_col), !is.na(prev_col))
Bac_full   <- setNames(as.numeric(R_genus[, bac_col ]), rownames(R_genus))
Prev_full  <- setNames(as.numeric(R_genus[, prev_col]),  rownames(R_genus))

## clusters (Guild parsed; optional in SEM)
Cmat <- residuals_qned_counts_clusters_objs[[1]]
rExact <- which(rownames(Cmat) == opt$guild_name)
rPick  <- if (length(rExact)) rExact[1] else which(startsWith(rownames(Cmat), opt$guild_name))[1]
Guild_full <- setNames(as.numeric(Cmat[rPick, ]), norm_id(colnames(Cmat)))

## beta
BETA <- as.matrix(residuals_qned_counts_beta_objs[[1]])
beta_row <- if (opt$beta_trait %in% c("beta__PD_PC1","Axis.1")) 5L else if (opt$beta_trait %in% c("beta__PD_PC2","Axis.2")) 6L else 5L
Beta_full <- setNames(as.numeric(BETA[beta_row, ]), norm_id(colnames(BETA)))

## alpha (optional)
Alpha_full <- NULL
if (!is.na(opt$alpha_rda)) {
  ALPHA <- tryCatch(as.matrix(residuals_qned_counts_alpha_objs[[1]]), error=function(e) NULL)
  if (!is.null(ALPHA)) {
    # prefer a rowname match; otherwise fall back to the first row
    alpha_row <- if (opt$alpha_trait %in% rownames(ALPHA)) which(rownames(ALPHA) == opt$alpha_trait)[1] else 1L
    Alpha_full <- setNames(as.numeric(ALPHA[alpha_row, ]), norm_id(colnames(ALPHA)))
  } else {
    msg("WARNING: 'residuals_qned_counts_alpha_objs' not found in --alpha_rda; skipping Alpha.")
  }
}


## phenotype
df$samples <- norm_id(df$samples)
Gluc_full  <- setNames(df[[opt$glucose_col]], df$samples)
## optional second outcome
BMI_full <- NULL
if (!is.na(opt$bmi_col)) {
  BMI_full <- setNames(df[[opt$bmi_col]], df$samples)
}

## ---- EARLY non-NA mask applied to EVERYTHING ----
if (!is.null(BMI_full)) {
  msg("Applying early non-NA filter for Glucose & BMI to lock the sample set...")
  has_both <- cohort_ids_all %in% names(Gluc_full)  & !is.na(Gluc_full [cohort_ids_all]) &
    cohort_ids_all %in% names(BMI_full) & !is.na(BMI_full[cohort_ids_all])
  drop_n <- sum(!has_both)
  if (drop_n > 0) msg("Dropping ", drop_n, " cohort IDs due to missing Glucose or BMI.")
  cohort_ids <- cohort_ids_all[has_both]
} else {
  msg("Applying early Glucose non-NA filter to lock the sample set for all steps...")
  has_gluc <- cohort_ids_all %in% names(Gluc_full) & !is.na(Gluc_full[cohort_ids_all])
  drop_n   <- sum(!has_gluc)
  if (drop_n > 0) msg("Dropping ", drop_n, " cohort IDs due to missing Glucose.")
  cohort_ids <- cohort_ids_all[has_gluc]
}

## Subset kernels to filtered cohort now
stop_if_missing("GRM", rownames(grm), cohort_ids); grm <- grm[cohort_ids, cohort_ids]
if (use_DamCage_mats) {
  stop_if_missing("Dam", rownames(Dam), cohort_ids); stop_if_missing("Cage", rownames(Cage), cohort_ids)
  Dam <- Dam[cohort_ids, cohort_ids]; Cage <- Cage[cohort_ids, cohort_ids]
}

## Pull vectors for filtered cohort
Bact     <- vec_from(Bac_full);  Prev     <- vec_from(Prev_full)
Guild    <- vec_from(Guild_full); Beta <- vec_from(Beta_full)
Alpha    <- if (!is.null(Alpha_full)) vec_from(Alpha_full) else NULL
Glucose  <- vec_from(Gluc_full)
BMI      <- if (!is.null(BMI_full)) vec_from(BMI_full) else NULL

## ================================================================
## 2) Instruments (all tiers), split sets
## ================================================================
msg("Reading mediation outputs and filtering to instruments")

if (isTRUE(DROP_GENES)) {
  med_bac  <- data.table()
  med_prev <- data.table()
  msg("DROP_GENES=TRUE: skipping instrument selection (no gene exogenous variables).")
} else {
  if (opt$iv_gold_tier == TRUE){
    med_bac  <- data.table::fread(opt$med_bac_tsv)[tier=="gold"]
    med_prev <- data.table::fread(opt$med_prev_tsv)[tier=="gold"]
  } else{
    med_bac  <- data.table::fread(opt$med_bac_tsv)
    med_prev <- data.table::fread(opt$med_prev_tsv)
  }

  med_bac  <- dplyr::distinct(med_bac,  snp, mediator, .keep_all=TRUE)
  med_prev <- dplyr::distinct(med_prev, snp, mediator, .keep_all=TRUE)

  mk_lab <- function(dt){
    if (!"snp" %in% names(dt)) dt$snp <- NA_character_
    dt$chr <- as.character(dt$chr); dt$pos <- as.integer(dt$pos)
    dt$snp_label <- ifelse(!is.na(dt$snp) & grepl("^chr", as.character(dt$snp)),
                           as.character(dt$snp),
                           paste0("chr", dt$chr, ":", dt$pos))
    
    dt
  }
  med_bac  <- mk_lab(med_bac)
  med_prev <- mk_lab(med_prev)

  if (!is.null(opt$putative_genes_xlsx) && file.exists(opt$putative_genes_xlsx)) {
    pg <- readxl::read_excel(opt$putative_genes_xlsx)
    med_bac$Bioprocess  <- pg$Bioprocess[match(med_bac$gene_name,  pg$Gene)]
    med_prev$Bioprocess <- pg$Bioprocess[match(med_prev$gene_name, pg$Gene)]
  }
}

if (!isTRUE(DROP_GENES)) {
  fwrite(rbindlist(list(med_bac,med_prev), use.names=TRUE, fill=TRUE),
         file=file.path(opt$outdir,"selected_instruments.tsv"), sep="\t")
  msg("Instruments selected (rows read): ", nrow(med_bac)+nrow(med_prev))
} else {
  msg("Instruments selected (rows read): 0")
}

## ================================================================
## 3) Genotypes from GDS (subset to the EARLY-filtered cohort)
## ================================================================
read_genotypes_from_gds <- function(gds_path, instr_dt, keep_ids) {
  msg("Opening GDS: ", gds_path)
  g <- snpgdsOpen(gds_path); on.exit(snpgdsClose(g), add=TRUE)
  snp.id  <- read.gdsn(index.gdsn(g,"snp.id"))
  snp.lab <- tryCatch(read.gdsn(index.gdsn(g,"snp.label")), error=function(e) NULL)
  if (is.null(snp.lab)) {
    snp.chr <- read.gdsn(index.gdsn(g,"snp.chromosome"))
    snp.pos <- read.gdsn(index.gdsn(g,"snp.position"))
    snp.lab <- paste0("chr", snp.chr, ":", snp.pos)
  }
  ann <- data.table(snp.id=snp.id, snp_label=as.character(snp.lab))
  M <- merge(data.table::copy(instr_dt), ann, by="snp_label", all.x=TRUE)
  if (anyNA(M$snp.id)) {
    miss <- M[is.na(snp.id)]
    if (nrow(miss)) {
      msg("WARNING: ", nrow(miss), " instrument SNPs not in GDS; skipping. Examples: ",
          paste(head(unique(miss$snp_label),3), collapse=", "))
    }
    M <- M[!is.na(snp.id)]
  }
  if (nrow(M) == 0L) {
    msg("WARNING: 0 instrument SNPs from this set were found in the GDS — returning empty dosage/meta.")
    return(list(dosage = data.table(id = keep_ids),
                meta   = data.table()))
  }
  g_sample_id <- read.gdsn(index.gdsn(g,"sample.id"))
  g_sample_lb <- tryCatch(read.gdsn(index.gdsn(g,"sample.label")), error=function(e) NULL)
  if (is.null(g_sample_lb)) g_sample_lb <- as.character(g_sample_id) else g_sample_lb <- as.character(g_sample_lb)
  keep_norm <- norm_id(keep_ids); g_norm <- norm_id(g_sample_lb)
  map_idx <- match(keep_norm, g_norm)
  if (any(is.na(map_idx))) {
    bad <- keep_ids[is.na(map_idx)]
    stop(sprintf("Genotypes GDS is missing %d cohort IDs. Examples: %s",
                 sum(is.na(map_idx)), paste(head(bad,10), collapse=", ")))
  }
  sel_sid <- g_sample_id[map_idx]
  G <- snpgdsGetGeno(g, sample.id=sel_sid, snp.id=M$snp.id, with.id=TRUE)
  X <- t(G$genotype); if (any(X==3L, na.rm=TRUE)) X[X==3L] <- NA_integer_
  colnames(X) <- make.names(M$gene_name, unique=TRUE)
  list(dosage = data.table(id=keep_ids, X),
       meta   = unique(M[, intersect(c("gene_name","snp","chr","pos","mediator","Bioprocess","snp_label","snp.id"), names(M)), with=FALSE]))
}

if (!isTRUE(DROP_GENES)) {
  stopifnot(file.exists(opt$selected_genes_rdata))
  g_bac  <- read_genotypes_from_gds(opt$selected_genes_rdata, med_bac,  cohort_ids)
  g_prev <- read_genotypes_from_gds(opt$selected_genes_rdata, med_prev, cohort_ids)
  
  GENO_ALL_dt <- as.data.table(Reduce(function(a,b) merge(a,b,by="id",all=TRUE), list(g_bac$dosage, g_prev$dosage)))
  GENO_ALL_dt <- GENO_ALL_dt[match(cohort_ids, GENO_ALL_dt$id), ]
  stopifnot(identical(GENO_ALL_dt$id, cohort_ids))
  fwrite(rbindlist(list(g_bac$meta,g_prev$meta), use.names=TRUE, fill=TRUE),
         file=file.path(opt$outdir,"selected_instruments_with_ids.tsv"), sep="\t")
  msg("Genotype matrix ready: n=", nrow(GENO_ALL_dt), " samples | p=", ncol(GENO_ALL_dt)-1, " instrument columns")
} else {
  # No gene instruments: keep an id-only design matrix
  GENO_ALL_dt <- data.table(id = cohort_ids)
  g_bac  <- list(meta = data.table())
  g_prev <- list(meta = data.table())
  msg("DROP_GENES=TRUE: skipping genotype reads; proceeding without host-gene instruments.")
}


if (!isTRUE(DROP_GENES)) {
  gene_bac_names  <- unique(med_bac$gene_name)
  gene_prev_names <- unique(med_prev$gene_name)
} else {
  gene_bac_names  <- character(0)
  gene_prev_names <- character(0)
}

genes_bac_only_names  <- setdiff(gene_bac_names,  gene_prev_names)
genes_prev_only_names <- setdiff(gene_prev_names, gene_bac_names)
genes_overlap_names   <- intersect(gene_bac_names, gene_prev_names)

all_gene_cols <- setdiff(colnames(GENO_ALL_dt), "id")
cols_for <- function(gset) {
  if (!length(gset)) return(character(0))
  pref <- unique(make.names(gset, unique=FALSE))
  mask <- Reduce(`|`, lapply(pref, function(p) startsWith(all_gene_cols, p)), init=rep(FALSE,length(all_gene_cols)))
  all_gene_cols[mask]
}

genes_bact_cols0 <- if (isTRUE(DROP_GENES)) character(0) else cols_for(genes_bac_only_names)
genes_prev_cols0 <- if (isTRUE(DROP_GENES)) character(0) else cols_for(genes_prev_only_names)
genes_both_cols0 <- if (isTRUE(DROP_GENES)) character(0) else cols_for(genes_overlap_names)
genes_union0     <- unique(c(genes_bact_cols0, genes_prev_cols0, genes_both_cols0))

msg(sprintf("Gene columns split: Bact-only=%d, Prev-only=%d, Overlap=%d; union=%d",
            length(genes_bact_cols0), length(genes_prev_cols0), length(genes_both_cols0), length(genes_union0)))


## ================================================================
## 4) Build my_data + residualize via sommer
## ================================================================
if (is.na(opt$metadata) || !file.exists(opt$metadata)) {
  stop("metadata file is required to build 'dam' and 'cage' factors. Provide --metadata with columns '",
       opt$meta_rfid_col, "', '", opt$meta_dam_col, "', '", opt$meta_cage_col, "'.")
}
msg("Reading metadata for dam/cage (as factors): ", opt$metadata)
Meta <- tryCatch(suppressWarnings(fread(opt$metadata)), error=function(e) NULL)
if (is.null(Meta)) Meta <- readr::read_delim(opt$metadata, delim="\t", guess_max=1e5, show_col_types=FALSE)
stopifnot(opt$meta_rfid_col %in% names(Meta),
          opt$meta_dam_col  %in% names(Meta),
          opt$meta_cage_col %in% names(Meta))
Meta[[opt$meta_rfid_col]] <- norm_id(Meta[[opt$meta_rfid_col]])

gene_cols_all <- setdiff(colnames(GENO_ALL_dt), "id")

# Start with the phenos/mediators only
my_data <- data.frame(
  sample   = cohort_ids,
  Bact     = Bact,
  Prev     = Prev,
  Guild    = Guild,
  Beta     = Beta,
  Glucose  = Glucose,
  stringsAsFactors = FALSE
)
if (!is.null(Alpha)) my_data$Alpha <- Alpha
if (!is.null(BMI))   my_data$BMI   <- BMI



# If we actually have gene columns, append them (row-aligned)
if (length(gene_cols_all) > 0) {
  geno_block <- as.data.frame(
    GENO_ALL_dt[match(cohort_ids, GENO_ALL_dt$id), ..gene_cols_all]
  )
  stopifnot(nrow(geno_block) == nrow(my_data))
  my_data <- cbind(my_data, geno_block)
}

meta_idx <- match(my_data$sample, Meta[[opt$meta_rfid_col]])
if (any(is.na(meta_idx))) {
  stop("metadata is missing cohort IDs. Examples: ",
       paste(head(my_data$sample[is.na(meta_idx)], 10), collapse=", "))
}
my_data$dam  <- factor(Meta[[opt$meta_dam_col ]][meta_idx])
my_data$cage <- factor(Meta[[opt$meta_cage_col]][meta_idx])

rownames(my_data) <- my_data$sample
my_data <- my_data[rownames(grm), , drop=FALSE]
stop_if_missing("GRM", rownames(grm), my_data$sample)

if (length(gene_cols_all)) {
  my_data[, gene_cols_all] <- lapply(my_data[, gene_cols_all, drop=FALSE], numify)
}
my_data <- na.omit(my_data)

msg("Fitting sommer models with random effects: GRM(sample) + dam + cage (no fixed effects)")
present_nodes <- intersect(c("Bact","Prev","Guild","Alpha","Beta","Glucose","BMI"), colnames(my_data))
mediators <- intersect(c("Bact","Prev","Guild","Alpha","Beta"), present_nodes)
OUTCOMES  <- intersect(c("Glucose","BMI"), present_nodes)
targets   <- c(mediators, OUTCOMES)

fit_list <- list()
random_summary_list <- list()

for (lh in targets) {
  msg("sommer fit: ", lh)
  lhs_var <- as.formula(paste0(lh, " ~ 1"))
  fit <- mmer(
    fixed  = lhs_var,
    random = ~
      vsr(sample, Gu = grm) +
      vsr(dam) +
      vsr(cage),
    rcov   = ~ units,
    data   = my_data,
    verbose = FALSE,
    dateWarning = FALSE, date.warning = FALSE
  )
  if (isTRUE(opt$print_mmer_summary)) {
    print(summary(fit))
  }
  fit_list[[lh]] <- fit
  random_summary_list[[lh]] <- summary(fit)
}

my_data2 <- my_data
for (lh in targets) {
  # residuals e = y - (Xβ + Zu)
  e  <- resid_vec_from_mmer(
    fit     = fit_list[[lh]],
    y       = my_data[[lh]],
    row_key = my_data$sample,
    trait   = lh
  )
  # fixed part Xβ (robust to beta.hat shape)
  Xb <- fixed_part_vec_from_mmer(
    fit         = fit_list[[lh]],
    formula_str = paste0(lh, " ~ 1"),
    data        = my_data
  )
  # remove random effects only: y_noRE = e + Xβ = y - Zu
  y_noRE <- e + Xb
  my_data2[[lh]] <- as.numeric(y_noRE)
}

my_data2$id <- my_data2$sample
DF <- my_data

msg("----- HEADS: before/after sommer (your 5 targets) -----")
print(utils::head(my_data[,  c("sample", mediators, OUTCOMES)], 6))
print(utils::head(my_data2[, c("sample", mediators, OUTCOMES)], 6))

Cmed <- suppressWarnings(cor(my_data2[, c(mediators, OUTCOMES)], use="pairwise.complete.obs"))
if (any(abs(Cmed[lower.tri(Cmed)]) > 0.95, na.rm = TRUE)) {
  msg("WARNING: Residualized targets show |r| > 0.95 among themselves. See residuals_debug_cor.tsv")
  data.table::fwrite(as.data.table(Cmed, keep.rownames = "var"),
                     file = file.path(opt$outdir, "residuals_debug_cor.tsv"), sep = "\t")
}

msg("Running lm() + VIF diagnostics (on residualized values)...")
lm_dir <- file.path(opt$outdir, "lm_vif"); dir.create(lm_dir, showWarnings = FALSE, recursive = TRUE)
vif_safe <- function(fm) {
  k <- length(coef(fm)) - 1L
  if (k <= 1) return(setNames(1, names(coef(fm))[-1]))
  tryCatch(car::vif(fm), error=function(e) setNames(rep(NA_real_, k), names(coef(fm))[-1]))
}
run_lm <- function(formula_str, name) {
  dat <- my_data2
  fm <- lm(as.formula(formula_str), data = dat)
  s  <- capture.output(print(summary(fm)))
  V  <- vif_safe(fm)
  writeLines(c(paste0("MODEL: ", formula_str), s, "", "VIF:", capture.output(print(V)) ),
             con = file.path(lm_dir, paste0("lm_", name, ".txt")))
  msg("lm() done: ", name, "  |  ", formula_str)
  print(summary(fm)); print(V)
}
run_lm("Guild ~ Prev + Bact",                     "guild_on_prev_bact")
run_lm("Beta ~ Guild + Prev + Bact",              "beta_on_guild_prev_bact")
#run_lm("Beta ~ Prev + Bact",                      "beta_on_prev_bact")
run_lm("Alpha ~ Guild + Prev + Bact",             "alpha_on_guild_prev_bact")
#run_lm("Alpha ~ Prev + Bact",                     "alpha_on_prev_bact")
run_lm("Glucose ~ Beta + Guild + Prev + Bact",    "glucose_on_beta_guild_prev_bact")
run_lm("Glucose ~ Alpha + Guild + Prev + Bact",    "glucose_on_alpha_guild_prev_bact")
run_lm("Glucose ~ Alpha + Beta + Guild + Prev + Bact",    "glucose_on_alpha_beta_guild_prev_bact")
#run_lm("Glucose ~ Beta + Prev + Bact",            "glucose_on_beta_prev_bact")
run_lm("BMI ~ Beta + Guild + Prev + Bact",        "BMI_on_beta_guild_prev_bact")
run_lm("BMI ~ Alpha + Guild + Prev + Bact",        "BMI_on_alpha_guild_prev_bact")
run_lm("BMI ~ Alpha + Beta + Guild + Prev + Bact",        "BMI_on_alpha_beta_guild_prev_bact")
#run_lm("BMI ~ Beta + Prev + Bact",                "BMI_on_beta_prev_bact")
msg("lm() + VIF outputs: ", lm_dir)

#msg("Analysing random effects (variance components)")
#mediators_all <- names(random_summary_list)
#vc_tab <- do.call(rbind, lapply(mediators_all, function(m) {
#  vc <- random_summary_list[[m]]$varcomp$VarComp
#  # define the names you expect (adjust if your labels differ)
#  wanted <- c("sample", "dam", "cage", "units")
#  # get the actual rownames
#  rn <- rownames(vc)
#  
#  # helper to extract safely
#  get_component <- function(name) {
#    idx <- grep(paste0("^", name, "$"), rn)
#    if (length(idx)) return(as.numeric(vc[idx, 1]))
#    warning(sprintf("Component '%s' not found in VarComp rows: %s", name, paste(rn, collapse=", ")))
#    return(NA_real_)
#  }
#  
#  samp  <- get_component("sample")
#  mat   <- get_component("dam")
#  coh   <- get_component("cage")
#  resid <- get_component("units")
#  totalR <- sum(c(samp, mat, coh), na.rm = TRUE)
#  
#  data.frame(
#    Mediator = m, GRM = samp, Maternal = mat, Cohousing= coh, Residual = resid, TotalRE = totalR,
#    stringsAsFactors = FALSE
#  )
#}))
#row.names(vc_tab) <- NULL
#vc_tab$Prop_GRM      <- vc_tab$GRM      / vc_tab$TotalRE
#vc_tab$Prop_Maternal <- vc_tab$Maternal / vc_tab$TotalRE
#vc_tab$Prop_Cohousing<- vc_tab$Cohousing/ vc_tab$TotalRE
#classify <- function(p) { if (is.na(p)) "Negligible" else if (p >= 0.5) "Major" else if (p >= 0.1) "Minor" else "Negligible" }
#vc_tab$Impact_GRM      <- sapply(vc_tab$Prop_GRM,      classify)
#vc_tab$Impact_Maternal <- sapply(vc_tab$Prop_Maternal, classify)
#vc_tab$Impact_Cohousing<- sapply(vc_tab$Prop_Cohousing,classify)
#data.table::fwrite(vc_tab, file=file.path(opt$outdir, "variance_components_by_target.tsv"), sep="\t")
#impact_counts <- data.frame(
#  Effect     = c("GRM","Maternal","Cohousing"),
#  Major      = c(sum(vc_tab$Impact_GRM      == "Major"),
#                 sum(vc_tab$Impact_Maternal == "Major"),
#                 sum(vc_tab$Impact_Cohousing== "Major")),
#  Minor      = c(sum(vc_tab$Impact_GRM      == "Minor"),
#                 sum(vc_tab$Impact_Maternal == "Minor"),
#                 sum(vc_tab$Impact_Cohousing== "Minor")),
#  Negligible = c(sum(vc_tab$Impact_GRM      == "Negligible"),
#                 sum(vc_tab$Impact_Maternal == "Negligible"),
#                 sum(vc_tab$Impact_Cohousing== "Negligible")),
#  stringsAsFactors = FALSE
#)
#impact_counts <- impact_counts[order(-impact_counts$Major, -impact_counts$Minor), ]
#data.table::fwrite(impact_counts, file=file.path(opt$outdir, "random_effects_importance.tsv"), sep="\t")
#most_idx <- which.max(vc_tab$TotalRE); most_med <- vc_tab$Mediator[most_idx]
#sink(file.path(opt$outdir, "random_effects_summary.txt"))
#cat("\nVariance components by mediator:\n")
#print(vc_tab[, c("Mediator","GRM","Maternal","Cohousing","Residual",
#                 "Impact_GRM","Impact_Maternal","Impact_Cohousing")], row.names = FALSE)
#cat("\n\nSummary of impacts:\n")
#apply(vc_tab, 1, function(r) {
#  cat(sprintf(" • %s: GRM[%s], Maternal[%s], Cohousing[%s]\n",
#              r["Mediator"], r["Impact_GRM"], r["Impact_Maternal"], r["Impact_Cohousing"]))
#})
#cat("\n\nThe mediator most influenced by the combined random effects is:",
#    most_med,
#    sprintf("(%.1f%% of its variance is random effects).",
#            vc_tab$TotalRE[most_idx] /
#              (vc_tab$TotalRE[most_idx] + vc_tab$Residual[most_idx]) * 100),
#    "\n\n")
#cat("\nRandom‐effect importance ranking:\n")
#print(impact_counts, row.names = FALSE)
#top <- impact_counts$Effect[1]; maj <- impact_counts$Major[1]; min <- impact_counts$Minor[1]
#cat(sprintf(
#  "\nThe most important random effect is **%s**, with %d Major and %d Minor impacts across your mediators.\n",
#  top, maj, min))
#cat("**Recommendation:** adjust for this random effect (use residuals from these mixed models) ",
#    "before SEM/BSEM to avoid confounding.\n")
#sink()

scale_cols <- c(gene_cols_all, mediators, OUTCOMES)
scale_cols <- intersect(colnames(my_data2), scale_cols)
data_matrix <- scale(as.data.frame(my_data2[, scale_cols, drop=FALSE]), center=TRUE, scale=TRUE)
my_data2_std <- cbind(my_data2[, c("sample","id"), drop=FALSE],
                      data_matrix,
                      my_data2[, c("dam","cage"), drop=FALSE])
data.table::fwrite(as.data.table(my_data2_std),
                   file=file.path(opt$outdir, "my_data2_standardized.tsv"), sep="\t")

DF_res <- my_data2

## ================================================================
## 6) Pre-SEM gene collapse (constants, duplicates, high corr, QR)
## ================================================================
collapse_log <- list()
mediator_cols_all <- intersect(
  c("Bact","Prev","Guild","Alpha","Beta","Glucose", if (!is.null(BMI)) "BMI"),
  colnames(DF_res)
)
nongene_cols <- c("id","sample","dam","cage", mediator_cols_all, "DamF","CageF","idF")

all_gene_cols <- setdiff(colnames(DF_res), nongene_cols)
gene_cols <- all_gene_cols

sdv <- vapply(DF_res[, gene_cols, drop=FALSE], function(x) sd(numify(x), na.rm=TRUE), numeric(1))
const_drop <- names(sdv)[sdv < 1e-12]
if (length(const_drop)) {
  collapse_log$constant <- const_drop
  DF_res <- DF_res[, setdiff(colnames(DF_res), const_drop), drop=FALSE]
  gene_cols <- setdiff(gene_cols, const_drop)
}
dup_drop <- character(0)
if (length(gene_cols) > 1) {
  key_of <- function(v) paste0(ifelse(is.na(v), "NA", sprintf("%.6f", numify(v))), collapse="|")
  keys <- vapply(gene_cols, function(nm) key_of(DF_res[[nm]]), character(1))
  dup_mask <- duplicated(keys)
  if (any(dup_mask)) {
    dup_drop <- gene_cols[dup_mask]
    collapse_log$duplicates <- dup_drop
    DF_res <- DF_res[, setdiff(colnames(DF_res), dup_drop), drop=FALSE]
    gene_cols <- setdiff(gene_cols, dup_drop)
  }
}
hc_drop <- character(0)
if (length(gene_cols) > 1) {
  C <- suppressWarnings(stats::cor(as.data.frame(lapply(DF_res[, gene_cols, drop=FALSE], numify)),
                                   use="pairwise.complete.obs"))
  C[is.na(C)] <- 0
  to_drop <- rep(FALSE, length(gene_cols)); names(to_drop) <- gene_cols
  for (i in seq_along(gene_cols)) {
    if (to_drop[i]) next
    nm_i <- gene_cols[i]
    idx <- which(abs(C[nm_i, ]) >= opt$drop_corr_threshold)
    idx <- idx[idx > i]
    if (length(idx)) to_drop[idx] <- TRUE
  }
  if (any(to_drop)) {
    hc_drop <- names(to_drop)[to_drop]
    collapse_log$highcorr <- hc_drop
    DF_res <- DF_res[, setdiff(colnames(DF_res), hc_drop), drop=FALSE]
    gene_cols <- setdiff(gene_cols, hc_drop)
  }
}
qr_dropped <- character(0)
mediator_cols_all <- intersect(
  c("Bact","Prev","Guild","Alpha","Beta","Glucose", if (!is.null(BMI)) "BMI"),
  colnames(DF_res)
)

all_gene_cols <- setdiff(
  colnames(DF_res),
  c("id","sample","dam","cage", mediator_cols_all, "DamF","CageF","idF")
)
SEM_MEDIATORS     <- c("Bact","Prev","Beta")   # default; may be updated below
SEM_TARGETS       <- c(SEM_MEDIATORS, "Glucose",if (!is.null(BMI)) "BMI")
if (isTRUE(opt$qr_drop) && length(gene_cols) > 0) {
  model_cols_now <- c(gene_cols, mediator_cols_all)
  DF_num <- DF_res
  DF_num[, model_cols_now] <- lapply(DF_num[, model_cols_now, drop=FALSE], numify)
  qr_dropped <- qr_drop_dependent_genes(DF_num, model_cols_now, mediator_cols_all, gene_cols)
  if (length(qr_dropped)) {
    collapse_log$qr_drop <- qr_dropped
    DF_res <- DF_res[, setdiff(colnames(DF_res), qr_dropped), drop=FALSE]
    gene_cols <- setdiff(gene_cols, qr_dropped)
  }
}
all_gene_cols <- setdiff(colnames(DF_res), c("id", mediator_cols_all, "DamF","CageF","idF"))
genes_bact_cols <- intersect(genes_bact_cols0, colnames(DF_res))
genes_prev_cols <- intersect(genes_prev_cols0, colnames(DF_res))
genes_both_cols <- intersect(genes_both_cols0, colnames(DF_res))
genes_union     <- unique(c(genes_bact_cols, genes_prev_cols, genes_both_cols))

model_cols_all <- unique(c(genes_union, mediator_cols_all))
Z <- DF_res
Z[, model_cols_all] <- lapply(Z[, model_cols_all, drop=FALSE], numify)
Z[, model_cols_all] <- scale(as.data.frame(Z[, model_cols_all, drop=FALSE]), center=TRUE, scale=TRUE)

save(DF, DF_res, Z, collapse_log, file=file.path(opt$outdir, "harmonized_data.RData"))
msg(sprintf("Final modeling matrix: n=%d | genes_kept=%d", nrow(Z), length(genes_union)))
if (length(collapse_log)) {
  if (!is.null(collapse_log$constant))   msg("Collapsed (constant): ", paste(collapse_log$constant, collapse=", "))
  if (!is.null(collapse_log$duplicates)) msg("Collapsed (duplicates): ", paste(collapse_log$duplicates, collapse=", "))
  if (!is.null(collapse_log$highcorr))   msg("Collapsed (|r|>=", opt$drop_corr_threshold, "): ", paste(collapse_log$highcorr, collapse=", "))
  if (!is.null(collapse_log$qr_drop))    msg("Collapsed (QR dependent): ", paste(collapse_log$qr_drop, collapse=", "))
}

miss_mode <- match.arg(opt$missing_mode, c("fiml","listwise"))
maybe_jitter <- function(df, cols) {
  if (opt$jitter_sd > 0) {
    for (nm in cols) df[[nm]] <- df[[nm]] + rnorm(nrow(df), 0, opt$jitter_sd)
  }
  df
}

## ======= SEM setup (component-driven; model comes from --model_file / --model_text) =======

# 1) Parse components (what to include)
COMP_ALL <- c("Bact","Prev","Guild","Beta","Alpha","Glucose","BMI")
COMP <- unique(trimws(unlist(strsplit(opt$components, "[, ]+"))))
COMP <- COMP[nzchar(COMP)]
if (!length(COMP)) stop("--components is empty; choose from: ", paste(COMP_ALL, collapse=", "))
COMP <- intersect(COMP, COMP_ALL)
msg("Components enabled: ", paste(COMP, collapse=", "))

# 2) Decide outcomes/mediators from selected components
SEM_OUTCOMES  <- intersect(c("Glucose","BMI"), COMP)
SEM_MEDIATORS <- setdiff(intersect(c("Bact","Prev","Guild","Beta","Alpha"), COMP), SEM_OUTCOMES)
SEM_TARGETS   <- c(SEM_MEDIATORS, SEM_OUTCOMES)

# 3) Get the model skeleton text
raw_model <- NA_character_
if (!is.na(opt$model_file) && file.exists(opt$model_file)) {
  raw_model <- readr::read_file(opt$model_file)
} else if (nzchar(opt$model_text)) {
  raw_model <- opt$model_text
} else {
  stop("Provide --model_file or --model_text with the SEM skeleton.")
}
cat("\nRaw model:\n")
cat(raw_model,"\n")
# NOTE: lavaan uses `~` for regressions and `~~` for (co)variances; we’ll parse & keep only vars in --components.

# 4) Parse model text -> adjacency (regressions), fixed covariances, and a topological order
parse_model_to_state <- function(txt, keep_vars) {
  lines <- unlist(strsplit(gsub("\r\n?", "\n", txt), "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  
  covs <- grep("~~", lines, value = TRUE, fixed = TRUE)
  regs <- setdiff(lines, covs)
  
  # keep only covariances with both sides present
  kept_cov <- character(0)
  for (li in covs) {
    parts <- trimws(strsplit(li, "~~", fixed = TRUE)[[1]])
    if (length(parts) == 2 && all(parts %in% keep_vars)) {
      kept_cov <- c(kept_cov, sprintf("%s ~~ %s", parts[1], parts[2]))
    }
  }
  
  # regressions: build adjacency list
  adj <- list()
  for (li in regs) {
    if (!grepl("~", li, fixed = TRUE)) next
    bits <- strsplit(li, "~", fixed = TRUE)[[1]]
    lhs <- trimws(bits[1]); rhs <- trimws(bits[2])
    if (!(lhs %in% keep_vars)) next
    toks <- trimws(unlist(strsplit(gsub(";", "", rhs), "\\+")))
    toks <- toks[nzchar(toks) & toks %in% keep_vars]
    adj[[lhs]] <- unique(c(adj[[lhs]], toks))
  }
  
  # topological order (fallback to canonical order if needed)
  nodes <- unique(c(names(adj), unlist(adj, use.names = FALSE)))
  el <- do.call(rbind, lapply(names(adj), function(to) if (length(adj[[to]])) cbind(adj[[to]], to) else NULL))
  topo <- character(0)
  if (!is.null(el) && nrow(el) > 0) {
    g <- igraph::graph_from_edgelist(el, directed = TRUE)
    topo <- tryCatch(names(igraph::topo_sort(g, mode="out")), error=function(e) character(0))
  }
  if (!length(topo)) {
    topo <- intersect(c("Bact","Prev","Guild","Beta","Alpha","Glucose","BMI"), nodes)
  }
  
  list(adj = adj, cov_pairs = unique(kept_cov), order = topo)
}

st <- parse_model_to_state(raw_model, keep_vars = SEM_TARGETS)
COV_FIXED <- st$cov_pairs           # skeleton-defined covariances (e.g. Bact ~~ Prev)
adj_tt    <- st$adj

COV_PAIRS <- unique(c(
  COV_FIXED,
  detect_resid_covariances(
    Z, adj_tt,
    candidates = SEM_TARGETS, alpha = 0.05, thr = 0.10
  )
))
targets_order <- unique(c(st$order, setdiff(SEM_TARGETS, st$order)))


# 5) Gene instruments mapped to exposures actually present
genes_for <- setNames(rep(list(character(0)), length(SEM_TARGETS)), SEM_TARGETS)
if ("Bact" %in% SEM_TARGETS) genes_for$Bact <- genes_bact_cols
if ("Prev" %in% SEM_TARGETS) genes_for$Prev <- genes_prev_cols

compute_exogenous_ok <- function(adj, genes_for, nodes) {
  nodes <- unique(nodes)
  lhs_nodes <- intersect(names(adj), nodes)
  rhs_nodes <- intersect(unique(unlist(adj, use.names = FALSE)), nodes)
  sort(setdiff(union(lhs_nodes, rhs_nodes), lhs_nodes))
}

EXOGENOUS_OK <- compute_exogenous_ok(adj_tt, genes_for, SEM_TARGETS)
msg("Exogenous (by model skeleton): ", paste(EXOGENOUS_OK, collapse = ", "))
dropped_gene_rhs <- setNames(lapply(SEM_TARGETS, function(x) character(0)), SEM_TARGETS)

# 6) SPD check over current SEM set, switch missing mode if needed
sem_model_cols <- unique(c(genes_union, SEM_TARGETS))
chk0 <- spd_check(Z, sem_model_cols)
if (!is.na(chk0$min_eig) && chk0$min_eig < 1e-8 && miss_mode=="fiml") {
  msg(sprintf("Covariance SPD min eigen=%.2e (SEM set) → switching to missing='listwise' for stability.", chk0$min_eig))
  miss_mode <- "listwise"
} else {
  msg(sprintf("Covariance SPD check (SEM set): min eigen=%s (neg=%s); using missing='%s'.",
              format(chk0$min_eig, digits=3, scientific=TRUE), chk0$n_neg, miss_mode))
}

# 7) Preflight per outcome (keeps your earlier early failures)
for (oc in SEM_OUTCOMES) {
  try(preflight_sem(df = Z, mediators = SEM_MEDIATORS, genes = genes_union,
                    outcome = oc, tag = paste0("prelavaan_", oc)), silent = FALSE)
}


## ================================================================
## 7) SEM + trimming (respect --include_guild)
## ================================================================

build_model_txt <- function(adj) {
  # gene terms (instrument blocks) that survived so far
  eff <- lapply(names(adj), function(lhs) setdiff(genes_for[[lhs]], dropped_gene_rhs[[lhs]]))
  names(eff) <- names(adj)
  
  reg_lines <- character(0)
  for (lhs in names(adj)) {
    # labeled TT edges so the trimmer can drop by label; genes stay unlabeled
    tt_terms <- if (length(adj[[lhs]])) {
      vapply(
        adj[[lhs]],
        function(rhs) {
          lab <- label_for(lhs, rhs, if (lhs %in% SEM_OUTCOMES) "g" else "d")
          sprintf("%s*%s", lab, rhs)   # <- now supplies BOTH %s
        },
        character(1)
      )
    } else character(0)
    
    
    gene_terms <- setdiff(eff[[lhs]], adj[[lhs]])
    rhs_all <- c(tt_terms, gene_terms)
    if (length(rhs_all)) reg_lines <- c(reg_lines, paste(lhs, "~", paste(rhs_all, collapse = " + ")))
  }
  
  cov_block <- if (length(COV_PAIRS)) paste0(COV_PAIRS, " ;") else character(0)
  
  paste(
    "# --- Directional regressions (component-driven) ---",
    paste(paste0(reg_lines, " ;"), collapse = "\n"),
    if (length(cov_block)) "# --- Residual covariances ---" else NULL,
    if (length(cov_block)) paste(cov_block, collapse = "\n") else NULL,
    sep = "\n"
  )
}

fit_lavaan <- function(txt, data, missing_mode) {
  txt <- gsub("\r\n?", "\n", txt, perl = TRUE)
  txt <- gsub("\\+\\s*(?=\\r?\\n|$)", "", txt, perl = TRUE)  # '+ <EOL>'
  txt <- gsub("\\+\\s*;", " ;", txt, perl = TRUE)            # '+ ;' -> ' ;'
  txt <- gsub("(?<=\\n)\\s*\\+\\s*", "", txt, perl = TRUE)   # line starts '+'
  sem(model = txt, data = data, estimator="MLR",
      missing = if (missing_mode=="fiml") "fiml" else "listwise",
      se="robust", rstarts=5, control=list(iter.max=2000))
}

paths_tbl <- function(fit) {
  parameterEstimates(fit, standardized=TRUE) %>%
    as_tibble() %>% dplyr::filter(op=="~") %>%
    mutate(fdr = p.adjust(pvalue, method="BH"),
           sig = fdr < opt$fdr_level,
           label = ifelse(label=="", paste0("u_", dplyr::row_number()), label),
           formula = paste(lhs, "~", rhs))
}

qr_rescue_once <- function() {
  if (!nrow(Z)) return(FALSE)
  cc_cols <- unique(c(genes_union, SEM_TARGETS))
  cc <- stats::complete.cases(Z[, cc_cols, drop=FALSE])
  if (sum(cc) < 2) return(FALSE)
  Xcols <- c(genes_union, SEM_TARGETS)
  dep <- qr_drop_dependent_genes(Z[cc, , drop=FALSE], Xcols, SEM_TARGETS, genes_union)
  if (!length(dep)) return(FALSE)
  
  msg("QR rescue: dropping linearly dependent gene columns on listwise-complete rows (SEM set): ",
      paste(dep, collapse=", "))
  genes_for$Bact     <<- setdiff(genes_for$Bact,     dep)
  genes_for$Prev     <<- setdiff(genes_for$Prev,     dep)
  for (oc in SEM_OUTCOMES) genes_for[[oc]] <<- setdiff(genes_for[[oc]], dep)
  
  return(TRUE)   # <-- missing before
}


safe_initial_fit <- function(model_txt_local) {
  fit_try <- try(fit_lavaan(model_txt_local, Z, miss_mode), silent=TRUE)
  if (!inherits(fit_try, "try-error")) return(fit_try)
  emsg <- conditionMessage(attr(fit_try, "condition"))
  if (grepl("not positive-definite", emsg, ignore.case=TRUE)) {
    msg("lavaan error: sample covariance not positive-definite. Attempting QR rescue...")
    if (qr_rescue_once()) {
      fit_try2 <- try(fit_lavaan(build_model_txt(adj_tt), Z, miss_mode), silent=TRUE)
      if (!inherits(fit_try2, "try-error")) return(fit_try2)
      emsg2 <- conditionMessage(attr(fit_try2, "condition"))
      if (grepl("not positive-definite", emsg2, ignore.case=TRUE) && opt$jitter_sd > 0) {
        msg("Still NPD; applying tiny jitter (sd=", opt$jitter_sd, ") to genes and refitting once...")
        set.seed(12345)
        Z <<- maybe_jitter(Z, genes_union)
        return(fit_lavaan(build_model_txt(adj_tt), Z, miss_mode))
      } else stop(emsg2)
    } else if (opt$jitter_sd > 0) {
      msg("No dependent genes found via QR; applying tiny jitter (sd=", opt$jitter_sd, ") and refitting once...")
      set.seed(12345)
      Z <<- maybe_jitter(Z, genes_union)
      return(fit_lavaan(build_model_txt(adj_tt), Z, miss_mode))
    } else stop(emsg)
  } else stop(emsg)
}

## Robust label-based remover with regex fallback; never leaves '+' before a newline
## --- drop labeled terms on "~" lines; keep lines separate; keep ';' if present
drop_terms_from_model_by_label <- function(model_txt, labels) {
  labs <- unique(trimws(labels))
  if (!length(labs)) return(list(text = model_txt, removed = 0L, emptied = character(0)))
  
  txt <- gsub("\r\n?", "\n", model_txt, perl = TRUE)
  lines <- unlist(strsplit(txt, "\n", fixed = TRUE))
  removed_count <- 0L
  emptied_lhs <- character(0)
  
  for (i in seq_along(lines)) {
    li <- lines[i]
    if (!grepl("~", li, fixed = TRUE)) next
    
    had_sc <- grepl(";\\s*$", li)                    # remember semicolon
    lhs    <- sub("\\s*~.*$", "", li)
    rhs    <- sub("^\\s*[^~]+~", "", li)
    rhs    <- sub(";\\s*$", "", rhs)                 # strip trailing ';' BEFORE tokenizing
    
    toks <- strsplit(rhs, "\\+", perl = TRUE)[[1]]
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    if (!length(toks)) next
    
    # tokens look like "label*var" or just "var"
    tok_lab <- ifelse(grepl("\\*", toks, fixed = TRUE), sub("\\*.*$", "", toks), "")
    keep    <- !(tok_lab %in% labs)
    
    removed_here <- sum(!keep)
    if (removed_here > 0) {
      removed_count <- removed_count + removed_here
      toks <- toks[keep]
      
      if (length(toks)) {
        rhs_new <- paste(toks, collapse = " + ")
        lines[i] <- paste0(lhs, " ~ ", rhs_new, if (had_sc) " ;" else "")
      } else {
        # leave an explicitly empty RHS; no semicolon — the cascade will rebuild this line
        lines[i] <- paste0(lhs, " ~ ")
        emptied_lhs <- c(emptied_lhs, trimws(lhs))
      }
    }
  }
  
  # Whole-text cleanup (defensive): kill '+' at EOL or before ';' and normalize spacing
  txt1 <- paste(lines, collapse = "\n")
  txt1 <- gsub("\\+\\s*(?=\\r?\\n|$)", "", txt1, perl = TRUE)     # "+ <EOL>"
  txt1 <- gsub("\\+\\s*;", " ;", txt1, perl = TRUE)               # "+ ;" -> " ;"
  txt1 <- gsub("(?<=\\n)\\s*\\+\\s*", "", txt1, perl = TRUE)      # line begins with '+'
  # normalize "~" and "+" spacing on each line
  lines <- unlist(strsplit(txt1, "\n", fixed = TRUE))
  for (i in seq_along(lines)) {
    if (grepl("~", lines[i], fixed = TRUE)) {
      lines[i] <- sub("^\\s*([^~]+?)\\s*~\\s*", "\\1 ~ ", lines[i], perl = TRUE)
      lines[i] <- gsub("\\s*\\+\\s*", " + ", lines[i], perl = TRUE)
      lines[i] <- sub("\\s*\\+\\s*$", "", lines[i], perl = TRUE)
      lines[i] <- sub("^\\s*\\+\\s*", "", lines[i], perl = TRUE)
    }
  }
  
  list(
    text    = paste(lines, collapse = "\n"),
    removed = removed_count,
    emptied = unique(emptied_lhs)
  )
}

## --- remove any existing ':=' lines
strip_defined_block <- function(txt) {
  lines <- unlist(strsplit(gsub("\r\n?", "\n", txt, perl = TRUE), "\n", fixed = TRUE))
  keep  <- !grepl("^\\s*\\w+\\s*:=", lines)
  paste(lines[keep], collapse = "\n")
}

## --- list names of defined parameters currently in the text
list_defined_names <- function(txt) {
  lines <- unlist(strsplit(gsub("\r\n?", "\n", txt, perl = TRUE), "\n", fixed = TRUE))
  nm <- sub("^\\s*(\\w+)\\s*:=.*$", "\\1", grep("^\\s*\\w+\\s*:=", lines, value = TRUE))
  trimws(nm)
}

## --- drop full ':=' lines by their LHS names
drop_defined_by_name <- function(txt, names) {
  if (!length(names)) return(txt)
  lines <- unlist(strsplit(gsub("\r\n?", "\n", txt, perl = TRUE), "\n", fixed = TRUE))
  pat <- paste0("^\\s*(", paste(paste0("\\Q", names, "\\E"), collapse="|"), ")\\s*:=")
  keep <- !grepl(pat, lines, perl = TRUE)
  paste(lines[keep], collapse = "\n")
}

## --- after dropping some INDs, remove references to missing INDs from TOTALs; drop empty TOTALs
prune_totals_against_existing_inds <- function(txt) {
  lines <- unlist(strsplit(gsub("\r\n?", "\n", txt, perl = TRUE), "\n", fixed = TRUE))
  defined <- list_defined_names(txt)
  have_ind <- defined[startsWith(defined, "IND_")]
  
  for (i in seq_along(lines)) {
    li <- lines[i]
    if (!grepl("^\\s*TOTAL_\\w+\\s*:=", li)) next
    rhs <- sub("^\\s*TOTAL_\\w+\\s*:=\\s*", "", li)
    toks <- trimws(strsplit(rhs, "\\+", perl = TRUE)[[1]])
    toks <- toks[nzchar(toks)]
    # keep DIR_* always; keep IND_* only if they still exist
    keep <- vapply(toks, function(t) startsWith(t, "DIR_") || t %in% have_ind, logical(1))
    toks <- toks[keep]
    if (length(toks)) {
      rhs_new <- paste(toks, collapse = " + ")
      lines[i] <- sub(":=.*$", paste0(":= ", rhs_new), li)
    } else {
      lines[i] <- ""  # drop empty TOTAL
    }
  }
  paste(lines[nzchar(lines)], collapse = "\n")
}

## --- EFFECTS (:=) utilities -----------------------------------------------

# build a fresh effects block from the current (trimmed) fit
# - IND_<gene>_<nn> := product of labels along each simple path gene -> outcome
# - TOTAL_<gene>    := sum( all IND_<gene>_*  + direct label (if any) )
# REPLACE your build_effects_block_from_fit() with this:
build_effects_block_from_fit <- function(fit, sources, outcome) {
  pe_dir <- parameterEstimates(fit, standardized = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::filter(op == "~", nzchar(label)) %>%
    dplyr::select(from = rhs, to = lhs, label)
  
  if (!nrow(pe_dir)) return(NULL)
  
  g <- igraph::graph_from_data_frame(pe_dir, directed = TRUE)
  if (!(outcome %in% igraph::V(g)$name)) return(NULL)
  
  # fast label lookup
  key_edge <- paste(pe_dir$from, pe_dir$to, sep = "->")
  lab_edge <- pe_dir$label
  lab_for  <- function(u, v) lab_edge[match(paste(u, v, sep = "->"), key_edge)]
  
  defs <- character()
  
  for (src in sources) {
    if (!(src %in% igraph::V(g)$name)) next
    
    # all simple paths src -> outcome
    aps <- igraph::all_simple_paths(g, from = src, to = outcome, mode = "out")
    
    # keep ONLY *indirect* paths (≥ 2 edges = ≥ 3 nodes)
    indirect_paths <- Filter(function(p) length(p) >= 3, aps)
    
    # build IND names with outcome scope
    ind_names <- character(0)
    if (length(indirect_paths)) {
      path_labs <- lapply(indirect_paths, function(p) {
        nodes <- names(p)
        vapply(seq_len(length(nodes) - 1L), function(k) lab_for(nodes[k], nodes[k + 1L]), character(1))
      })
      keep_idx <- which(vapply(path_labs, function(x) all(nzchar(x)), logical(1)))
      if (length(keep_idx)) {
        ord <- order(vapply(path_labs[keep_idx], paste, collapse = "*", character(1)))
        keep_idx <- keep_idx[ord]
      }
      for (j in keep_idx) {
        labs <- path_labs[[j]]
        nm   <- sprintf("IND_%s_%s_%02d", outcome, src, length(ind_names) + 1L)
        ind_names <- c(ind_names, nm)
        defs <- c(defs, sprintf("%s := %s", nm, paste(labs, collapse = " * ")))
      }
    }
    
    # direct label (if any) is NOT an IND; include it only in TOTAL once
    dir_lab <- lab_for(src, outcome)
    total_terms <- c(ind_names, if (!is.na(dir_lab) && nzchar(dir_lab)) dir_lab else NULL)
    
    if (length(total_terms)) {
      defs <- c(defs, sprintf("TOTAL_%s_%s := %s", outcome, src, paste(total_terms, collapse = " + ")))
    }
  }
  
  if (!length(defs)) return(NULL)
  paste("# --- Indirect and total effects ---", paste(defs, collapse = "\n"), sep = "\n")
}


# drop ':=' rows by name on the LHS (exact match before ':=')
drop_defined_params_by_name <- function(model_txt, lhs_names) {
  lhs_names <- unique(lhs_names); if (!length(lhs_names)) return(list(text = model_txt, dropped = character(0)))
  model_txt <- gsub("\r\n?", "\n", model_txt, perl = TRUE)
  lines <- unlist(strsplit(model_txt, "\n", fixed = TRUE))
  keep <- rep(TRUE, length(lines)); dropped <- character(0)
  for (nm in lhs_names) {
    idx <- grep(paste0("^\\s*", nm, "\\s*:="), lines, perl = TRUE)
    if (length(idx)) { keep[idx] <- FALSE; dropped <- c(dropped, nm) }
  }
  list(text = paste(lines[keep], collapse = "\n"), dropped = unique(dropped))
}

# after dropping IND_* rows, remove their mentions inside TOTAL_* and drop empty TOTALs
prune_totals_after_ind_drop <- function(model_txt, ind_removed) {
  if (!length(ind_removed)) return(model_txt)
  model_txt <- gsub("\r\n?", "\n", model_txt, perl = TRUE)
  lines <- unlist(strsplit(model_txt, "\n", fixed = TRUE))
  for (i in seq_along(lines)) {
    if (!grepl("^\\s*TOTAL_", lines[i])) next
    lhs <- sub("\\s*:=.*$", "", lines[i])
    rhs <- sub("^\\s*[^:]+:=\\s*", "", lines[i])
    toks <- trimws(strsplit(rhs, "\\+", perl = TRUE)[[1]])
    toks <- toks[nzchar(toks)]
    toks <- toks[!(toks %in% ind_removed)]
    if (length(toks)) {
      rhs_new <- paste(toks, collapse = " + ")
      lines[i] <- paste0(lhs, " := ", rhs_new)
    } else {
      lines[i] <- NA_character_  # drop empty TOTAL
    }
  }
  paste(lines[!is.na(lines)], collapse = "\n")
}

## >>> when a target’s equation empties, cascade-remove it from downstream targets
enforce_after_line_empties <- function(model_txt, adj, order_vec, exogenous_ok = character(0)) {
  lines <- unlist(strsplit(model_txt, "\n", fixed = TRUE))
  
  # auto-whitelist nodes that only appear on RHS and have no "~" line
  rhs_nodes    <- unique(unlist(adj, use.names = FALSE))
  exogenous_ok <- unique(c(exogenous_ok, setdiff(rhs_nodes, names(adj))))
  
  lhs_present <- unique(sub("\\s*~.*$", "", grep("~", lines, value = TRUE)))
  
  is_empty_present <- function(lhs) {
    i <- grep(paste0("^\\s*", lhs, "\\s*~"), lines)
    if (!length(i)) return(FALSE)
    rhs <- trimws(sub("^\\s*[^~]+~", "", lines[i[1]]))
    !nzchar(gsub("\\++|\\s+|;", "", rhs))
  }
  
  empties_present <- lhs_present[vapply(lhs_present, is_empty_present, logical(1))]
  missing_lhs     <- setdiff(order_vec, lhs_present)
  
  # do NOT treat legitimate exogenous nodes as empties
  empties_present <- setdiff(empties_present, exogenous_ok)
  missing_lhs     <- setdiff(missing_lhs,     exogenous_ok)
  
  empties <- unique(c(empties_present, missing_lhs))
  if (!length(empties)) return(list(txt = model_txt, adj = adj, empties = character(0)))
  
  # remove each emptied/missing node as a predictor from all downstream equations
  for (src in empties) {
    to_idx <- match(src, order_vec)
    if (is.na(to_idx)) next
    downstream <- order_vec[seq.int(to_idx, length(order_vec))]
    for (tgt in downstream) if (!is.null(adj[[tgt]])) adj[[tgt]] <- setdiff(adj[[tgt]], src)
  }
  # prevent re-adding genes for those nodes (uses your globals)
  for (src in empties) {
    if (!is.null(genes_for[[src]]))        genes_for[[src]]        <<- character(0)
    if (!is.null(dropped_gene_rhs[[src]])) dropped_gene_rhs[[src]] <<- unique(dropped_gene_rhs[[src]])
  }
  
  list(txt = build_model_txt(adj), adj = adj, empties = empties)
}


## Build initial model and fit
msg("COV_PAIRS type: ", typeof(COV_PAIRS), " | class: ", paste(class(COV_PAIRS), collapse=","))
stopifnot(is.null(COV_PAIRS) || is.character(COV_PAIRS))

if (is.function(COV_PAIRS)) {
  warning("COV_PAIRS was a function; resetting to empty character vector.")
  COV_PAIRS <- character(0)
}
COV_PAIRS <- base::unique(c(COV_PAIRS, detect_resid_covariances(
  Z, adj_tt,
  candidates = intersect(c("Bact","Prev","Guild","Beta","Alpha", SEM_OUTCOMES), SEM_TARGETS),
  alpha = 0.05, thr = 0.10
)))

model_txt <- build_model_txt(adj_tt)
writeLines(model_txt, con=file.path(opt$outdir, "model_initial.lav"))
fit <- safe_initial_fit(model_txt)
log_lavaan(fit, iter = 0L, stage = "init", tag = "initial")
if (!lavInspect(fit, "converged")) {
  if (isTRUE(opt$jitter_sd > 0)) {
    msg("Initial fit failed; applying small jitter to genes and refitting once...")
    set.seed(12345)
    Z <- maybe_jitter(Z, genes_union)
    fit <- fit_lavaan(model_txt, Z, miss_mode)
    log_lavaan(fit, iter = 0L, stage = "init", tag = "post_jitter")
  }
  if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
    msg("Initial fit did not converge under FIML; refitting once with listwise...")
    fit <- fit_lavaan(model_txt, Z, "listwise")
    log_lavaan(fit, iter = 0L, stage = "init", tag = "fallback_listwise")
  }
}


## --- Trimming loop (phased: edges -> IND -> TOTAL) -------------------------

## helpers used only in this section
cascade_remove <- function(adj, from, to, order_vec) {
  adj[[to]] <- setdiff(adj[[to]], from)  # only remove from -> to
  adj
}

rebuild_after_gene_drop <- function(adj_now) {
  # Persist gene-level drops per LHS so rebuilds never re-add them
  for (lhs in names(dropped_gene_rhs)) {
    if (length(dropped_gene_rhs[[lhs]]) > 0) {
      genes_for[[lhs]] <<- setdiff(genes_for[[lhs]], dropped_gene_rhs[[lhs]])
    }
  }
  build_model_txt(adj_now)
}

rebuild_after_tt_drop <- function(adj_now) {
  exposures_in_dag <- setdiff(names(adj_now), SEM_OUTCOMES)
  exogenous_now    <- setdiff(SEM_TARGETS, names(adj_now))  # e.g., Bact, Prev
  
  COV_PAIRS <<- unique(c(
    COV_FIXED,                     # <- always keep skeleton covariances
    COV_PAIRS,                     # <- keep whatever we already had
    detect_resid_covariances(
      Z, adj_now,
      candidates = unique(c(exposures_in_dag, SEM_OUTCOMES, exogenous_now)),
      alpha = 0.05, thr = 0.10, forbid_parent_child = TRUE
    )
  ))
  
  build_model_txt(adj_now)
}





iter <- 0L
stage <- "edges"   # then "indirect", then "total"
droplog_path <- file.path(opt$outdir, "lavaan_iter_drop_log.tsv")
if (file.exists(droplog_path)) unlink(droplog_path)

## one-time prune for empty/missing LHS, then refit
pruned0 <- enforce_after_line_empties(model_txt, adj_tt, targets_order, EXOGENOUS_OK)
if (length(pruned0$empties)) {
  msg("[SEM] Nodes removed due to empty/missing equations: ", paste(pruned0$empties, collapse = ", "))
  adj_tt <- pruned0$adj
  exposures_in_dag <- setdiff(names(adj_tt), SEM_OUTCOMES)
  COV_PAIRS <- detect_resid_covariances(
    Z, adj_tt, candidates = unique(c(exposures_in_dag, SEM_OUTCOMES)),
    alpha = 0.05, thr = 0.10
  )
  model_txt <- build_model_txt(adj_tt)
  fit <- fit_lavaan(model_txt, Z, miss_mode)
  log_lavaan(fit, iter = 0L, stage = "init", tag = "post_prune")
}



repeat {
  iter <- iter + 1L
  if (iter >= opt$max_iter) { msg("[SEM] Reached --max_iter."); break }
  if (iter > 20)          { msg("[SEM] Safety stop at 20 iterations."); break }
  
  writeLines(model_txt, con = file.path(opt$outdir, sprintf("model_iter_%02d.lav", iter)))
  msg(sprintf("[SEM] Iteration %d (stage: %s)...", iter, stage))
  msg(sprintf("[SEM] Model BEFORE trimming (iter %d):\n%s", iter, model_txt))
  
  if (!lavInspect(fit, "converged")) {
    msg("[SEM] Model not converged; stopping trimming.")
    break
  }
  
  did_anything <- FALSE
  
  ## ====================== PHASE 1: EDGES (labels) =========================
  if (stage == "edges") {
    pe <- paths_tbl(fit)  # lhs, rhs, label, pvalue, fdr, sig
    
    ## 1) HOST-GENE --> target drops (by label), then rebuild from state
    is_gene_rhs <- pe$rhs %in% genes_union
    kill_gene <- pe %>%
      dplyr::filter(is_gene_rhs & !sig & nzchar(label)) %>%
      arrange(desc(fdr))
    
    if (nrow(kill_gene)) {
      to_drop_labels <- unique(kill_gene$label)
      msg("[SEM] Gene-labels to drop: ", paste(to_drop_labels, collapse = ", "))
      
      # persist dropped genes per LHS so rebuild never re-adds them
      by_lhs <- split(kill_gene$rhs, kill_gene$lhs)
      for (nm in names(by_lhs)) {
        dropped_gene_rhs[[nm]] <- unique(c(dropped_gene_rhs[[nm]], by_lhs[[nm]]))
      }
      
      # log
      data.table::fwrite(
        data.table(iter=iter, stage="edges", type="gene_label",
                   dropped=paste(to_drop_labels, collapse=",")),
        file=droplog_path, sep="\t", append=file.exists(droplog_path)
      )
      
      # rebuild and refit
      model_txt <- rebuild_after_gene_drop(adj_tt)
      fit <- fit_lavaan(model_txt, Z, miss_mode)
      log_lavaan(fit, iter, stage, tag = "post_gene_drop")
      if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
        msg("[SEM] Refit with listwise after gene-drops...")
        fit <- fit_lavaan(model_txt, Z, "listwise")
        log_lavaan(fit, iter, stage, tag = "fallback_listwise")
      }
      msg(sprintf("[SEM] Model AFTER gene-drop rebuild (iter %d):\n%s", iter, model_txt))
      did_anything <- TRUE
    }
    
    ## 2) TARGET -> TARGET drops (by label), update DAG, rebuild
    pe <- paths_tbl(fit)
    is_target_rhs <- pe$rhs %in% SEM_MEDIATORS
    kill_tt <- pe %>%
      dplyr::filter(is_target_rhs & !sig & nzchar(label)) %>%
      arrange(desc(fdr))
    
    if (nrow(kill_tt)) {
      to_drop_tt <- unique(kill_tt$label)
      msg("[SEM] TT-labels to drop: ", paste(to_drop_tt, collapse = ", "))
      
      # Update DAG adjacency by removing edges (from -> to) *downstream* and cascading
      dropped_pairs <- kill_tt %>% dplyr::transmute(from=rhs, to=lhs)
      for (k in seq_len(nrow(dropped_pairs))) {
        adj_tt <- cascade_remove(adj_tt, dropped_pairs$from[k], dropped_pairs$to[k], targets_order)
      }
      
      # log
      data.table::fwrite(
        data.table(iter=iter, stage="edges", type="target_target_label",
                   dropped=paste(to_drop_tt, collapse=",")),
        file=droplog_path, sep="\t", append=file.exists(droplog_path)
      )
      
      # rebuild and refit
      model_txt <- rebuild_after_tt_drop(adj_tt)  # recomputes COV_PAIRS internally
      fit <- fit_lavaan(model_txt, Z, miss_mode)
      log_lavaan(fit, iter, stage, tag = "post_gene_drop")
      if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
        msg("[SEM] Refit with listwise after TT-drops...")
        fit <- fit_lavaan(model_txt, Z, "listwise")
        log_lavaan(fit, iter, stage, tag = "fallback_listwise")
      }
      msg(sprintf("[SEM] Model AFTER TT-drop rebuild (iter %d):\n%s", iter, model_txt))
      did_anything <- TRUE
    }
    
    # If nothing else to drop at the edges layer, lock edges and inject ':=' effects:
    if (!did_anything) {
      # Prune once if there are empties, then rebuild/refit
      EXOGENOUS_OK <- compute_exogenous_ok(adj_tt, genes_for, SEM_TARGETS)
      msg("Exogenous (by model skeleton): ", paste(EXOGENOUS_OK, collapse = ", "))
      final_prune  <- enforce_after_line_empties(model_txt, adj_tt, targets_order, EXOGENOUS_OK)
      if (length(final_prune$empties)) {
        COV_PAIRS <- unique(c(
          COV_FIXED,
          COV_PAIRS,
          detect_resid_covariances(
            Z, adj_tt,
            candidates = unique(c("Bact","Prev", SEM_OUTCOMES)),
            alpha = 0.05, thr = 0.10
          )
        ))
        model_txt <- build_model_txt(adj_tt)
        fit <- fit_lavaan(model_txt, Z, miss_mode)
        msg("[SEM] Final pre-lock prune of empty/missing nodes: ",
            paste(final_prune$empties, collapse = ", "))
        adj_tt    <- final_prune$adj
        model_txt <- final_prune$txt
        fit <- fit_lavaan(model_txt, Z, miss_mode)
        log_lavaan(fit, iter, stage, tag = "post_gene_drop")
        if (!lavInspect(fit, "converged") && miss_mode == "fiml") {
          fit <- fit_lavaan(model_txt, Z, "listwise")
          log_lavaan(fit, iter, stage, tag = "fallback_listwise")
        }
        msg(sprintf("[SEM] Model AFTER prune (iter %d):\n%s", iter, model_txt))
      }
      
      msg("[SEM] No more non-significant labeled paths. Locking edges and building ':=' effects...")
      blocks <- character(0)
      
      # named exogenous by current DAG (nodes that appear only on RHS)
      # You already computed EXOGENOUS_OK above via compute_exogenous_ok()
      exo_named <- intersect(EXOGENOUS_OK, c("Bact","Prev"))
      
      for (oc in SEM_OUTCOMES) {
        # gene instruments are exogenous
        if (length(genes_union)) {
          bg <- build_effects_block_from_fit(fit, genes_union, outcome = oc)
          if (!is.null(bg)) blocks <- c(blocks, bg)
        }
        # named exogenous (Bact/Prev) only
        if (length(exo_named)) {
          bex <- build_effects_block_from_fit(fit, exo_named, outcome = oc)
          if (!is.null(bex)) blocks <- c(blocks, bex)
        }
      }
      
      
      model_txt <- strip_defined_block(model_txt)
      if (length(blocks)) model_txt <- paste(model_txt, paste(blocks, collapse = "\n"), sep = "\n")
      
      writeLines(model_txt, con = file.path(opt$outdir, sprintf("model_iter_%02d_with_effects.lav", iter)))
      fit <- fit_lavaan(model_txt, Z, miss_mode)
      log_lavaan(fit, iter, stage, tag = "post_gene_drop")
      if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
        msg("[SEM] Refit with listwise after injecting effects...")
        fit <- fit_lavaan(model_txt, Z, "listwise")
        log_lavaan(fit, iter, stage, tag = "fallback_listwise")
      }
      
      stage <- "indirect"
      next
    }
    
    ## ====================== PHASE 2: IND (defined params) ====================
  } else if (stage == "indirect") {
    pe_all <- parameterEstimates(fit, standardized = TRUE)
    if (nrow(pe_all)) {
      pe_ind <- pe_all %>%
        as_tibble() %>%
        dplyr::filter(op == ":=", grepl("^IND_", lhs)) %>%
        dplyr::mutate(fdr = p.adjust(pvalue, method="BH"), sig = fdr < opt$fdr_level)
      
      kill_ind <- pe_ind %>% dplyr::filter(!sig)
      if (nrow(kill_ind)) {
        ind_names <- unique(kill_ind$lhs)
        msg("[SEM] IND to drop (BH-FDR): ", paste(ind_names, collapse=", "))
        
        data.table::fwrite(
          data.table(iter=iter, stage="indirect", type="IND_defined",
                     dropped=paste(ind_names, collapse=",")),
          file=droplog_path, sep="\t", append=file.exists(droplog_path)
        )
        
        # Drop those ':=' rows by name and repair TOTAL_* rows referencing them
        dd <- drop_defined_params_by_name(model_txt, ind_names)
        model_txt <- prune_totals_after_ind_drop(dd$text, dd$dropped)
        
        fit <- fit_lavaan(model_txt, Z, miss_mode)
        log_lavaan(fit, iter, stage, tag = "post_gene_drop")
        if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
          msg("[SEM] Refit with listwise after IND drops...")
          fit <- fit_lavaan(model_txt, Z, "listwise")
          log_lavaan(fit, iter, stage, tag = "fallback_listwise")
        }
        
        msg(sprintf("[SEM] Model AFTER IND-drop rebuild (iter %d):\n%s", iter, model_txt))
        did_anything <- TRUE
      }
    }
    
    if (!did_anything) {
      msg("[SEM] No more non-significant IND effects. Proceeding to TOTAL.")
      stage <- "total"
      next
    }
    
    ## ====================== PHASE 3: TOTAL (defined params) ==================
  } else if (stage == "total") {
    pe_all <- parameterEstimates(fit, standardized = TRUE)
    if (nrow(pe_all)) {
      pe_tot <- pe_all %>%
        as_tibble() %>%
        dplyr::filter(op == ":=", grepl("^TOTAL_", lhs)) %>%
        dplyr::mutate(fdr = p.adjust(pvalue, method="BH"), sig = fdr < opt$fdr_level)
      
      kill_tot <- pe_tot %>% dplyr::filter(!sig)
      if (nrow(kill_tot)) {
        tot_names <- unique(kill_tot$lhs)
        msg("[SEM] TOTAL to drop (BH-FDR): ", paste(tot_names, collapse=", "))
        
        data.table::fwrite(
          data.table(iter=iter, stage="total", type="TOTAL_defined",
                     dropped=paste(tot_names, collapse=",")),
          file=droplog_path, sep="\t", append=file.exists(droplog_path)
        )
        
        dd <- drop_defined_params_by_name(model_txt, tot_names)
        model_txt <- dd$text
        
        fit <- fit_lavaan(model_txt, Z, miss_mode)
        log_lavaan(fit, iter, stage, tag = "post_gene_drop")
        if (!lavInspect(fit, "converged") && miss_mode=="fiml") {
          msg("[SEM] Refit with listwise after TOTAL drops...")
          fit <- fit_lavaan(model_txt, Z, "listwise")
          log_lavaan(fit, iter, stage, tag = "fallback_listwise")
        }
        
        msg(sprintf("[SEM] Model AFTER TOTAL-drop rebuild (iter %d):\n%s", iter, model_txt))
        did_anything <- TRUE
      }
    }
    
    if (!did_anything) {
      msg("[SEM] Trimming complete (edges, IND, TOTAL).")
      break
    }
  }
}


## ---- Ensure ':=' effects exist on the FINAL model (if possible) ----
if (lavInspect(fit, "converged")) {
  pe_dir <- parameterEstimates(fit, standardized = TRUE) %>%
    as_tibble() %>% dplyr::filter(op == "~", nzchar(label)) %>%
    dplyr::select(from = rhs, to = lhs, label)
  if (nrow(pe_dir)) {
    g <- igraph::graph_from_data_frame(pe_dir, directed = TRUE)
    avail <- igraph::V(g)$name
    eff_blocks <- character(0)
    
    # recompute exogenous from final adjacency
    EXOGENOUS_OK <- compute_exogenous_ok(adj_tt, genes_for, SEM_TARGETS)
    exo_named    <- intersect(EXOGENOUS_OK, c("Bact","Prev"))
    
    for (oc in SEM_OUTCOMES) if (oc %in% avail) {
      ebg <- build_effects_block_from_fit(fit, genes_union, outcome = oc)
      ebx <- if (length(exo_named)) build_effects_block_from_fit(fit, exo_named, outcome = oc) else NULL
      if (!is.null(ebg)) eff_blocks <- c(eff_blocks, ebg)
      if (!is.null(ebx)) eff_blocks <- c(eff_blocks, ebx)
    }
    
    model_txt <- strip_defined_block(model_txt)
    if (length(eff_blocks)) model_txt <- paste(model_txt, paste(eff_blocks, collapse = "\n"), sep = "\n")
  }
}

## Save model comparisons
fit_first <- try(fit_lavaan(readr::read_file(file.path(opt$outdir, "model_initial.lav")), Z, miss_mode), silent=TRUE)
fit_last  <- fit
fm_idx <- c("cfi","tli","rmsea","srmr","aic","bic")
F1 <- tryCatch(fitMeasures(fit_first, fm_idx), error=function(e) setNames(rep(NA_real_,length(fm_idx)), fm_idx))
FL <- tryCatch(fitMeasures(fit_last,  fm_idx), error=function(e) setNames(rep(NA_real_,length(fm_idx)), fm_idx))
data.table::fwrite(
  tibble::tibble(Model=c("Initial_SEM","Final_SEM"), CFI=c(F1["cfi"],FL["cfi"]), TLI=c(F1["tli"],FL["tli"]),
                 RMSEA=c(F1["rmsea"],FL["rmsea"]), SRMR=c(F1["srmr"],FL["srmr"]),
                 AIC=c(F1["aic"],FL["aic"]), BIC=c(F1["bic"],FL["bic"])),
  file=file.path(opt$outdir,"lavaan_fit_comparison.tsv"), sep="\t"
)
save(model_txt, fit_first, fit_last, collapse_log, file=file.path(opt$outdir,"lavaan_models.RData"))
writeLines(model_txt, con=file.path(opt$outdir,"model_final_trimmed.lav"))

## Mermaid (SEM)
msg("[SEM] Generating Mermaid diagram code...")
pe_final <- tryCatch(parameterEstimates(fit_last, standardized = TRUE), error = function(e) NULL)

if (!is.null(pe_final)) {
  pe_final <- as_tibble(pe_final) %>%
    dplyr::filter(op %in% c("~","~~")) %>%
    # drop variances (lhs==rhs) so we only keep covariances for '~~'
    dplyr::filter(op == "~" | (op == "~~" & lhs != rhs)) %>%
    # if you want CORRELATIONS on covariances (recommended for readability):
    dplyr::mutate(lbl = dplyr::case_when(
      op == "~"  ~ sprintf('"%0.2f"', round(std.all, 2)),  # standardized path coef
      op == "~~" ~ sprintf('"%0.2f"', round(std.all, 2))   # correlation (std.all)
      # if you prefer UNSTANDARDIZED covariance values, swap std.all -> est above
    ))
  
  ## de-duplicate symmetric covariances
  pe_cov <- pe_final %>%
    dplyr::filter(op == "~~") %>%
    dplyr::mutate(pair = paste(pmin(lhs, rhs), pmax(lhs, rhs), sep = "~~")) %>%
    dplyr::group_by(pair) %>%
    dplyr::slice(1) %>% dplyr::ungroup()
  
  pe_dir <- pe_final %>% dplyr::filter(op == "~")
  pe_final2 <- dplyr::bind_rows(pe_dir, pe_cov) %>%
    dplyr::select(lhs, rhs, op, lbl)
  
  nodes <- unique(c(pe_final2$lhs, pe_final2$rhs))
  node_group <- function(x) dplyr::case_when(
    x %in% genes_union               ~ "host_gene",
    x %in% c("Bact","Prev","Guild")  ~ "bacteria",
    x %in% c("Beta","Alpha")         ~ "community",
    x %in% c("Glucose","BMI")        ~ "phenotype",
    TRUE                             ~ "other"
  )
  N <- tibble::tibble(id = nodes, label = nodes, class = node_group(nodes)) %>%
    dplyr::filter(class != "other")
  
  edge_line <- function(lhs, op, rhs, lbl) {
    # regressions use -->|label| ; covariances use ---|label|
    sprintf("%s %s|%s| %s", rhs, ifelse(op == "~~", "---", "-->"), lbl, lhs)
  }
  E <- pe_final2 %>% dplyr::transmute(line = edge_line(lhs, op, rhs, lbl)) %>%
    dplyr::pull(line)
  
  mermaid <- c(
    "graph TB",
    "classDef bacteria    fill:#E0FFFF,stroke:#333;",
    "classDef community   fill:#98FB98,stroke:#333;",
    "classDef host_gene   fill:#FFE4E1,stroke:#333;",
    "classDef phenotype   fill:#FFD700,stroke:#333;",
    paste0(N$id, "[", N$label, "]:::", N$class),
    E
  )
  writeLines(paste(mermaid, collapse = "\n"),
             con = file.path(opt$outdir, "mermaid_final.mmd"))
  msg("[SEM] Mermaid code written to mermaid_final.mmd")
} else {
  msg("[SEM] Skipping Mermaid: no parameterEstimates available.")
}

## ---- blavaan validation (AFTER final lavaan model is finished) ----
if (lavInspect(fit_last, "converged")) {
  msg("[SEM] Refitting final trimmed model in blavaan... (weakly informative priors)")
  model_for_bayes <- model_txt  # already includes ':=' effects if available
  fit_bsem <- tryCatch(
    bsem(model    = model_for_bayes, data = Z,
         dp       = dpriors(alpha="normal(0,0.5)", beta="normal(0,0.5)"),
         cp       = "srs", target="stan", n.chains=4, burnin=2000, sample=5000,
         control  = list(adapt_delta=0.98, max_treedepth=14)),
    error=function(e){ msg("blavaan fit failed: ", e$message); NULL }
  )
  if (!is.null(fit_bsem)) {
    sum_bsem=summary(fit_bsem,
                     standardized = TRUE,
                     fit.measures = TRUE,
                     rsquare = TRUE,
                     ci = TRUE)
    write(sum_bsem,file=file.path(opt$outdir,"blavaan_fit_summary.txt"))
    save(fit_bsem, file=file.path(opt$outdir,"blavaan_fit.RData"))
    print(summary(fit_bsem, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE, ci=TRUE))
    pe_bsem <- parameterEstimates(fit_bsem, standardized = TRUE, ci = TRUE) %>%
      as_tibble() %>% dplyr::filter(op %in% c("~", ":="))
    lo_col <- dplyr::case_when(
      "ci.lower" %in% names(pe_bsem) ~ "ci.lower",
      "pi.lower" %in% names(pe_bsem) ~ "pi.lower",
      TRUE ~ NA_character_
    )
    hi_col <- dplyr::case_when(
      "ci.upper" %in% names(pe_bsem) ~ "ci.upper",
      "pi.upper" %in% names(pe_bsem) ~ "pi.upper",
      TRUE ~ NA_character_
    )
    if (!is.na(lo_col) && !is.na(hi_col)) {
      pe_bsem <- pe_bsem %>% dplyr::mutate(
        sig95 = (!!rlang::sym(lo_col) > 0 & !!rlang::sym(hi_col) > 0) |
          (!!rlang::sym(lo_col) < 0 & !!rlang::sym(hi_col) < 0)
      )
    } else {
      warning("Credible interval columns not found in parameterEstimates(); sig95 set to NA.")
      pe_bsem <- pe_bsem %>% dplyr::mutate(sig95 = NA)
    }
    data.table::fwrite(pe_bsem, file = file.path(opt$outdir, "blavaan_paths.tsv"), sep = "\t")
  } else msg("[SEM] Skipping blavaan outputs due to earlier error.")
} else {
  msg("[SEM] Skipping blavaan: lavaan did not converge; fix SEM first (see logs).")
}

## ================================================================
## Pairwise direct / indirect / total effect extractor
## - Works from a fit object (lavaan/blavaan) or from blavaan_paths.tsv
## - Optional 'through' filter (require one or more intermediates to be on the path)
## - Writes a neat TSV of path-by-path effects and a one-line summary
## ================================================================
suppressWarnings({
  library(dplyr); library(tibble); library(stringr); library(data.table); library(igraph)
})

# -- core worker: compute effects from a parameter-estimates table (pe)
.sem_pair_effects_from_pe <- function(pe,
                                      start, end,
                                      type = c("std.all","est"),
                                      through = NULL, through_all = FALSE) {
  type <- match.arg(type)
  req_cols <- c("lhs","op","rhs","est")
  stopifnot(all(req_cols %in% names(pe)))
  if (type == "std.all" && !("std.all" %in% names(pe))) {
    warning("std.all not found; falling back to 'est'.")
    type <- "est"
  }
  # keep only directed regressions
  Edf <- pe %>%
    dplyr::filter(.data$op == "~") %>%
    dplyr::transmute(from = trimws(.data$rhs),
                     to   = trimws(.data$lhs),
                     w    = if (type == "est") .data$est else .data$std.all,
                     #p    = dplyr::coalesce(.data$pvalue, NA_real_)
                     )
  if (!nrow(Edf)) {
    return(list(
      ok = FALSE, reason = "No '~' rows in parameter table.",
      direct = NA_real_, #direct_p = NA_real_,
      indirect = NA_real_, total = NA_real_,
      paths = tibble(Path = character(), Effect = numeric())
    ))
  }
  # build graph
  nodes <- unique(c(Edf$from, Edf$to, start, end))
  g <- igraph::graph_from_data_frame(Edf[,c("from","to")], directed = TRUE, vertices = nodes)
  igraph::E(g)$weight <- Edf$w
  #igraph::E(g)$pval   <- Edf$p
  
  # direct
  direct   <- NA_real_
  direct_p <- NA_real_
  if (igraph::are.connected(g, start, end)) {
    eid <- igraph::get.edge.ids(g, c(start, end))
    if (eid != 0) direct <- igraph::E(g)$weight[eid]
  }
  # pull p-value from parameter table (lavaan PE)
  dir_row <- pe %>% dplyr::filter(op == "~", lhs == end, rhs == start)
  if (nrow(dir_row)) direct_p <- dplyr::coalesce(dir_row$pvalue[1], NA_real_)
  
  
  # all simple paths start -> end
  if (!(start %in% igraph::V(g)$name) || !(end %in% igraph::V(g)$name)) {
    return(list(
      ok = FALSE, reason = sprintf("Either '%s' or '%s' not in graph.", start, end),
      direct = direct, #direct_p = direct_p,
      indirect = NA_real_, total = if (is.finite(direct)) direct else NA_real_,
      paths = tibble(Path = character(), Effect = numeric())
    ))
  }
  aps <- igraph::all_simple_paths(g, from = start, to = end, mode = "out")
  
  # optional 'through' filtering
  if (!is.null(through) && length(through)) {
    through <- unique(through)
    aps <- Filter(function(p) {
      pn <- names(p)
      if (!length(pn)) return(FALSE)
      if (through_all) all(through %in% pn) else any(through %in% pn)
    }, aps)
  }
  
  # compute path products
  path_tab <- tibble(Path = character(), Effect = numeric(), Length = integer())
  for (p in aps) {
    nodes <- names(p)
    if (length(nodes) < 2) next
    eff <- 1
    ok  <- TRUE
    for (k in seq_len(length(nodes) - 1L)) {
      eid <- igraph::get.edge.ids(g, c(nodes[k], nodes[k+1]))
      if (eid == 0) { ok <- FALSE; break }
      eff <- eff * igraph::E(g)$weight[eid]
    }
    if (ok) {
      path_tab <- add_row(path_tab,
                          Path = paste(nodes, collapse = " → "),
                          Effect = eff,
                          Length = length(nodes) - 1L)
    }
  }
  # indirect = sum of paths with ≥ 2 edges (Length >= 2)
  indirect <- path_tab %>% dplyr::filter(.data$Length >= 2) %>% dplyr::pull(.data$Effect) %>% sum(na.rm = TRUE)
  # total
  total <- (if (is.finite(direct)) direct else 0) + indirect
  
  list(
    ok = TRUE,
    direct   = if (is.finite(direct)) direct else NA_real_,
    direct_p = if (is.finite(direct_p)) direct_p else NA_real_,
    indirect = if (is.finite(indirect)) indirect else NA_real_,
    total    = if (is.finite(total)) total else NA_real_,
    paths    = path_tab %>% dplyr::arrange(dplyr::desc(abs(.data$Effect)))
  )
  
}

# -- wrapper: take a lavaan/blavaan fit object
sem_pair_effects_from_fit <- function(fit, start, end,
                                      type = c("std.all","est"),
                                      through = NULL, through_all = FALSE) {
  stopifnot(!missing(fit))
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  .sem_pair_effects_from_pe(pe, start, end, type = match.arg(type),
                            through = through, through_all = through_all)
}

# -- wrapper: take the blavaan_paths.tsv (optionally use ':=' rows if present)
sem_pair_effects_from_file <- function(tsv_path, start, end,
                                       type = c("std.all","est"),
                                       through = NULL, through_all = FALSE) {
  stopifnot(file.exists(tsv_path))
  D <- data.table::fread(tsv_path, sep = "\t", data.table = FALSE)
  res <- .sem_pair_effects_from_pe(D, start, end, type = match.arg(type),
                                   through = through, through_all = through_all)
  
  # Attach any defined parameters already present (lavaan ':=')
  Ddef <- D %>% dplyr::filter(.data$op == ":=")
  if (nrow(Ddef)) {
    tot_name <- sprintf("TOTAL_%s_%s", end, start)
    ind_re   <- sprintf("^IND_%s_%s_", end, start)
    # pick column for estimate
    est_col <- if (match.arg(type) == "std.all" && "std.all" %in% names(Ddef)) "std.all" else "est"
    
    tot_row <- Ddef %>% dplyr::filter(.data$lhs == tot_name) %>% dplyr::slice(1)
    ind_rows <- Ddef %>% dplyr::filter(grepl(ind_re, .data$lhs))
    
    res$defined_total   <- if (nrow(tot_row)) tot_row[[est_col]] else NA_real_
    #res$defined_total_p <- if (nrow(tot_row) && "pvalue" %in% names(tot_row)) tot_row$pvalue else NA_real_
    res$defined_indirect <- if (nrow(ind_rows)) sum(ind_rows[[est_col]], na.rm = TRUE) else NA_real_
    res$defined_individual_inds <- if (nrow(ind_rows)) {
      ind_rows %>% dplyr::transmute(name = .data$lhs,
                                    Effect = .data[[est_col]],
                                    #pvalue = dplyr::coalesce(.data$pvalue, NA_real_)
                                    )
    } else {
      tibble(name = character(), Effect = numeric()#, pvalue = numeric()
             )
    }
  }
  res
}

# -- convenience writer
write_pair_effects <- function(res, start, end, outdir, suffix = "") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  base <- paste0("pair_effects_", start, "_to_", end, suffix)
  # paths table
  data.table::fwrite(as.data.table(res$paths),
                     file = file.path(outdir, paste0(base, "_paths.tsv")),
                     sep = "\t")
  # one-liner summary
  summ <- tibble::tibble(
    start = start, end = end,
    direct  = res$direct,
    direct_p = if (!is.null(res$direct_p)) res$direct_p else NA_real_,
    indirect = res$indirect,
    total    = res$total,
    defined_total   = if (!is.null(res$defined_total)) res$defined_total else NA_real_,
    defined_total_p = if (!is.null(res$defined_total_p)) res$defined_total_p else NA_real_,
    defined_indirect = if (!is.null(res$defined_indirect)) res$defined_indirect else NA_real_
  )
  
  data.table::fwrite(as.data.table(summ),
                     file = file.path(outdir, paste0(base, "_summary.tsv")),
                     sep = "\t")
  invisible(summ)
}

### ------------------ Examples (safe to leave; no side effects) ------------------
### 1) From final lavaan fit (completely standardized):
#if (exists("fit_last")) {
#  ex1 <- sem_pair_effects_from_fit(fit_last, start = "Bact", end = "Glucose", type = "std.all")
#  write_pair_effects(ex1, "Bact", "Glucose", opt$outdir, suffix = "_lavaan")
#}

## 2) From blavaan_paths.tsv, requiring the path to include 'Beta' (through_all=FALSE -> any-of):
#bp_tsv <- file.path(opt$outdir, "blavaan_paths.tsv")
bp_tsv <- file.path(opt$outdir, "blavaan_paths.tsv")
#if (file.exists(bp_tsv)) {
#  ex2 <- sem_pair_effects_from_file(bp_tsv, start = "Bact", end = "Glucose",
#                                    type = "std.all", #through = c("Beta"), 
#                                    through_all = TRUE)
#  paths=ex2$paths
#  #write_pair_effects(ex2, "Prev", "Glucose", bp_tsv, suffix = "_blavaan")
#}

### 3) Another pair (intermediate to intermediate), e.g., Prev -> Alpha:
#if (exists("fit_last")) {
#  ex3 <- sem_pair_effects_from_fit(fit_last, start = "Prev", end = "Alpha", type = "std.all")
#  write_pair_effects(ex3, "Prev", "Alpha", opt$outdir, suffix = "_lavaan")
#}
## ===================================================================


msg("All done. Outputs in: ", opt$outdir)
