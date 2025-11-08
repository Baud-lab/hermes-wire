#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GENESIS)
  library(MVMR)
})

## ---------------- CLI ----------------
opt_list <- list(
  make_option("--selected_genes_rdata", type="character"),
  make_option("--snps_tsv",             type="character"),
  make_option("--trait",                type="character"),
  make_option("--resid_rda",            type="character"),
  make_option("--beta_rda",             type="character"),
  make_option("--alpha_rda",            type="character"),
  make_option("--clusters_rda",         type="character"),
  make_option("--phen_rda",             type="character"),
  make_option("--phen_glucose_col",     type="character", default="Glucose"),
  make_option("--phen_bmi_col",         type="character", default="BMI"),  # already present
  make_option("--dam_rds",              type="character"),
  make_option("--cage_rds",             type="character"),
  make_option("--loco_rds",             type="character"),     # RDS: chr -> GRM prefix
  make_option("--scan_lookup",          type="character"),     # TSV with scan_id, sample_label
  make_option("--genus",                type="character"),
  make_option("--other_genus",          type="character"),
  make_option("--trait_prefix",         type="character", default="g__"),
  make_option("--beta_trait",           type="character", default="beta__PD_PC1"),
  make_option("--alpha_trait",           type="character", default="alpha__PD_q2"),
  make_option("--guild_trait",         type="character", default="cluster__3"),
  make_option("--mr_enable",           type="logical",   default=TRUE),
  make_option("--mr_fstat_min",        type="double",    default=10),
  make_option("--mr_min_snps", type="integer", default=3, help="min IVs per outcome for MVMR"),
  make_option(
    "--outcomes",
    type    = "character",
    default = "genus_other",
    help    = "Comma-separated subset of outcomes to test (genus_other,beta,guild,glucose,bmi) or 'all'."
  ),
  make_option("--sig_metric", type="character", default="q",
              help="Use 'p' or 'q' to call significance for alpha/beta/direct/indirect"),
  make_option("--fdr_level",  type="double",    default=0.10,
              help="FDR threshold when --sig_metric=q"),
  make_option("--p_level",    type="double",    default=0.05,
              help="Nominal p-value threshold when --sig_metric=p"),
  make_option("--out",        type="character", default="mediation_results.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ---- robust boolean + numeric coercion ----
to_bool <- function(x) {
  if (isTRUE(x)) return(TRUE)
  vx <- tolower(as.character(x))
  vx %in% c("true","t","1","yes","y")
}
opt$mr_enable    <- to_bool(opt$mr_enable)
opt$mr_fstat_min <- as.numeric(opt$mr_fstat_min)

msg     <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
norm_id <- function(x) sub("_.*$", "", as.character(x))
`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a

## ---------------- sanity/files ----------------
required <- c("selected_genes_rdata","snps_tsv","resid_rda","beta_rda","alpha_rda","clusters_rda",
              "phen_rda","dam_rds","cage_rds","loco_rds","scan_lookup")
miss <- required[!file.exists(unlist(opt[required]))]
if (length(miss)) stop("Missing required file(s): ", paste(miss, collapse=", "))

## ---------------- lookup (to harmonize IDs) ----------------
LK <- fread(opt$scan_lookup)  # must have: scan_id, sample_label
if (!all(c("scan_id","sample_label") %in% names(LK)))
  stop("--scan_lookup must have columns: scan_id, sample_label")
LK[, scan_id := as.character(scan_id)]
LK[, sample_label_trim := norm_id(sample_label)]
lab2scan <- setNames(LK$scan_id,           LK$sample_label_trim)
scan2lab <- setNames(LK$sample_label_trim, LK$scan_id)

## ---------------- load selected genes ----------------
load(opt$selected_genes_rdata)   # needs: genotypes, correlations
if (!exists("genotypes"))   stop("selected_genes_rdata lacks object 'genotypes'")
if (!exists("correlations")) stop("selected_genes_rdata lacks object 'correlations'")

stopifnot("sample" %in% names(genotypes))
genotypes$sample <- norm_id(genotypes$sample)
genotypes <- genotypes[!duplicated(genotypes$sample), , drop = FALSE]

# SNPs to test for this fold
SNPs <- fread(opt$snps_tsv)
if (!("snp" %in% names(SNPs))) stop("--snps_tsv must contain a column named 'snp' (chr:pos)")
snps_fold <- unique(SNPs$snp)

## correlations: expected columns
need_cols <- c("Gene_name","gene","chr","pos","variant.id")
correlations$gene <- paste(correlations$chr,correlations$pos,sep=":")
if (!all(need_cols %in% names(correlations)))
  stop("correlations missing required columns: ", paste(setdiff(need_cols, names(correlations)), collapse=", "))
corr_fold <- correlations[correlations$gene %in% snps_fold, ]
if (nrow(corr_fold) == 0L) stop("No rows in 'correlations' match the fold SNPs")
gene_names <- intersect(corr_fold$Gene_name, setdiff(names(genotypes), "sample"))
if (!length(gene_names)) stop("No genotype columns in 'genotypes' match 'Gene_name' for this fold")
msg("Fold genes (after matching genotypes): ", length(gene_names))

## ---------------- residuals / beta / clusters / phenotype ----------------
if (opt$trait=="Taxa"){
  load(opt$resid_rda)
  genus_layer_idx <- length(residuals_qned_counts_objs) - 1L
  R_genus <- t(residuals_qned_counts_objs[[genus_layer_idx]])
} else if (opt$trait=="Guild") {
  load(opt$clusters_rda)
  R_genus <- t(residuals_qned_counts_clusters_objs[[1]])
} else if (opt$trait=="Beta") {
  load(opt$beta_rda)
  R_genus <- t(residuals_qned_counts_beta_objs[[1]])
} else if (opt$trait=="Alpha") {
  load(opt$alpha_rda)
  R_genus <- t(residuals_qned_counts_alpha_objs[[1]])
} else {
  stop("Unknown --trait: ", opt$trait)
}
rownames(R_genus) <- sapply(strsplit(rownames(R_genus), "_"), `[`, 1)

pick_col <- function(M, key) {
  cn <- colnames(M)
  exact <- which(cn == key)
  if (length(exact)) return(cn[exact[1]])
  st <- which(startsWith(cn, key))
  if (length(st)) return(cn[st[1]])
  NA_character_
}

## mediator + (optional) other-genus
med_key <- paste0(opt$trait_prefix, opt$genus)
med_col <- pick_col(R_genus, med_key)
if (is.na(med_col)) stop("Mediator trait not found in residuals: ", med_key)
mediator_v <- setNames(as.numeric(R_genus[, med_col]), rownames(R_genus))

genus_other_v <- NULL
if (opt$trait=="Taxa") {
  other_key <- paste0(opt$trait_prefix, opt$other_genus)
  other_col <- pick_col(R_genus, other_key)
  if (is.na(other_col)) stop("Other-genus trait not found in residuals: ", other_key)
  genus_other_v <- setNames(as.numeric(R_genus[, other_col]), rownames(R_genus))
  msg("Mediator col='", med_col, "' | Other-genus col='", other_col, "'")
} else {
  msg("Mediator col='", med_col, "'")
}

## beta
load(opt$beta_rda)
B <- as.matrix(residuals_qned_counts_beta_objs[[1]])
beta_row <- if (opt$beta_trait %in% c("beta__PD_PC1","Axis.1")) 5L else if (opt$beta_trait %in% c("beta__PD_PC2","Axis.2")) 6L else 5L
beta_v <- setNames(as.numeric(B[beta_row, ]), norm_id(colnames(B)))
msg("Beta row used: ", beta_row)

## alpha
load(opt$alpha_rda)
A <- as.matrix(residuals_qned_counts_alpha_objs[[1]])
alpha_row <- if (opt$alpha_trait %in% c("alpha__PD_q2")) 5L else 5L
alpha_v <- setNames(as.numeric(A[alpha_row, ]), norm_id(colnames(A)))
msg("Alpha row used: ", alpha_row)

## guild
load(opt$clusters_rda)
C <- residuals_qned_counts_clusters_objs[[1]]
rExact <- which(rownames(C) == opt$guild_trait)
rPick  <- if (length(rExact)) rExact[1] else which(startsWith(rownames(C), opt$guild_trait))[1]
if (!length(rPick) || is.na(rPick)) stop("Could not locate guild row starting with: ", opt$guild_trait)
guild_v <- setNames(as.numeric(C[rPick, ]), norm_id(colnames(C)))
msg("Guild row picked: ", rownames(C)[rPick])

## outcomes list (now includes BMI)
OUTCOMES_ALL <- c("genus_other","beta","alpha","guild","glucose","bmi")

parse_outcomes <- function(x, trait) {
  x <- tolower(trimws(as.character(x %||% "all")))
  if (!nzchar(x) || identical(x, "all")) {
    v <- OUTCOMES_ALL
  } else {
    v <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    v <- trimws(v)
    v[v %in% c("prevotella","other","other_genus","genus-other")] <- "genus_other"
    v[v %in% c("pd_p1","beta1","beta_p1")]                  <- "beta"
    v[v %in% c("pd_q2","alphaPD2","alpha_q2")]                  <- "alpha"
    v[v %in% c("cluster3","guild3")]                               <- "guild"
    v[v %in% c("wt","bodyBMI","bw","BMI")]                           <- "bmi"
    v[v %in% c("glucose","Glucose","gluc")]                           <- "glucose"
    v <- unique(intersect(v, OUTCOMES_ALL))
  }
  if (trait != "Taxa") {
    v <- setdiff(v, "genus_other")  # remove meaningless outcome outside Taxa runs
  }
  v
}

OUTCOMES <- parse_outcomes(opt$outcomes, opt$trait)
msg("Testing outcomes: ", paste(OUTCOMES, collapse = ", "))

## --- buffers for MVMR summary data (per outcome Y) ---
MVMR_BUF <- list()   # keys: outcome; values: data.table rows per SNP
OUTCOMES_MVMR <- setdiff(OUTCOMES, c("genus_other"))  # MVMR targets (no X->X)


## phenotypes (glucose and/or BMI) from same df
load(opt$phen_rda)
if (!exists("df") || !("samples" %in% names(df)))
  stop("--phen_rda must expose data.frame 'df' with a 'samples' column")

need_phen <- c(if ("glucose" %in% OUTCOMES) opt$phen_glucose_col,
               if ("bmi"  %in% OUTCOMES) opt$phen_bmi_col)
if (length(need_phen)) {
  miss_p <- need_phen[!need_phen %in% names(df)]
  if (length(miss_p))
    stop("--phen_rda$df is missing required column(s): ", paste(miss_p, collapse=", "))
}

glucose_v <- if ("glucose" %in% OUTCOMES)
  setNames(as.numeric(df[[opt$phen_glucose_col]]), norm_id(df$samples)) else NULL

BMI_v  <- if ("bmi" %in% OUTCOMES)
  setNames(as.numeric(df[[opt$phen_bmi_col]]),  norm_id(df$samples)) else NULL

## ---------------- SRMs & LOCO ----------------
Dam      <- readRDS(opt$dam_rds)
Cage     <- readRDS(opt$cage_rds)
loco_map <- readRDS(opt$loco_rds)   # chr -> prefix
msg("LOCO prefixes available: ", paste(names(loco_map), collapse=", "))

read_grm_bin <- function(prefix) {
  idf <- paste0(prefix, ".grm.id")
  bbf <- paste0(prefix, ".grm.bin")
  if (!file.exists(idf) || !file.exists(bbf))
    stop("Missing GRM files for prefix: ", prefix)
  ids <- fread(idf, header=FALSE)
  n <- nrow(ids); ntri <- n*(n+1)/2
  con <- file(bbf, "rb"); on.exit(close(con), add=TRUE)
  vals <- readBin(con, what="numeric", n=ntri, size=4)
  M <- matrix(0, n, n)
  k <- 1L
  for (j in 1:n) for (i in 1:j) { M[j,i] <- vals[k]; M[i,j] <- vals[k]; k <- k+1L }
  rn <- as.character(ids$V2)
  rownames(M) <- rn; colnames(M) <- rn
  M
}
rename_to_labels <- function(M_scan) {
  labs <- scan2lab[rownames(M_scan)]
  keep <- !is.na(labs)
  M <- M_scan[keep, keep, drop=FALSE]
  labs <- labs[keep]
  rownames(M) <- labs; colnames(M) <- labs
  M
}
.GRMs <- new.env(parent=emptyenv())
get_loco <- function(chr_str) {
  key <- sub("^chr","", as.character(chr_str))
  if (!(key %in% names(loco_map))) stop("LOCO map has no entry for chr=", key)
  if (!exists(key, envir=.GRMs)) {
    Ms <- read_grm_bin(loco_map[[key]])
    assign(key, rename_to_labels(Ms), envir=.GRMs)
  }
  get(key, envir=.GRMs)
}

## ---------------- Build harmonized DF in LABEL domain ----------------
DF <- as.data.table(genotypes)     # columns: sample (label) + gene cols
setnames(DF, "sample", "id")
DF[, id := norm_id(id)]

add_col <- function(dt, vec, nm) {
  if (is.null(vec)) {
    set(dt, j=nm, value = NA_real_)
  } else {
    set(dt, j=nm, value = as.numeric(vec[dt$id]))
  }
}

# Always mediator; outcomes added conditionally
add_col(DF, mediator_v, "mediator")
if ("genus_other" %in% OUTCOMES) add_col(DF, genus_other_v, "genus_other")
if ("beta"        %in% OUTCOMES) add_col(DF, beta_v,        "beta")
if ("alpha"        %in% OUTCOMES) add_col(DF, alpha_v,        "alpha")
if ("guild"      %in% OUTCOMES) add_col(DF, guild_v,      "guild")
if ("glucose"     %in% OUTCOMES) add_col(DF, glucose_v,     "glucose")
if ("bmi"      %in% OUTCOMES) add_col(DF, BMI_v,      "bmi")

core_cols <- c("mediator", OUTCOMES)
msg("Base DF rows (genotypes backbone): ", nrow(DF),
    " | NA counts core: ",
    paste(sprintf("%s=%d", core_cols, sapply(core_cols, function(x) sum(is.na(DF[[x]])))), collapse=" | "))
msg("Dam rows present on backbone: ",
    sum(DF$id %in% rownames(Dam)), " / ", nrow(DF),
    " | Cage rows present: ",
    sum(DF$id %in% rownames(Cage)), " / ", nrow(DF))

## ---------------- helpers ----------------
coef_from_nm <- function(nm, term) {
  if (!is.null(nm$fixef) && is.data.frame(nm$fixef)) {
    if (!term %in% rownames(nm$fixef)) return(c(NA_real_, NA_real_, NA_real_))
    est <- as.numeric(nm$fixef[term, "Est"])
    se  <- as.numeric(nm$fixef[term, "SE"])
    p   <- if ("pval" %in% colnames(nm$fixef)) {
      as.numeric(nm$fixef[term, "pval"])
    } else if (is.finite(se) && se > 0) {
      2 * pnorm(abs(est / se), lower.tail = FALSE)
    } else NA_real_
    return(c(est, se, p))
  }
  b <- NULL; V <- NULL
  if (!is.null(nm$beta))     { b <- nm$beta }
  if (!is.null(nm$varbeta))  { V <- nm$varbeta }
  if (!is.null(nm$betaCov))  { V <- nm$betaCov }
  if (!is.null(nm$vcov.beta)){ V <- nm$vcov.beta }
  if (is.null(b) || is.null(V)) return(c(NA_real_, NA_real_, NA_real_))
  se <- sqrt(diag(as.matrix(V)))
  idx <- which(names(b) == term)
  if (!length(idx)) return(c(NA_real_, NA_real_, NA_real_))
  beta <- as.numeric(b[idx]); s <- as.numeric(se[idx])
  p <- if (is.finite(s) && s > 0) 2 * pnorm(abs(beta / s), lower.tail = FALSE) else NA_real_
  c(beta, s, p)
}

classify_conclusion <- function(has_med, has_dir, has_assoc, ind_est, dir_est) {
  if (has_med && !has_dir) return("Full Mediation")
  if (has_med &&  has_dir)
    return(ifelse((ind_est < 0 & dir_est > 0) | (ind_est > 0 & dir_est < 0),
                  "Partial Mediation - Reverse", "Partial Mediation"))
  if (!has_med && has_dir && has_assoc)
    return(ifelse(dir_est > 0, "Pleiotropy - Positive", "Pleiotropy - Negative"))
  if (!has_med && has_dir && !has_assoc)
    return(ifelse(dir_est > 0, "Only Direct Effect - Positive", "Only Direct Effect - Negative"))
  "No Effect"
}

safe_r2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) < 3 || length(unique(yhat[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  cor(y[ok], yhat[ok])^2
}

## ---------------- main loop ----------------
ACC <- list(); k <- 0L

gene_names <- intersect(gene_names, setdiff(names(DF), c("id", core_cols)))
msg("Usable gene columns in DF for this fold: ", length(gene_names))

for (gname in gene_names) {
  row_corr <- corr_fold[corr_fold$Gene_name == gname, ]
  if (!nrow(row_corr)) next
  chr_str <- as.character(row_corr$chr[1])
  snp_tag <- as.character(row_corr$gene[1])
  varid   <- as.character(row_corr$variant.id[1])
  pos     <- as.integer(row_corr$pos[1])
  
  Kin_full <- tryCatch(get_loco(chr_str), error=function(e) { msg("SKIP gene ", gname, " LOCO: ", e$message); NULL })
  if (is.null(Kin_full)) next
  
  gg <- DF[[gname]]
  keepA <- DF$id[!is.na(DF$mediator) & !is.na(gg)]
  keepA <- Reduce(intersect, list(keepA, rownames(Dam), rownames(Cage), rownames(Kin_full)))
  n_alpha <- length(keepA)
  
  msg(sprintf("Gene=%s | chr=%s | snp=%s | pos=%d | varid=%s | n_alpha=%d (pre-outcome)",
              gname, chr_str, snp_tag, pos, varid, n_alpha))
  if (n_alpha < 20) next
  
  KinA  <- Kin_full[keepA, keepA, drop=FALSE]
  DamA  <- Dam [keepA, keepA, drop=FALSE]
  CageA <- Cage[keepA, keepA, drop=FALSE]
  
  ## --- alpha path: mediator ~ G + (Kin + Dam + Cage) ---
  Yk_A <- DF$mediator[match(keepA, DF$id)]
  Gk_A <- DF[[gname]][match(keepA, DF$id)]
  xA   <- data.frame(Y = Yk_A, G = Gk_A, row.names = keepA, check.names = FALSE)
  
  nmA <- tryCatch(
    GENESIS::fitNullModel(
      xA,
      outcome = "Y",
      covars  = c("G"),
      cov.mat = list(Kin = KinA, Dam = DamA, Cage = CageA),
      family  = "gaussian",
      verbose = FALSE
    ),
    error = function(e) { msg("  fitNullModel(alpha) failed for ", gname, ": ", e$message); NULL }
  )
  if (is.null(nmA)) next
  ca <- coef_from_nm(nmA, "G")
  alpha    <- ca[1]; alpha_se <- ca[2]; p_alpha <- ca[3]
  
  ## IV diagnostics
  z_alpha  <- if (is.finite(alpha_se) && alpha_se > 0) alpha / alpha_se else NA_real_
  F_alpha  <- if (is.finite(z_alpha)) z_alpha^2 else NA_real_
  R2_alpha <- safe_r2(Yk_A, alpha * Gk_A)
  
  ## --- alpha for the other genus (needed for two-exposure MVMR) ---
  alpha_other <- alpha_other_se <- F_alpha_other <- R2_alpha_other <- NA_real_
  if (opt$trait == "Taxa" && !is.null(genus_other_v)) {
    Yk_A2 <- genus_other_v[keepA]
    xA2   <- data.frame(Y = Yk_A2, G = Gk_A, row.names = keepA, check.names = FALSE)
    nmA2 <- tryCatch(
      GENESIS::fitNullModel(
        xA2, outcome = "Y", covars = c("G"),
        cov.mat = list(Kin = KinA, Dam = DamA, Cage = CageA),
        family = "gaussian", verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(nmA2)) {
      ca2 <- coef_from_nm(nmA2, "G")
      alpha_other    <- ca2[1]; alpha_other_se <- ca2[2]
      z2             <- if (is.finite(alpha_other_se) && alpha_other_se > 0) alpha_other/alpha_other_se else NA_real_
      F_alpha_other  <- if (is.finite(z2)) z2^2 else NA_real_
      R2_alpha_other <- safe_r2(Yk_A2, alpha_other * Gk_A)
    }
  }
  
  for (yn in OUTCOMES) {
    keepB <- keepA[!is.na(DF[[yn]][match(keepA, DF$id)])]
    nB <- length(keepB); if (nB < 20) { msg("  outcome=", yn, " | n=", nB, " <20 -> skip"); next }
    
    KinB  <- Kin_full[keepB, keepB, drop=FALSE]
    DamB  <- Dam [keepB, keepB, drop=FALSE]
    CageB <- Cage[keepB, keepB, drop=FALSE]
    
    Mk <- DF$mediator[match(keepB, DF$id)]
    Gk <- DF[[gname]][match(keepB, DF$id)]
    Yk <- DF[[yn]][match(keepB, DF$id)]
    
    ## total effect c via Y ~ G (same SRMs/GRM)
    xC <- data.frame(Y = Yk, G = Gk, row.names = keepB, check.names = FALSE)
    nmC <- tryCatch(
      GENESIS::fitNullModel(
        xC,
        outcome = "Y",
        covars  = c("G"),
        cov.mat = list(Kin = KinB, Dam = DamB, Cage = CageB),
        family  = "gaussian",
        verbose = FALSE
      ),
      error = function(e) { msg("  fitNullModel(total c) failed for outcome ", yn, " gene ", gname, ": ", e$message); NULL }
    )
    if (is.null(nmC)) next
    cc      <- coef_from_nm(nmC, "G")
    c_total <- cc[1]; c_total_se <- cc[2]; p_total <- cc[3]
    ## --- collect summary stats for MVMR ---
    if (opt$trait == "Taxa" && yn %in% OUTCOMES_MVMR) {
      if (is.null(MVMR_BUF[[yn]])) MVMR_BUF[[yn]] <- data.table()
      MVMR_BUF[[yn]] <- rbind(
        MVMR_BUF[[yn]],
        data.table(
          snp = snp_tag,
          bx_anchor = alpha,          bxse_anchor = alpha_se,          F_anchor = F_alpha,          R2_anchor = R2_alpha,
          bx_other  = alpha_other,    bxse_other  = alpha_other_se,    F_other  = F_alpha_other,    R2_other  = R2_alpha_other,
          by        = c_total,        byse        = c_total_se
        ),
        fill = TRUE
      )
    }
    z_total <- if (is.finite(c_total_se) && c_total_se > 0) c_total / c_total_se else NA_real_
    F_total <- if (is.finite(z_total)) z_total^2 else NA_real_
    R2_total <- safe_r2(Yk, c_total * Gk)
    steiger_ok <- is.finite(R2_alpha) && is.finite(R2_total) && (R2_alpha > R2_total)
    
    ## outcome model with mediator: Y ~ M + G
    xB <- data.frame(Y = Yk, M = Mk, G = Gk, row.names = keepB, check.names = FALSE)
    nmB <- tryCatch(
      GENESIS::fitNullModel(
        xB,
        outcome = "Y",
        covars  = c("M","G"),
        cov.mat = list(Kin = KinB, Dam = DamB, Cage = CageB),
        family  = "gaussian",
        verbose = FALSE
      ),
      error = function(e) { msg("  fitNullModel(b,c') failed for outcome ", yn, " gene ", gname, ": ", e$message); NULL }
    )
    if (is.null(nmB)) next
    
    b  <- coef_from_nm(nmB, "M")
    c_ <- coef_from_nm(nmB, "G")
    
    beta      <- b[1];  beta_se   <- b[2];  p_beta   <- b[3]
    direct    <- c_[1]; direct_se <- c_[2]; p_direct <- c_[3]
    
    indirect    <- alpha * beta
    indirect_se <- sqrt((alpha^2)*(beta_se^2) + (beta^2)*(alpha_se^2))
    p_indirect  <- if (is.finite(indirect_se) && indirect_se > 0) 2*pnorm(abs(indirect/indirect_se), lower.tail=FALSE) else NA_real_
    total       <- direct + indirect
    prop_med    <- if (!is.na(total) && total != 0) 100*indirect/total else NA_real_
    
    has_mediation     <- !is.na(p_indirect) && p_indirect < 0.05
    has_direct_effect <- !is.na(p_direct)   && p_direct   < 0.05
    has_assoc         <- !is.na(p_alpha)    && p_alpha    < 0.05
    conclusion <- classify_conclusion(has_mediation, has_direct_effect, has_assoc, indirect, direct)
    
    ## Wald ratio
    iv_theta <- if (is.finite(alpha) && alpha != 0) c_total / alpha else NA_real_
    iv_theta_se <- if (is.finite(iv_theta) && is.finite(c_total_se) && is.finite(alpha_se) && alpha != 0 && is.finite(c_total) && c_total != 0) {
      abs(iv_theta) * sqrt( (c_total_se / c_total)^2 + (alpha_se / alpha)^2 )
    } else NA_real_
    p_iv_theta <- if (is.finite(iv_theta_se) && iv_theta_se > 0) 2*pnorm(abs(iv_theta/iv_theta_se), lower.tail=FALSE) else NA_real_
    
    k <- k + 1L
    ACC[[k]] <- data.table(
      gene_name   = gname,
      snp         = snp_tag,
      chr         = as.character(chr_str),
      pos         = pos,
      variant_id  = varid,
      outcome     = yn,
      mediator    = opt$genus,
      n_alpha     = n_alpha,
      n_outcome   = nB,
      alpha       = alpha,   alpha_se   = alpha_se,   p_alpha   = p_alpha,
      F_alpha     = F_alpha,
      R2_alpha    = R2_alpha,
      beta        = beta,    beta_se    = beta_se,    p_beta    = p_beta,
      direct      = direct,  direct_se  = direct_se,  p_direct  = p_direct,
      indirect    = indirect,indirect_se= indirect_se,p_indirect= p_indirect,
      total       = total,   prop_med   = prop_med,
      c_total     = c_total, c_total_se = c_total_se, p_total   = p_total,
      F_total     = F_total, R2_total   = R2_total,  steiger_ok = steiger_ok,
      iv_theta    = iv_theta, iv_theta_se = iv_theta_se, p_iv_theta = p_iv_theta,
      conclusion  = conclusion
    )
    
    msg(sprintf("  outcome=%-10s | n=%-3d | alpha=% .3g (p=% .2g) F=%.1f R2=%.3f | "
                , yn, nB, alpha, p_alpha, F_alpha %||% NaN, R2_alpha %||% NaN))
    msg(sprintf("                | c=% .3g (p=% .2g) F=%.1f R2=%.3f | c'=% .3g (p=% .2g) | "
                , c_total, p_total, F_total %||% NaN, R2_total %||% NaN, direct, p_direct))
    msg(sprintf("                | beta=% .3g (p=% .2g) | ind=% .3g (p=% .2g) | Steiger=%s | Wald θ=% .3g (p=% .2g) | %s",
                beta, p_beta, indirect, p_indirect,
                if (isTRUE(steiger_ok)) "OK" else "FAIL",
                iv_theta %||% NaN, p_iv_theta %||% NaN, conclusion))
  }
}

if (!length(ACC)) {
  msg("No rows produced. Common causes: too few subjects after SRMs/LOCO, or null model failures.")
  quit(save="no", status=2)
}

RES <- rbindlist(ACC, fill=TRUE)

# Multiple testing control by outcome
RES[, q_alpha    := p.adjust(p_alpha,    method="BH"), by=.(outcome)]
RES[, q_beta     := p.adjust(p_beta,     method="BH"), by=.(outcome)]
RES[, q_direct   := p.adjust(p_direct,   method="BH"), by=.(outcome)]
RES[, q_indirect := p.adjust(p_indirect, method="BH"), by=.(outcome)]
RES[, q_total    := p.adjust(p_total,    method="BH"), by=.(outcome)]
RES[, q_iv_theta := p.adjust(p_iv_theta, method="BH"), by=.(outcome)]

use_q   <- tolower(opt$sig_metric) == "q"
pthr    <- if (!is.null(opt$p_level))  opt$p_level  else 0.05
qdash   <- if (!is.null(opt$fdr_level)) opt$fdr_level else 0.10

RES[, `:=`(
  sig_alpha    = if (use_q) q_alpha    < qdash else p_alpha    < pthr,
  sig_beta     = if (use_q) q_beta     < qdash else p_beta     < pthr,
  sig_direct   = if (use_q) q_direct   < qdash else p_direct   < pthr,
  sig_indirect = if (use_q) q_indirect < qdash else p_indirect < pthr,
  sig_total    = if (use_q) q_total    < qdash else p_total    < pthr,
  sig_iv_theta = if (use_q) q_iv_theta < qdash else p_iv_theta < pthr
)]

# Narrative conclusion with chosen metric
classify <- function(has_med, has_dir, has_assoc, ind_est, dir_est) {
  if (has_med && !has_dir) return("Full Mediation")
  if (has_med &&  has_dir)
    return(ifelse((ind_est < 0 & dir_est > 0) | (ind_est > 0 & dir_est < 0),
                  "Partial Mediation - Reverse", "Partial Mediation"))
  if (!has_med && has_dir && has_assoc)
    return(ifelse(dir_est > 0, "Pleiotropy - Positive", "Pleiotropy - Negative"))
  if (!has_med && has_dir && !has_assoc)
    return(ifelse(dir_est > 0, "Only Direct Effect - Positive", "Only Direct Effect - Negative"))
  "No Effect"
}
RES[, conclusion := mapply(classify, sig_indirect, sig_direct, sig_alpha, indirect, direct)]

# IV credibility flags
RES[, strong_iv := is.finite(F_alpha) & (F_alpha >= 10)]
RES[, iv_ok := strong_iv & isTRUE(steiger_ok) & sig_alpha & !sig_direct]

# Decision breadcrumb
RES[, sig_rule := if (use_q) sprintf("FDR<%.2f", qdash) else sprintf("p<%.3f", pthr)]

## =================== MVMR (two-exposure) per outcome ===================
MVMR_ROWS <- list()
if (isTRUE(opt$mr_enable) && opt$trait == "Taxa" && length(MVMR_BUF)) {
  anchor <- opt$genus
  other  <- opt$other_genus
  
  for (yn in names(MVMR_BUF)) {
    # --- gate rows that have everything MVMR needs ---
    need <- c("bx_anchor","bxse_anchor","bx_other","bxse_other","by","byse")
    DTm0 <- MVMR_BUF[[yn]][complete.cases(MVMR_BUF[[yn]][, ..need])]
    
    # Log strength distribution
    qq <- quantile(pmax(DTm0$F_anchor, DTm0$F_other), c(.1,.25,.5,.75,.9,.95), na.rm=TRUE)
    msg(sprintf("[MVMR:%s] max(F_anchor,F_other) quantiles: %s", yn, paste(round(qq,2), collapse=", ")))
    msg(sprintf("[MVMR:%s] raw=%d | complete=%d", yn, nrow(MVMR_BUF[[yn]]), nrow(DTm0)))
    
    thr <- ifelse(is.finite(opt$mr_fstat_min), opt$mr_fstat_min, 10)
    mr_min_snps <- if (!is.null(opt$mr_min_snps)) as.integer(opt$mr_min_snps) else 3L
    
    KEEP10 <- with(DTm0, (is.finite(F_anchor) & F_anchor >= thr) | (is.finite(F_other) & F_other >= thr))
    DTm <- DTm0[KEEP10]
    using_fallback <- FALSE
    if (nrow(DTm) < mr_min_snps) {
      KEEP5 <- with(DTm0, (is.finite(F_anchor) & F_anchor >= 5) | (is.finite(F_other) & F_other >= 5))
      DTm <- DTm0[KEEP5]
      if (nrow(DTm) < mr_min_snps && nrow(DTm0) >= mr_min_snps) {
        DTm <- DTm0[order(pmax(F_anchor, F_other), decreasing=TRUE)][1:mr_min_snps]
      }
      using_fallback <- TRUE
    }
    msg(sprintf("[MVMR:%s] passF>=%.1f=%d | using_fallback=%s", yn, thr, sum(KEEP10, na.rm=TRUE), using_fallback))
    if (nrow(DTm) < mr_min_snps) {
      msg(sprintf("[MVMR:%s] insufficient SNPs after gating (n=%d < %d) -> skip", yn, nrow(DTm), mr_min_snps))
      next
    }
    
    # --- build MVMR input (current API) ---
    r_input <- tryCatch(
      MVMR::format_mvmr(
        BXGs  = as.matrix(DTm[, .(bx_anchor, bx_other)]),
        BYG   = DTm$by,
        seBXGs= as.matrix(DTm[, .(bxse_anchor, bxse_other)]),
        seBYG = DTm$byse,
        RSID  = DTm$snp
      ),
      error = function(e) { msg(sprintf("[MVMR:%s] format_mvmr error: %s", yn, e$message)); NULL }
    )
    if (is.null(r_input)) next
    
    # helpers to extract values robustly across return types
    get_col <- function(obj, col) {
      if (is.null(obj)) return(NA_real_)
      if (inherits(obj, "error")) return(NA_real_)
      if (is.data.frame(obj) || is.matrix(obj)) {
        cn <- colnames(obj); try_names <- c(col, gsub(" ", ".", col), gsub("[ .]", "", col))
        hit <- try_names[try_names %in% cn]
        if (length(hit)) return(as.numeric(obj[, hit[1]]))
        return(NA_real_)
      }
      if (is.list(obj) && !is.null(obj[[col]])) return(as.numeric(obj[[col]]))
      if (is.numeric(obj)) {
        nm <- names(obj)
        if (!is.null(nm) && col %in% nm) return(as.numeric(obj[[col]]))
        return(as.numeric(obj))
      }
      NA_real_
    }
    
    # --- run strength / IVW / Q with gencov=0 to silence warnings ---
    st  <- tryCatch(MVMR::strength_mvmr(r_input, gencov=0),   error=function(e) { msg(sprintf("[MVMR:%s] strength_mvmr error: %s", yn, e$message)); e })
    ivw <- tryCatch(MVMR::ivw_mvmr(r_input, gencov=0),        error=function(e) { msg(sprintf("[MVMR:%s] ivw_mvmr error: %s", yn, e$message)); e })
    pq  <- tryCatch(MVMR::pleiotropy_mvmr(r_input, gencov=0), error=function(e) { msg(sprintf("[MVMR:%s] pleiotropy_mvmr error: %s", yn, e$message)); e })
    
    # strength_mvmr often returns a 1x2 matrix with rowname like "F-statistic"
    F1 <- F2 <- NA_real_
    if (is.matrix(st) || is.data.frame(st)) {
      # try to use column names "exposure1","exposure2"; otherwise take the first two entries
      if (!is.null(colnames(st)) && all(c("exposure1","exposure2") %in% colnames(st))) {
        F1 <- as.numeric(st[1, "exposure1"])
        F2 <- as.numeric(st[1, "exposure2"])
      } else {
        vals <- as.numeric(st)[1:min(2, length(st))]
        if (length(vals) >= 1) F1 <- vals[1]
        if (length(vals) >= 2) F2 <- vals[2]
      }
    } else if (is.numeric(st)) {
      vals <- as.numeric(st)[1:min(2, length(st))]
      if (length(vals) >= 1) F1 <- vals[1]
      if (length(vals) >= 2) F2 <- vals[2]
    }
    
    # ivw_mvmr returns a 2x4 table (rows=exposures; cols ~ "Estimate","Std. Error","t value","Pr(>|t|)")
    beta_vec <- se_vec <- p_vec <- rep(NA_real_, 2L)
    if (is.data.frame(ivw) || is.matrix(ivw)) {
      # row order is exposure1, exposure2
      # columns may be "Std. Error" (printed) -> in object it's often "Std..Error" or similar; use get_col
      ests <- get_col(ivw, "Estimate")
      ses  <- get_col(ivw, "Std. Error"); if (all(is.na(ses))) ses <- get_col(ivw, "Std.Error")
      pval <- get_col(ivw, "Pr(>|t|)");  if (all(is.na(pval))) pval <- get_col(ivw, "p.value")
      # If columns come back as vectors of length 2, use them
      if (length(ests) >= 2) beta_vec <- as.numeric(ests[1:2])
      if (length(ses)  >= 2) se_vec   <- as.numeric(ses[1:2])
      if (length(pval) >= 2) p_vec    <- as.numeric(pval[1:2])
      # compute p’s if missing and SE available
      for (j in 1:2) if (is.na(p_vec[j]) && is.finite(beta_vec[j]) && is.finite(se_vec[j]) && se_vec[j] > 0)
        p_vec[j] <- 2*pnorm(abs(beta_vec[j]/se_vec[j]), lower.tail=FALSE)
    } else if (is.numeric(ivw)) {
      # very old versions may return a named numeric vector — try best effort
      beta_vec[1] <- ivw[1]; if (length(ivw) >= 2) beta_vec[2] <- ivw[2]
    }
    
    # Q-statistic: sometimes returns a small object with elements; sometimes just prints.
    Qstat <- Qp <- NA_real_
    if (is.list(pq)) {
      Qstat <- as.numeric(pq$Qstat %||% pq$Q %||% NA_real_)
      Qp    <- as.numeric(pq$Qstat_pval %||% pq$Qp %||% NA_real_)
    } else if (is.numeric(pq) && length(pq) >= 2) {
      Qstat <- as.numeric(pq[1]); Qp <- as.numeric(pq[2])
    }
    
    MVMR_ROWS[[length(MVMR_ROWS)+1L]] <- data.table(
      analysis      = "mvmr",
      anchor        = anchor,
      exp1          = anchor,  exp2 = other,
      outcome       = yn,
      n_snps        = nrow(DTm),
      F_cond_exp1   = F1,
      F_cond_exp2   = F2,
      beta_exp1_dir = beta_vec[1],  se_exp1_dir = se_vec[1],  p_exp1_dir = p_vec[1],
      beta_exp2_dir = beta_vec[2],  se_exp2_dir = se_vec[2],  p_exp2_dir = p_vec[2],
      Q_mvmr        = Qstat,
      Qp_mvmr       = Qp,
      weak_iv_gate  = using_fallback
    )
    
  }
}

## tag mediation rows so the meta script can tell them apart (back-compat if missing)
## tag mediation rows
if (!"analysis" %in% names(RES)) RES[, analysis := "med"]

## collapse MVMR rows (if any)
MVMR_DT <- if (length(MVMR_ROWS)) rbindlist(MVMR_ROWS, fill=TRUE) else data.table()
if (nrow(MVMR_DT)) MVMR_DT[, analysis := "mvmr"]

## ---------------- Conclusion labels ----------------
# For mediation/univariable rows (per SNP × outcome)
call_conclusion_med <- function(sig_iv, iv_ok, sig_total, sig_direct) {
  if (isTRUE(iv_ok) && isTRUE(sig_iv))      return("True causal")
  if (isTRUE(sig_total) && (!isTRUE(iv_ok) || isTRUE(sig_direct))) return("Correlated pleiotropy")
  return("Not causal")
}
RES[, conclusion_simple :=
      mapply(call_conclusion_med, sig_iv_theta, iv_ok, sig_total, sig_direct)]

# For MVMR (per outcome): anchor vs other_genus -> Yn
call_conclusion_mvmr <- function(p1, F1, p2, F2, Qp) {
  strong1 <- is.finite(F1) && F1 >= 10
  strong2 <- is.finite(F2) && F2 >= 10
  sig1 <- is.finite(p1) && p1 < 0.05
  sig2 <- is.finite(p2) && p2 < 0.05
  pleio <- is.finite(Qp) && Qp < 0.05
  if (sig1 && sig2 && strong1 && strong2) return("Bidirectional")
  if (sig1 && strong1 && !sig2)           return("True causal")
  if (sig2 && strong2 && !sig1)           return("Not causal")     # favors other_genus, not anchor
  if (pleio || (!strong1 && sig1) || (!strong2 && sig2)) return("Correlated pleiotropy")
  return("Not causal")
}
if (nrow(MVMR_DT)) {
  MVMR_DT[, conclusion_simple :=
            mapply(call_conclusion_mvmr, p_exp1_dir, F_cond_exp1, p_exp2_dir, F_cond_exp2, Qp_mvmr)]
}

## Unify + write once
ALL <- rbindlist(list(RES, MVMR_DT), fill=TRUE)
fwrite(ALL, file=opt$out, sep="\t")
msg("[med:genesis_loco + mvmr] wrote: ", opt$out, " | rows=", nrow(ALL))
quit(save="no")

