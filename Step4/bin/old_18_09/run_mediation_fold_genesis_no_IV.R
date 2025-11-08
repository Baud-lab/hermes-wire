#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GENESIS)
})

## ---------------- CLI ----------------
opt_list <- list(
  make_option("--selected_genes_rdata", type="character"),
  make_option("--snps_tsv",             type="character"),
  make_option("--resid_rda",            type="character"),
  make_option("--beta_rda",             type="character"),
  make_option("--clusters_rda",         type="character"),
  make_option("--phen_rda",             type="character"),
  make_option("--phen_glucose_col",     type="character", default="glucose_at_dissection"),
  make_option("--dam_rds",              type="character"),
  make_option("--cage_rds",             type="character"),
  make_option("--loco_rds",             type="character"),     # RDS: chr -> GRM prefix
  make_option("--scan_lookup",          type="character"),     # NEW: TSV with scan_id, sample_label
  make_option("--genus",                type="character"),
  make_option("--other_genus",          type="character"),
  make_option("--trait_prefix",         type="character", default="g__"),
  make_option("--beta_trait",           type="character", default="beta__TD_PCoA1"),
  make_option("--guild3_trait",         type="character", default="cluster__3"),
  make_option(
    "--outcomes",
    type    = "character",
    default = "genus_other",
    help    = "Comma-separated subset of outcomes to test (genus_other,beta,guild3,glucose) or 'all'."
  ),
  make_option("--sig_metric", type="character", default="q",
              help="Use 'p' or 'q' to call significance for alpha/beta/direct/indirect"),
  make_option("--fdr_level",  type="double",    default=0.10,
              help="FDR threshold when --sig_metric=q"),
  make_option("--p_level",    type="double",    default=0.05,
              help="Nominal p-value threshold when --sig_metric=p"),
  make_option("--out",                  type="character", default="mediation_results.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

msg     <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
norm_id <- function(x) sub("_.*$", "", as.character(x))

## ---------------- sanity/files ----------------
required <- c("selected_genes_rdata","snps_tsv","resid_rda","beta_rda","clusters_rda",
              "phen_rda","dam_rds","cage_rds","loco_rds","scan_lookup")
miss <- required[!file.exists(unlist(opt[required]))]
if (length(miss)) stop("Missing required file(s): ", paste(miss, collapse=", "))

## ---------------- lookup (to harmonize IDs) ----------------
LK <- fread(opt$scan_lookup)  # must have: scan_id, sample_label
if (!all(c("scan_id","sample_label") %in% names(LK)))
  stop("--scan_lookup must have columns: scan_id, sample_label")
LK[, scan_id := as.character(scan_id)]
LK[, sample_label_trim := norm_id(sample_label)]
lab2scan <- setNames(LK$scan_id,          LK$sample_label_trim)
scan2lab <- setNames(LK$sample_label_trim, LK$scan_id)

## ---------------- load selected genes ----------------
load(opt$selected_genes_rdata)   # needs: genotypes, correlations
if (!exists("genotypes"))   stop("selected_genes_rdata lacks object 'genotypes'")
if (!exists("correlations")) stop("selected_genes_rdata lacks object 'correlations'")

stopifnot("sample" %in% names(genotypes))
genotypes$sample <- norm_id(genotypes$sample)
genotypes <- genotypes[!duplicated(genotypes$sample), , drop = FALSE]  # rows, not columns

# SNPs to test for this fold
SNPs <- fread(opt$snps_tsv)
if (!("snp" %in% names(SNPs))) stop("--snps_tsv must contain a column named 'snp' (chr:pos)")
snps_fold <- unique(SNPs$snp)

## correlations: expected columns
need_cols <- c("Gene_name","gene","chr","pos","variant.id")
if (!all(need_cols %in% names(correlations)))
  stop("correlations missing required columns: ", paste(setdiff(need_cols, names(correlations)), collapse=", "))
corr_fold <- correlations[correlations$gene %in% snps_fold, ]
if (nrow(corr_fold) == 0L) stop("No rows in 'correlations' match the fold SNPs")
gene_names <- intersect(corr_fold$Gene_name, setdiff(names(genotypes), "sample"))
if (!length(gene_names)) stop("No genotype columns in 'genotypes' match 'Gene_name' for this fold")
msg("Fold genes (after matching genotypes): ", length(gene_names))

## ---------------- residuals / beta / clusters / phenotype ----------------
load(opt$resid_rda)
if (!exists("residuals_qned_counts_objs"))
  stop("Expected 'residuals_qned_counts_objs' in --resid_rda")
genus_layer_idx <- length(residuals_qned_counts_objs) - 1L
R_genus <- t(residuals_qned_counts_objs[[genus_layer_idx]])
rownames(R_genus) <- norm_id(rownames(R_genus))

pick_col <- function(M, key) {
  cn <- colnames(M)
  exact <- which(cn == key)
  if (length(exact)) return(cn[exact[1]])
  st <- which(startsWith(cn, key))
  if (length(st)) return(cn[st[1]])
  NA_character_
}
med_key   <- paste0(opt$trait_prefix, opt$genus)
other_key <- paste0(opt$trait_prefix, opt$other_genus)
med_col   <- pick_col(R_genus, med_key)
other_col <- pick_col(R_genus, other_key)
if (is.na(med_col))   stop("Mediator trait not found in residuals: ", med_key)
if (is.na(other_col)) stop("Other-genus trait not found in residuals: ", other_key)

mediator_v    <- setNames(as.numeric(R_genus[, med_col]),   rownames(R_genus))
genus_other_v <- setNames(as.numeric(R_genus[, other_col]), rownames(R_genus))
msg("Mediator col='", med_col, "' | Other-genus col='", other_col, "'")

load(opt$beta_rda)
if (!exists("residuals_qned_counts_beta_objs"))
  stop("Expected 'residuals_qned_counts_beta_objs' in --beta_rda")
B <- as.matrix(residuals_qned_counts_beta_objs[[1]])
beta_row <- if (opt$beta_trait %in% c("beta__TD_PCoA1","Axis.1")) 5L else if (opt$beta_trait %in% c("beta__TD_PCoA2","Axis.2")) 6L else 5L
beta_v <- setNames(as.numeric(B[beta_row, ]), norm_id(colnames(B)))
msg("Beta row used: ", beta_row)

load(opt$clusters_rda)
if (!exists("residuals_qned_counts_clusters_objs"))
  stop("Expected 'residuals_qned_counts_clusters_objs' in --clusters_rda")
C <- residuals_qned_counts_clusters_objs[[1]]
rExact <- which(rownames(C) == opt$guild3_trait)
rPick  <- if (length(rExact)) rExact[1] else which(startsWith(rownames(C), opt$guild3_trait))[1]
guild3_v <- setNames(as.numeric(C[rPick, ]), norm_id(colnames(C)))
msg("Guild row picked: ", rownames(C)[rPick])

load(opt$phen_rda)
if (!exists("df") || !all(c("samples", opt$phen_glucose_col) %in% names(df)))
  stop("--phen_rda must expose data.frame 'df' with 'samples' and '", opt$phen_glucose_col, "'")
glucose_v <- setNames(df[[opt$phen_glucose_col]], norm_id(df$samples))

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
  rn <- as.character(ids$V2)                  # scan IDs
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
    assign(key, rename_to_labels(Ms), envir=.GRMs)  # convert to LABEL IDs
  }
  get(key, envir=.GRMs)
}

## ---------------- Build harmonized DF in LABEL domain ----------------
DF <- as.data.table(genotypes)     # columns: sample (label) + gene cols
setnames(DF, "sample", "id")
DF[, id := norm_id(id)]

add_col <- function(dt, vec, nm) set(dt, j=nm, value = as.numeric(vec[dt$id]))
add_col(DF, mediator_v,    "mediator")
add_col(DF, genus_other_v, "genus_other")
add_col(DF, beta_v,        "beta")
add_col(DF, guild3_v,      "guild3")
add_col(DF, glucose_v,     "glucose")

core_cols <- c("mediator","genus_other","beta","guild3","glucose")
msg("Base DF rows (genotypes backbone): ", nrow(DF),
    " | NA counts core: ",
    paste(sprintf("%s=%d", core_cols, sapply(core_cols, function(x) sum(is.na(DF[[x]])))), collapse=" | "))
msg("Dam rows present on backbone: ",
    sum(DF$id %in% rownames(Dam)), " / ", nrow(DF),
    " | Cage rows present: ",
    sum(DF$id %in% rownames(Cage)), " / ", nrow(DF))

## ---------------- helpers ----------------
coef_from_nm <- function(nm, term) {
  # Preferred: modern GENESIS returns a data.frame at $fixef
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
  
  # Fallbacks: older objects
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

OUTCOMES_ALL <- c("genus_other","beta","guild3","glucose")

parse_outcomes <- function(x) {
  x <- tolower(trimws(as.character(x %||% "all")))
  if (!nzchar(x) || identical(x, "all")) return(OUTCOMES_ALL)
  v <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  v <- trimws(v)
  
  # friendly synonyms
  v[v %in% c("prevotella","other","other_genus","genus-other")] <- "genus_other"
  v[v %in% c("td_pcoa1","beta1","beta_pcoa1")]                  <- "beta"
  v[v %in% c("cluster3","guild")]                               <- "guild3"
  
  unique(intersect(v, OUTCOMES_ALL))
}
`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a


## ---------------- main loop ----------------
OUTCOMES <- parse_outcomes(opt$outcomes)
if (!length(OUTCOMES)) {
  stop("No valid outcomes specified. Choose among: ",
       paste(OUTCOMES_ALL, collapse = ", "))
}
msg("Testing outcomes: ", paste(OUTCOMES, collapse = ", "))
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
  
  # LOCO GRM for this gene's chromosome (now in LABEL domain)
  Kin_full <- tryCatch(get_loco(chr_str), error=function(e) { msg("SKIP gene ", gname, " LOCO: ", e$message); NULL })
  if (is.null(Kin_full)) next
  
  # Base set in LABEL space
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
  Yk <- DF$mediator[match(keepA, DF$id)]
  Gk <- DF[[gname]][match(keepA, DF$id)]
  xA <- data.frame(Y = Yk, G = Gk, row.names = keepA, check.names = FALSE)
  
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
  alpha <- ca[1]; alpha_se <- ca[2]; p_alpha <- ca[3]
  
  for (yn in OUTCOMES) {
    keepB <- keepA[!is.na(DF[[yn]][match(keepA, DF$id)])]
    nB <- length(keepB); if (nB < 20) { msg("  outcome=", yn, " | n=", nB, " <20 -> skip"); next }
    
    KinB  <- Kin_full[keepB, keepB, drop=FALSE]
    DamB  <- Dam [keepB, keepB, drop=FALSE]
    CageB <- Cage[keepB, keepB, drop=FALSE]
    
    ## --- outcome paths: Y ~ M + G + (Kin + Dam + Cage) ---
    Mk <- DF$mediator[match(keepB, DF$id)]
    Gk <- DF[[gname]][match(keepB, DF$id)]
    Yk <- DF[[yn]][match(keepB, DF$id)]
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
      beta        = beta,    beta_se    = beta_se,    p_beta    = p_beta,
      direct      = direct,  direct_se  = direct_se,  p_direct  = p_direct,
      indirect    = indirect,indirect_se= indirect_se,p_indirect= p_indirect,
      total       = total,   prop_med   = prop_med,
      conclusion  = conclusion
    )
    
    msg(sprintf("  outcome=%-10s | n=%-3d | alpha=% .3g (p=% .2g) | beta=% .3g (p=% .2g) | c'=% .3g (p=% .2g) | ind=% .3g (p=% .2g) | total=% .3g | %%med=%.2f | %s",
                yn, nB, alpha, p_alpha, beta, p_beta, direct, p_direct, indirect, p_indirect, total, prop_med, conclusion))
  }
}

if (!length(ACC)) {
  msg("No rows produced. Common causes: too few subjects after SRMs/LOCO, or null model failures.")
  quit(save="no", status=2)
}

RES <- rbindlist(ACC, fill=TRUE)
# compute q-values (by outcome, matching your other scripts)
RES[, q_alpha    := p.adjust(p_alpha,    method="BH"), by=.(outcome)]
RES[, q_beta     := p.adjust(p_beta,     method="BH"), by=.(outcome)]
RES[, q_direct   := p.adjust(p_direct,   method="BH"), by=.(outcome)]
RES[, q_indirect := p.adjust(p_indirect, method="BH"), by=.(outcome)]

use_q   <- tolower(opt$sig_metric) == "q"
pthr    <- if (!is.null(opt$p_level))  opt$p_level  else 0.05
qdash   <- if (!is.null(opt$fdr_level)) opt$fdr_level else 0.10

RES[, `:=`(
  sig_alpha    = if (use_q) q_alpha    < qdash else p_alpha    < pthr,
  sig_beta     = if (use_q) q_beta     < qdash else p_beta     < pthr,
  sig_direct   = if (use_q) q_direct   < qdash else p_direct   < pthr,
  sig_indirect = if (use_q) q_indirect < qdash else p_indirect < pthr
)]

# re-make 'conclusion' using the chosen metric
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

RES[, conclusion :=
      mapply(classify, sig_indirect, sig_direct, sig_alpha, indirect, direct)]

# (optional) keep a breadcrumb about the decision rule used
RES[, sig_rule := if (use_q) sprintf("FDR<%.2f", qdash) else sprintf("p<%.3f", pthr)]

fwrite(RES, file=opt$out, sep="\t")
msg("[med:genesis_loco] wrote: ", opt$out, " | rows=", nrow(RES))
