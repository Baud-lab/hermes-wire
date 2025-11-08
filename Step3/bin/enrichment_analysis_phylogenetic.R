#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# PGLS-based function enrichment with phylogenetic correction (Pagel's λ)
# Thin CLI added; core analysis logic & outputs preserved.
#
# Key invariants vs your original:
#  - REPORT_MODEL == "ML" => Direction == Direction_ML (enforced).
#  - Sign-stability diagnostics (strict & available) retained.
#  - Robust flag: q<0.10, ΔAIC_used<=2, and sign_stable_avail==TRUE.
#  - Saves BOTH res_full (all columns) and res_filt (decision columns) per (var, group, func).
#  - Volcano PDF saved; per-(var,group,func) .RData saved; CSV of significant rows saved.
#
# Required objects in inputs:
#  - Heritability .RData must contain all_VCs_full (as in your original).
#  - Function matrices .RData must contain func_matrix_<FUNCTION>, e.g. func_matrix_eggNOG_OGs
#  - Taxonomy file: GTDB (Genome \t Taxonomy) or Greengenes2 (Genome \t Taxonomy \t ...).
#
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ape)
  library(nlme)
  library(ggplot2)
  library(data.table)
})

# ------------------------- Minimal CLI -------------------------
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k, d=NULL) {
  i <- which(args == paste0("--",k))
  if (!length(i)) return(d)
  val <- args[i+1]
  if (is.na(val)) d else val
}

taxonomy_path  <- arg("taxonomy")             # e.g., taxonomy_GTDB207.txt
phylo_path     <- arg("phylo")                # e.g., phylotree_GTDB207.tree
herit_path     <- arg("herit")                # e.g., .../Heritability/heritability.RData
func_mats_rda  <- arg("func_mats")            # e.g., .../func_matrices_v4.RData
func_enog_rda  <- arg("func_enog", func_mats_rda) # optional extra .RData (if separate eggNOG OGs file)
subset_id      <- arg("subset", "ALL")
method         <- arg("method", "Shallow")

REPORT_MODEL   <- arg("report", "ML")         # "ML" or "AIC"
branch_length_transform <- arg("branchlen", "none") # "none" | "grafen"
outdir         <- arg("outdir", ".")
groups_csv     <- arg("groups", "g__Prevotella,g__Bacteroides")
functions_csv  <- arg("functions", "eggNOG_OGs")
variables_csv  <- arg("variables", "Heritability") # keep original default

# Output names (preserve your original naming where possible)
pdf_out        <- arg("pdf", file.path(outdir, "enrichment_heritability_PGLS_volcano_ALL.pdf"))

# Optional post-hoc extras (off by default)
run_extras     <- tolower(arg("run_extras", "false")) %in% c("1","true","yes")
extras_bac_xlsx <- arg("extras_bac_xlsx", "")
extras_prev_xlsx <- arg("extras_prev_xlsx", "")
extras_pdf_heat <- arg("extras_pdf_heat", file.path(outdir, "heatmaps_enrichment_ALL.pdf"))
extras_pdf_qq   <- arg("extras_pdf_qq",   file.path(outdir, "qqplot_bac_prev_genomes.pdf"))

if (is.null(taxonomy_path) || is.null(phylo_path) || is.null(herit_path) || is.null(func_mats_rda)) {
  stop("Required: --taxonomy --phylo --herit --func_mats [--func_enog --subset --method --report --branchlen --outdir --groups --functions --variables --pdf]")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
GROUPS   <- strsplit(groups_csv, ",")[[1]]    |> trimws()
FUNCS    <- strsplit(functions_csv, ",")[[1]] |> trimws()
VARIABLES<- strsplit(variables_csv, ",")[[1]] |> trimws()

# ------------------------- Helpers (unchanged logic) -------------------------
read_gzipped   <- function(f) if (grepl("\\.gz$", f)) gzfile(f) else f
read_tab_file  <- function(f, header=FALSE) read.delim(read_gzipped(f), header=header, check.names=FALSE)
read_csv_file  <- function(f) read.csv(read_gzipped(f), stringsAsFactors = FALSE)

# rank-based inverse-normal transform (keeps NA)
invrank <- function(x) qnorm((rank(x, na.last="keep", ties.method='random') - 0.5) / sum(!is.na(x)))

is_constant_column <- function(col) length(unique(col)) == 1

ensure_branch_lengths <- function(tr, method = c("none", "grafen")) {
  method <- match.arg(method)
  if (is.null(tr$edge.length)) {
    if (method == "grafen") tr <- compute.brlen(tr, method = "Grafen") else tr$edge.length <- rep(1, nrow(tr$edge))
  }
  if (any(tr$edge.length <= 0)) tr$edge.length[tr$edge.length <= 0] <- 1e-8
  tr
}

make_species_tree <- function(master_tree, species_vec, species2genome, genome2species, branchlen_method = "none") {
  g <- unname(species2genome[species_vec]); g <- g[!is.na(g)]
  common <- intersect(master_tree$tip.label, g)
  if (length(common) < 3) return(NULL)
  sub <- keep.tip(master_tree, common)
  tip_species <- genome2species[sub$tip.label]
  keep_idx <- !duplicated(tip_species)
  if (sum(keep_idx) < 3) return(NULL)
  sub <- keep.tip(sub, sub$tip.label[keep_idx])
  sub$tip.label <- tip_species[keep_idx]
  ensure_branch_lengths(sub, method = branchlen_method)
}

# Extract λ on [0,1]; if unconstrained value returned, back-transform via logistic
get_fitted_lambda <- function(fit) {
  lam <- tryCatch(as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)[1]), error = function(e) NA_real_)
  if (!is.finite(lam) || lam < 0 || lam > 1) {
    lam_u <- tryCatch(as.numeric(coef(fit$modelStruct$corStruct, unconstrained = TRUE)[1]), error = function(e) NA_real_)
    if (is.finite(lam_u)) lam <- 1/(1+exp(-lam_u))
  }
  if (!is.finite(lam)) NA_real_ else min(max(lam, 0), 1)
}

# Small helper: robust sign that treats |x|<eps as 0 to avoid spurious flips
safe_sign <- function(x, eps = 1e-12) {
  out <- rep(NA_integer_, length(x))
  out[is.finite(x) & x >  eps] <-  1L
  out[is.finite(x) & x < -eps] <- -1L
  out[is.finite(x) & abs(x) <= eps] <- 0L
  out
}

# PGLS fits -------------------------------------------------------------
fit_pgls_ml <- function(tr, df2, dep, pred) {
  corML <- try(ape::corPagel(value = 0.5, phy = tr, fixed = FALSE, form = ~ species), silent = TRUE)
  if (inherits(corML, "try-error")) return(NULL)
  fit <- try(nlme::gls(stats::as.formula(paste(dep, "~", pred)), data = df2,
                       correlation = corML, method = "ML"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  tt <- summary(fit)$tTable
  list(beta = as.numeric(tt[pred, "Value"]),
       p    = as.numeric(tt[pred, "p-value"]),
       t    = as.numeric(tt[pred, "t-value"]),
       AIC  = AIC(fit),
       lambda = get_fitted_lambda(fit),
       fit = fit)
}

fit_pgls_lam0 <- function(tr, df2, dep, pred) {
  cor0 <- try(ape::corPagel(value = 0, phy = tr, fixed = TRUE, form = ~ species), silent = TRUE)
  if (!inherits(cor0, "try-error")) {
    fit0 <- try(nlme::gls(stats::as.formula(paste(dep, "~", pred)), data = df2,
                          correlation = cor0, method = "ML"), silent = TRUE)
    if (!inherits(fit0, "try-error")) {
      tt <- summary(fit0)$tTable
      return(list(beta = as.numeric(tt[pred, "Value"]), p = as.numeric(tt[pred, "p-value"]),
                  t = as.numeric(tt[pred, "t-value"]), AIC = AIC(fit0)))
    }
  }
  # Fallback: OLS if GLS failed (still useful as sensitivity)
  fit0_lm <- try(lm(stats::as.formula(paste(dep, "~", pred)), data = df2), silent = TRUE)
  if (inherits(fit0_lm, "try-error")) return(NULL)
  tt <- summary(fit0_lm)$coefficients
  list(beta = as.numeric(tt[pred, "Estimate"]), p = as.numeric(tt[pred, "Pr(>|t|)"]),
       t = as.numeric(tt[pred, "t value"]), AIC = AIC(fit0_lm))
}

fit_pgls_lam1 <- function(tr, df2, dep, pred) {
  cor1 <- try(ape::corPagel(value = 1, phy = tr, fixed = TRUE, form = ~ species), silent = TRUE)
  if (inherits(cor1, "try-error")) return(NULL)
  fit1 <- try(nlme::gls(stats::as.formula(paste(dep, "~", pred)), data = df2,
                        correlation = cor1, method = "ML"), silent = TRUE)
  if (inherits(fit1, "try-error")) return(NULL)
  tt <- summary(fit1)$tTable
  list(beta = as.numeric(tt[pred, "Value"]), p = as.numeric(tt[pred, "p-value"]),
       t = as.numeric(tt[pred, "t-value"]), AIC = AIC(fit1))
}

run_pgls_one <- function(tree_species, df, dep = "heritability", pred = "presence") {
  # df must have species, dep(=y), pred(=0/1)
  df <- df %>% filter(!is.na(.data[[dep]]), !is.na(.data[[pred]]))
  tab <- table(df[[pred]])
  if (length(tab) < 2 || any(tab < 3)) return(NULL)
  spp <- intersect(tree_species$tip.label, df$species)
  if (length(spp) < 5) return(NULL)
  tr <- keep.tip(tree_species, spp)
  df2 <- df[match(tr$tip.label, df$species), , drop = FALSE]
  rownames(df2) <- as.character(df2$species)
  df2$species <- factor(df2$species, levels = tr$tip.label)
  
  ml <- fit_pgls_ml(tr, df2, dep, pred)
  if (is.null(ml)) return(NULL)
  
  f0 <- fit_pgls_lam0(tr, df2, dep, pred)
  f1 <- fit_pgls_lam1(tr, df2, dep, pred)
  
  AICs <- c(ML = ml$AIC, L0 = if (!is.null(f0)) f0$AIC else NA, L1 = if (!is.null(f1)) f1$AIC else NA)
  best <- which.min(AICs)
  dAIC <- AICs - min(AICs, na.rm = TRUE)
  preferred <- names(AICs)[best]
  
  data.frame(
    pgls_beta = ml$beta, pgls_p = ml$p, t_obs = ml$t,
    lambda_ML = ml$lambda, AIC_ML = ml$AIC,
    beta_lam0 = if (!is.null(f0)) f0$beta else NA_real_,
    p_lam0    = if (!is.null(f0)) f0$p    else NA_real_,
    AIC_lam0  = if (!is.null(f0)) f0$AIC  else NA_real_,
    beta_lam1 = if (!is.null(f1)) f1$beta else NA_real_,
    p_lam1    = if (!is.null(f1)) f1$p    else NA_real_,
    AIC_lam1  = if (!is.null(f1)) f1$AIC  else NA_real_,
    dAIC_ML = dAIC["ML"], dAIC_lam0 = dAIC["L0"], dAIC_lam1 = dAIC["L1"],
    model_preferred = preferred
  )
}

# ------------------------- Load data -------------------------
cat("[load] Reading taxonomy info\n")
taxonomy <- read_tab_file(taxonomy_path, FALSE)
if (grepl("GTDB", taxonomy_path, ignore.case = TRUE) || ncol(taxonomy) == 2) {
  colnames(taxonomy) <- c("Genome","Taxonomy")
  taxonomy$Taxonomy <- gsub("; ",";", taxonomy$Taxonomy)
  taxonomy$Taxonomy <- gsub(" ","_", taxonomy$Taxonomy)
} else if (grepl("Greengenes2", taxonomy_path, ignore.case = TRUE) || ncol(taxonomy) >= 2) {
  colnames(taxonomy)[1:2] <- c("Genome","Taxonomy")
  taxonomy$Taxonomy <- gsub("\\s(?!.*\\s)", "_", taxonomy$Taxonomy, perl = TRUE)
  taxonomy$Taxonomy <- gsub(" ", "", taxonomy$Taxonomy, perl = TRUE)
} else stop("Taxonomy should be GTDB (shotgun) or Greengenes2 (16S) format.")

taxonomy$Phylum  <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 2)
taxonomy$Class   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 3)
taxonomy$Order   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 4)
taxonomy$Family  <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 5)
taxonomy$Genus   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 6)
taxonomy$Species <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 7)

species2genome <- setNames(taxonomy$Genome, taxonomy$Species)
genome2species <- setNames(taxonomy$Species, taxonomy$Genome)

cat("[load] Reading phylogeny: ", phylo_path, "\n")
phy_master <- read.tree(phylo_path)
phy_master <- ensure_branch_lengths(phy_master, method = branch_length_transform)

cat("[load] Loading heritability\n")
load(herit_path)  # must define all_VCs_full
if (!exists("all_VCs_full")) stop("Object 'all_VCs_full' not found in --herit file.")

cat("[load] Loading function matrices\n")
load(func_mats_rda)
if (!isTRUE(identical(func_enog_rda, func_mats_rda))) {
  if (nzchar(func_enog_rda) && file.exists(func_enog_rda)) load(func_enog_rda)
}

# ------------------------- Prepare heritability DF -------------------------
filtered_VCs <- all_VCs_full[grepl(subset_id, all_VCs_full$subset_id), ]
filtered_VCs <- filtered_VCs[grepl("Species", filtered_VCs$rank), ]
filtered_VCs <- filtered_VCs %>% mutate(
  All_Sig = dplyr::case_when(
    Significance %in% c("Significant (Bonferroni)", "Significant (FDR)") ~ "Heritable",
    TRUE ~ "Non-heritable"
  )
)

filtered_VCs$Genome       <- taxonomy$Genome[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Phylum       <- taxonomy$Phylum[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Class        <- taxonomy$Class[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Order        <- taxonomy$Order[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Family       <- taxonomy$Family[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Genus        <- taxonomy$Genus[match(filtered_VCs$trait, taxonomy$Species)]
filtered_VCs$Heritability <- as.numeric(filtered_VCs$var_Ad) * 100
filtered_VCs$Maternal     <- as.numeric(filtered_VCs$var_M) * 100
filtered_VCs$Cohousing    <- as.numeric(filtered_VCs$var_C) * 100

# ------------------------- Main loop -------------------------
cat("[run] Starting loops\n")
DO_CLASSICAL_TTEST <- TRUE
pdf(pdf_out)

make_group_list <- function(gs) lapply(gs, function(g) list(name=g, genera=c(g)))
GROUPS_LIST <- make_group_list(GROUPS)

for (var in VARIABLES) {
  for (grp in GROUPS_LIST) {
    group_name <- grp$name
    cat("\n[run] Var=", var, " | Group=", group_name, "\n", sep="")
    
    # subset species for this group
    this_VCs <- if (is.null(grp$genera)) filtered_VCs else filtered_VCs[filtered_VCs$Genus %in% grp$genera, ]
    if (nrow(this_VCs) < 6) { cat("  -> skip: too few species (", nrow(this_VCs), ")\n", sep=""); next }
    
    spp_group <- unique(this_VCs$trait)
    tr_species <- make_species_tree(phy_master, spp_group, species2genome, genome2species,
                                    branchlen_method = branch_length_transform)
    if (is.null(tr_species) || length(tr_species$tip.label) < 5) {
      cat("  -> skip: tree has <5 tips for this group\n")
      next
    }
    
    for (func in FUNCS) {
      cat("  [func] ", func, "\n", sep="")
      
      # 1) Function presence matrix (genomes as columns)
      mat <- get(paste0("func_matrix_", func))
      if (func == "eggNOG_OGs") {
        message(paste0("Before removing ENOG50 IDs: ", length(rownames(mat))))
        # (Your original had an optional COG-only filter commented out)
        message(paste0("After removing ENOG50 IDs: ", length(rownames(mat))))
      }
      
      # normalize '.1' suffix to match
      x <- data.frame(genome = colnames(mat)[grepl("^GCA_", colnames(mat)) & !grepl("\\.1$", colnames(mat))])
      x$fixed <- paste0(x$genome, ".1")
      y <- data.frame(genome = colnames(mat)[!colnames(mat) %in% x$genome])
      y$fixed <- y$genome
      all_genomes <- rbind(x, y)
      colnames(mat) <- all_genomes$fixed[match(colnames(mat), all_genomes$genome)]
      
      # Map genomes -> Species then transpose to species × functions
      colnames(mat) <- taxonomy$Species[match(colnames(mat), taxonomy$Genome)]
      mat <- t(mat)
      
      # 2) Align rows (species) to this group
      mat <- mat[match(this_VCs$trait, rownames(mat)), , drop = FALSE]
      mat <- mat[!is.na(rownames(mat)), , drop = FALSE]
      this_VCs2 <- this_VCs[this_VCs$trait %in% rownames(mat), ]
      cat("Genomes: ", nrow(mat), "\n", sep="")
      if (nrow(this_VCs2) < 6) { cat("    -> skip: <6 species after align\n"); next }
      
      # 3) Drop all-zero and constant columns
      message(paste0("Before dropping all-zero/constant: ", ncol(mat)))
      n0_pre <- ncol(mat)
      mat <- mat[, colSums(mat) != 0, drop = FALSE]
      mat <- mat[, apply(mat, 2, function(z) !is_constant_column(z)), drop = FALSE]
      if (!ncol(mat)) { cat("    -> skip: all functions dropped (all-zero/constant)\n"); next }
      cat("    kept functions after zero/constant filter: ", ncol(mat), " (from ", n0_pre, ")\n", sep="")
      
      # 4) Presence/absence + minimum group sizes (>=3 per group, or 10% rule)
      presence <- ifelse(mat > 0, 1L, 0L)
      
      # Save presence matrices for extras (genome homogeneity tests)
      if (group_name == "g__Bacteroides") {
        assign("mat_bac", presence, envir = .GlobalEnv)
      } else if (group_name == "g__Prevotella") {
        assign("mat_pre", presence, envir = .GlobalEnv)
      }
      
      ntot <- nrow(mat)
      min_n <- max(3, round(ntot*0.1, digits = 0))
      keep_cols <- apply(presence, 2, function(col) sum(col == 0, na.rm = TRUE) >= min_n && sum(col == 1, na.rm = TRUE) >= min_n)
      presence <- presence[, keep_cols, drop = FALSE]
      if (!ncol(presence)) { cat("    -> skip: none pass min n0/n1 >=", min_n, "\n"); next }
      cat("    functions passing n0>=", min_n, " & n1>=", min_n, ": ", ncol(presence), "\n", sep="")
      
      # 5) Outcome (transformed & raw for interpretability)
      characteristic_df <- data.frame(
        species = this_VCs2$trait,
        y_raw   = this_VCs2[[var]],
        y       = invrank(this_VCs2[[var]])
      )
      
      # 6) Merge presence + outcome
      presence_df <- as.data.frame(presence, check.names = FALSE)
      presence_df$species <- rownames(presence)
      df_base <- dplyr::left_join(presence_df[, c("species"), drop = FALSE], characteristic_df, by = "species")
      df_all  <- cbind(df_base, presence_df[, setdiff(colnames(presence_df), "species"), drop = FALSE])
      
      # 7) Fit one PGLS per function
      res_list <- lapply(colnames(presence), function(fname) {
        df <- df_all[, c("species", "y", "y_raw", fname)]
        colnames(df) <- c("species", "heritability", "heritability_raw", "presence")
        
        out <- run_pgls_one(tr_species, df, dep = "heritability", pred = "presence")
        if (is.null(out)) return(NULL)
        
        # Optional t-test on transformed y
        if (DO_CLASSICAL_TTEST && length(unique(df$presence)) >= 2) {
          tt <- try(t.test(heritability ~ presence, data = df), silent = TRUE)
          if (!inherits(tt, "try-error")) out$ttest_p <- tt$p.value
        }
        
        # Raw-scale summaries (Heritability %) for interpretability
        grp0 <- df$heritability_raw[df$presence == 0]
        grp1 <- df$heritability_raw[df$presence == 1]
        out$mean_raw_absent  <- suppressWarnings(mean(grp0, na.rm = TRUE))
        out$mean_raw_present <- suppressWarnings(mean(grp1, na.rm = TRUE))
        out$diff_raw_present_minus_absent <- out$mean_raw_present - out$mean_raw_absent
        
        # Ns
        tab <- table(df$presence)
        out$n  <- nrow(df)
        out$n0 <- as.integer(tab["0"]); if (is.na(out$n0)) out$n0 <- 0L
        out$n1 <- as.integer(tab["1"]); if (is.na(out$n1)) out$n1 <- 0L
        
        out$Function <- fname
        out
      })
      
      res <- dplyr::bind_rows(res_list)
      if (!nrow(res)) { cat("    -> skip: no valid PGLS rows\n"); next }
      
      # 8) Model-choice fields (AIC-preferred vs ML)
      res <- res %>%
        mutate(
          beta_aic = dplyr::case_when(
            model_preferred == "ML" ~ pgls_beta,
            model_preferred == "L0" ~ beta_lam0,
            model_preferred == "L1" ~ beta_lam1,
            TRUE ~ pgls_beta
          ),
          p_aic = dplyr::case_when(
            model_preferred == "ML" ~ pgls_p,
            model_preferred == "L0" ~ p_lam0,
            model_preferred == "L1" ~ p_lam1,
            TRUE ~ pgls_p
          )
        )
      
      # 9) Decide the reporting model
      if (toupper(REPORT_MODEL) == "AIC") {
        res <- res %>% mutate(beta_used = beta_aic, p_used = p_aic, model_used = model_preferred)
      } else {
        res <- res %>% mutate(beta_used = pgls_beta, p_used = pgls_p, model_used = "ML")
      }
      
      # 10) Directions (first compute ML-only diagnostic)
      res <- res %>%
        mutate(
          Direction_ML = if_else(pgls_beta > 0, paste0("High ", var), paste0("Low ", var)),
          Direction    = if_else(beta_used  > 0, paste0("High ", var), paste0("Low ", var))
        )
      
      # Enforce invariant: if REPORT_MODEL == "ML", then Direction == Direction_ML
      if (toupper(REPORT_MODEL) == "ML") {
        mismatch <- which(res$Direction != res$Direction_ML)
        if (length(mismatch)) {
          res$Direction[mismatch] <- res$Direction_ML[mismatch]
          res$direction_autofixed <- FALSE
          res$direction_autofixed[mismatch] <- TRUE
          cat("    [warn] Direction auto-fixed to ML for ", length(mismatch), " row(s)\n", sep="")
        } else {
          res$direction_autofixed <- FALSE
        }
      } else {
        res$direction_autofixed <- FALSE
      }
      
      # 12) Sign stability diagnostics
      sML <- safe_sign(res$pgls_beta)
      sL0 <- safe_sign(res$beta_lam0)
      sL1 <- safe_sign(res$beta_lam1)
      
      # strict: need all three and all equal
      res$sign_stable_strict <- (sML %in% c(-1L,0L,1L)) & (sL0 %in% c(-1L,0L,1L)) & (sL1 %in% c(-1L,0L,1L)) &
        (sML == sL0) & (sML == sL1)
      
      # available: compare ML only to those that exist
      agree_ML_L0 <- !is.na(sML) & !is.na(sL0) & (sML == sL0)
      agree_ML_L1 <- !is.na(sML) & !is.na(sL1) & (sML == sL1)
      res$sign_stable_avail <- dplyr::case_when(
        !is.na(sL0) & !is.na(sL1) ~ agree_ML_L0 & agree_ML_L1,
        !is.na(sL0) &  is.na(sL1) ~ agree_ML_L0,
        is.na(sL0) & !is.na(sL1) ~ agree_ML_L1,
        TRUE                     ~ NA
      )
      
      # Use ΔAIC for the REPORTED model (not always ML)
      res$dAIC_used <- dplyr::case_when(
        res$model_used == "ML" ~ res$dAIC_ML,
        res$model_used == "L0" ~ res$dAIC_lam0,
        res$model_used == "L1" ~ res$dAIC_lam1,
        TRUE ~ res$dAIC_ML
      )
      
      # 13) P guards + multiple testing on p_used
      res <- res %>%
        mutate(
          p_used = pmin(pmax(p_used, .Machine$double.xmin), 1),
          logP   = -log10(p_used),
          bonf_p = p.adjust(p_used, method = "bonferroni"),
          qvalue = p.adjust(p_used, method = "fdr"),
          Significance = dplyr::case_when(
            dAIC_used <=2 & qvalue < 0.10 & bonf_p < 0.05 ~ "Significant (Bonferroni)",
            dAIC_used <=2 & qvalue < 0.10                 ~ "Significant (FDR)",
            TRUE                           ~ "Non-significant"
          ),
          # Robustness flag: q<0.10, ΔAIC_used<=2, and sign_stable_avail==TRUE
          robust_flag = (!is.na(qvalue) & qvalue < 0.10) &
            (is.na(dAIC_used) | dAIC_used <= 2) &
            (isTRUE(sign_stable_avail))
        )
      
      # 14) Conclusion string (coherent with reported model)
      res <- res %>%
        mutate(
          lambda_str = ifelse(is.finite(lambda_ML), sprintf("%.2f", lambda_ML), "NA"),
          delta_str  = ifelse(is.finite(dAIC_used), sprintf("%.1f", dAIC_used), "NA"),
          beta_str   = ifelse(is.finite(beta_used), format(signif(beta_used, 3), trim = TRUE), "NA"),
          q_str      = ifelse(is.finite(qvalue), formatC(qvalue, format = "e", digits = 2), "NA"),
          Conclusion = paste0(
            ifelse(robust_flag, "Robust PGLS", "PGLS"),
            " (λ=", lambda_str, "; model=", model_used,
            ifelse(is.na(delta_str), "", paste0(", ΔAIC=", delta_str)),
            "): presence → ", Direction, " for ", var,
            " (β=", beta_str, ", q=", q_str, "; n0=", n0, ", n1=", n1,
            "; mean_raw_absent=", sprintf("%.2f", mean_raw_absent),
            ", mean_raw_present=", sprintf("%.2f", mean_raw_present), ")."
          )
        )
      
      # 15) Volcano plot (reported model); Y is -log10(p_used)
      max_lp  <- max(res$logP[is.finite(res$logP)], na.rm = TRUE)
      min_b   <- min(res$beta_used, na.rm = TRUE); max_b <- max(res$beta_used, na.rm = TRUE)
      pad     <- 0.1 * max(abs(c(min_b, max_b)))
      gtitle  <- paste0(gsub("g__","",group_name), " / ", var, " / ", func)
      plt <- ggplot(res, aes(x = beta_used, y = logP)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = 0) +
        geom_point(aes(color = qvalue < 0.10)) +
        scale_color_manual(values = c("darkgrey", "red"), labels = c("Not Significant", "Significant")) +
        theme_bw() +
        xlim(min_b - pad, max_b + pad) +
        ylim(0, max_lp + 1.5) +
        labs(title = gtitle,
             x = "Coefficient (presence -> heritability)",
             y = "-log10(p-value)", color = "") +
        theme(legend.position = "bottom",
              panel.grid.minor = element_line(color = "grey90"),
              axis.line = element_line(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 11), legend.title = element_text(size = 12),
              plot.title = element_text(hjust = 0.5, size = 13))
      print(plt)
      
      # 16) Save both full and filtered tables per-(var,group,func)
      res_full <- res
      
      res_filt <- res %>%
        transmute(
          Function,
          model_used,
          beta_used, p_used, qvalue, bonf_p,
          Direction, Significance,
          lambda_ML, dAIC_used,
          n, n0, n1, mean_raw_absent, mean_raw_present,
          Conclusion
        )
      
      safe_group <- gsub("[^A-Za-z0-9_]+", "_", group_name)
      out_rdata <- file.path(outdir, paste0("enrichment_PGLS_", var, "_", safe_group, "_", func, ".RData"))
      save(res_full, res_filt, file = out_rdata)
      cat("    [save] ", out_rdata, "  (rows=", nrow(res_full), ")\n", sep="")
      
      # Sanity log: how many significant & robust?
      cat("    [sum] sig(FDR<0.10): ", sum(res$Significance!="Non-significant", na.rm=TRUE),
          " | robust: ", sum(res$robust_flag, na.rm=TRUE), "\n", sep="")
      
      # Save significant filtered CSV (same as your original naming convention)
      res_filt_sig <- res_filt[res_filt$Significance!="Non-significant",]
      comb_name <- paste(grp$name, var, func, sep="_")
      csv_out <- file.path(outdir, paste0("res_", comb_name, ".csv"))
      data.table::fwrite(res_filt_sig, file = csv_out, sep = ",", quote = TRUE)
    }
  }
}

dev.off()
cat("[DONE] PGLS enrichment — invariant Direction under ML; res_full + res_filt saved; volcano PDF done.\n")

# ------------------------- Optional post-hoc EXTRAS (original snippets) -------------------------
if (run_extras) {
  cat("[extras] Running optional heatmaps and genome-homogeneity diagnostics\n")
  ok_lib <- function(pkg) { suppressWarnings(suppressPackageStartupMessages(require(pkg, character.only=TRUE))) }
  
  # Heatmaps (requires pheatmap, readxl); expects the two Excel files
  if (ok_lib("pheatmap") && ok_lib("readxl") && nzchar(extras_bac_xlsx) && nzchar(extras_prev_xlsx) &&
      file.exists(extras_bac_xlsx) && file.exists(extras_prev_xlsx)) {
    try({
      enrich_bac <- readxl::read_excel(extras_bac_xlsx)
      bioprocess <- enrich_bac$Bioprocess
      enrich_bac$`Bacteroides - High heritability`[is.na(enrich_bac$`Bacteroides - High heritability`)] <- 0
      enrich_bac$`Bacteroides - Low heritability`[is.na(enrich_bac$`Bacteroides - Low heritability`)] <- 0
      enrich_bac <- as.matrix(enrich_bac[,c(-1)])
      rownames(enrich_bac) <- bioprocess
      enrich_bac <- t(enrich_bac)
      
      enrich_prev <- readxl::read_excel(extras_prev_xlsx)
      bioprocess  <- enrich_prev$Bioprocess
      enrich_prev$`Prevotella - High heritability`[is.na(enrich_prev$`Prevotella - High heritability`)] <- 0
      enrich_prev$`Prevotella - Low heritability`[is.na(enrich_prev$`Prevotella - Low heritability`)] <- 0
      enrich_prev <- as.matrix(enrich_prev[,c(-1)])
      rownames(enrich_prev) <- bioprocess
      enrich_prev <- t(enrich_prev)
      
      col_fun <- colorRampPalette(c("#313695", "white", "#A50026"))
      pdf(extras_pdf_heat)
      pheatmap::pheatmap(enrich_bac,
                         color = col_fun(100), cluster_rows = FALSE, cluster_cols = TRUE,
                         clustering_method = "ward.D2", show_rownames = TRUE, show_colnames = TRUE,
                         fontsize_row = 10, fontsize_col = 10, main = "Number of eggNOG OGs by Bioprocess (Bacteroides)"
      )
      pheatmap::pheatmap(enrich_prev,
                         color = col_fun(100), cluster_rows = FALSE, cluster_cols = TRUE,
                         clustering_method = "ward.D2", show_rownames = TRUE, show_colnames = TRUE,
                         fontsize_row = 10, fontsize_col = 10, main = "Number of eggNOG OGs by Bioprocess (Prevotella)"
      )
      dev.off()
    }, silent = TRUE)
  } else {
    cat("[extras] Heatmaps skipped (missing packages/files)\n")
  }
  
  # Genome homogeneity / Jaccard & PERMDISP (requires vegan, jaccard)
  if (ok_lib("vegan")) {
    try({
      if (exists("mat_bac") && exists("mat_pre")) {
        # Align feature sets across genera
        align_binary_mats <- function(A, B) {
          all_cols <- union(colnames(A), colnames(B))
          fill0 <- function(M) {
            out <- matrix(0L, nrow = nrow(M), ncol = length(all_cols),
                          dimnames = list(rownames(M), all_cols))
            out[, colnames(M)] <- as.matrix(M)
            storage.mode(out) <- "integer"
            out
          }
          list(A = fill0(A), B = fill0(B))
        }
        al <- align_binary_mats(mat_bac, mat_pre)
        A  <- al$A; B <- al$B
        
        # Drop invariant columns across all genomes
        X <- rbind(A, B)
        keep <- apply(X, 2, function(z) length(unique(z)) > 1)
        X <- X[, keep, drop = FALSE]
        
        # Jaccard distances
        D <- vegan::vegdist(X, method = "jaccard")
        
        grp <- factor(c(rep("Bacteroides", nrow(A)), rep("Prevotella", nrow(B))))
        set.seed(42)
        bd <- vegan::betadisper(D, grp, type = "centroid")
        anova_bd <- anova(bd)
        perm_bd  <- vegan::permutest(bd, permutations = 9999)
        
        disp_df <- data.frame(
          genome = names(bd$distances),
          genus  = grp,
          dist_to_centroid = as.numeric(bd$distances)
        )
        by_genus <- aggregate(dist_to_centroid ~ genus, data = disp_df,
                              FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
        by_genus <- do.call(data.frame, by_genus)
        welch_disp <- t.test(dist_to_centroid ~ genus, data = disp_df, var.equal = FALSE)
        
        if (!ok_lib("jaccard")) {
          cat("[extras] 'jaccard' package not available; skipping pairwise jaccard.test\n")
        } else {
          compute_pairwise_jaccard <- function(mat, method = "bootstrap"){
            n <- nrow(mat)
            pvals <- matrix(NA, n, n, dimnames = list(rownames(mat), rownames(mat)))
            stats <- matrix(NA, n, n, dimnames = list(rownames(mat), rownames(mat)))
            for(i in 1:(n - 1)){
              for(j in (i+1):n){
                out <- jaccard::jaccard.test(mat[i,], mat[j,], method = method, verbose = FALSE)
                pvals[i,j] <- out$pvalue; stats[i,j] <- out$statistics
                pvals[j,i] <- out$pvalue; stats[j,i] <- out$statistics
              }
            }
            list(pvalues = pvals, statistics = stats)
          }
          res_bac_mca <- compute_pairwise_jaccard(mat_bac, method = "bootstrap")
          res_pre_mca <- compute_pairwise_jaccard(mat_pre, method = "bootstrap")
          pvec_bac_mca <- na.omit(as.vector(res_bac_mca$pvalues[upper.tri(res_bac_mca$pvalues)]))
          pvec_pre_mca <- na.omit(as.vector(res_pre_mca$pvalues[upper.tri(res_pre_mca$pvalues)]))
          qvec_bac_mca <- p.adjust(pvec_bac_mca, method = "BH")
          qvec_pre_mca <- p.adjust(pvec_pre_mca, method = "BH")
          cat(sprintf("[extras] Mean BH<0.05 proportion: Bacteroides=%.3f, Prevotella=%.3f\n",
                      mean(qvec_bac_mca < 0.05), mean(qvec_pre_mca < 0.05)))
        }
        
        # QQ-plots of distance vectors (optional)
        pdf(extras_pdf_qq)
        dists_bac <- as.vector(as.matrix(D)[1:nrow(A), 1:nrow(A)])
        dists_bac <- dists_bac[upper.tri(matrix(0, nrow=nrow(A), ncol=nrow(A)))]
        dists_pre <- as.vector(as.matrix(D)[(nrow(A)+1):nrow(X), (nrow(A)+1):nrow(X)])
        dists_pre <- dists_pre[upper.tri(matrix(0, nrow=nrow(B), ncol=nrow(B)))]
        qqnorm(dists_pre, main = "Q-Q Plot: Prevotella distances"); qqline(dists_pre, col = "red")
        qqnorm(dists_bac, main = "Q-Q Plot: Bacteroides distances"); qqline(dists_bac, col = "red")
        dev.off()
        
        cat("\n--- PERMDISP (betadisper -> permutest) ---\n")
        print(anova_bd)
        print(perm_bd)
        print(by_genus)
        cat("\nWelch test on per-genome centroid distances (descriptive):\n")
        print(welch_disp)
      } else {
        cat("[extras] Skipped genome-homogeneity analyses: mat_bac/mat_pre not found\n")
      }
    }, silent = TRUE)
  } else {
    cat("[extras] Vegan not available; homogeneity diagnostics skipped\n")
  }
}

# ------------------------------------------------------------------
# End of script
# ------------------------------------------------------------------
