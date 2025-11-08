#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ape)
  library(nlme)
  # library(ggplot2)  # keep disabled for now
})

## ------------------------- Helper utilities -------------------------

read_gzipped   <- function(f) if (grepl("\\.gz$", f)) gzfile(f) else f
read_tab_file  <- function(f, header=FALSE) read.delim(read_gzipped(f), header=header, check.names=FALSE)
read_csv_file  <- function(f) read.csv(read_gzipped(f), stringsAsFactors = FALSE)

ensure_branch_lengths <- function(tr, method = c("none", "grafen")) {
  method <- match.arg(method)
  if (is.null(tr$edge.length)) {
    if (method == "grafen") tr <- compute.brlen(tr, method = "Grafen") else tr$edge.length <- rep(1, nrow(tr$edge))
  }
  if (any(tr$edge.length <= 0)) tr$edge.length[tr$edge.length <= 0] <- 1e-8
  tr
}

# Build a species-labeled subtree from a genome-labeled master tree via taxonomy maps
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

canon <- function(s) sub("^\\s+|\\s+$","", gsub("__+","__", as.character(s)))

# Collapse duplicate (genus,gene,species) rows by inverse-variance weighting when available
collapse_duplicates <- function(dt) {
  keycols <- c("genus","gene","species")
  miss <- setdiff(keycols, names(dt))
  if (length(miss)) stop("collapse_duplicates: missing required columns: ", paste(miss, collapse=","))
  if (!("beta" %in% names(dt))) stop("collapse_duplicates: missing 'beta'")
  if (!("se_beta" %in% names(dt))) dt[, se_beta := NA_real_]
  
  dt[, {
    b  <- beta
    se <- se_beta
    if (all(is.finite(se)) && all(se > 0)) {
      w   <- 1/(se^2)
      b2  <- sum(w * b) / sum(w)
      se2 <- sqrt(1/sum(w))
    } else {
      b2  <- mean(b, na.rm=TRUE)
      se2 <- NA_real_
    }
    .(beta = b2,
      se_beta = se2,
      snp = if ("snp" %in% names(.SD)) snp[1] else NA_character_,
      chr = if ("chr" %in% names(.SD)) as.character(chr[1]) else NA_character_,
      pos = if ("pos" %in% names(.SD)) as.integer(pos[1]) else NA_integer_)
  }, by = keycols]
}

# Robustly fetch slope & p from lm/gls (column names differ)
safe_coef <- function(fit, term) {
  cf <- try(coef(summary(fit)), silent=TRUE)
  if (inherits(cf, "try-error") || is.null(dim(cf)) || !(term %in% rownames(cf))) {
    return(list(slope=NA_real_, p=NA_real_))
  }
  col_est <- if ("Estimate" %in% colnames(cf)) "Estimate" else if ("Value" %in% colnames(cf)) "Value" else NA_character_
  col_p   <- if ("Pr(>|t|)" %in% colnames(cf)) "Pr(>|t|)" else if ("p-value" %in% colnames(cf)) "p-value" else NA_character_
  if (is.na(col_est) || is.na(col_p)) return(list(slope=NA_real_, p=NA_real_))
  list(slope = unname(cf[term, col_est]),
       p     = unname(cf[term, col_p]))
}

## ------------------------- CLI options -------------------------

opt_list <- list(
  make_option("--beta_cf",        type="character"),
  make_option("--herit_rdata",    type="character"),
  make_option("--taxonomy",       type="character"),
  make_option("--phylotree",      type="character"),
  make_option("--tax_ranks",      type="character", default="Phylum,Class,Order,Family,Genus,Species"),
  make_option("--min_species",    type="integer",   default=6),
  make_option("--branch",         type="character", default="genus_sentinels"),
  make_option("--q",              type="double",    default=0.10),
  make_option("--outdir",         type="character", default="."),
  make_option("--log",            type="character", default=NULL,
              help="Optional log file; omit or set to '' or 'NULL' to disable")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$beta_cf), !is.null(opt$herit_rdata))

## ------------------------- logging -------------------------

log_con <- NULL
if (!is.null(opt$log) && nzchar(opt$log) && !identical(opt$log, "NULL")) {
  dir.create(dirname(opt$log), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(opt$log, open="wt")
  sink(log_con, split = TRUE)
  on.exit({ try(sink(), silent=TRUE); try(close(log_con), silent=TRUE) }, add = TRUE)
}
logf <- function(...) message(sprintf("[correlate_h2] %s", sprintf(...)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(opt$outdir, "correlation_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
writeLines("This file ensures the correlation_plots/* glob has a match even if no plots are drawn.",
           file.path(plots_dir, "README.txt"))

## ------------------------- read beta_cf -------------------------

logf("Reading beta_cf: %s", opt$beta_cf)
beta <- fread(opt$beta_cf)
logf("beta_cf rows=%d, cols=%d; names=%s", nrow(beta), ncol(beta), paste(names(beta), collapse=","))

# normalize common header variants when present
rename_if_present <- function(dt, mapping) {
  for (i in seq_len(nrow(mapping))) {
    o <- mapping$old[i]; n <- mapping$new[i]
    if (o %in% names(dt) && !identical(o, n)) setnames(dt, o, n)
  }
}
mapping <- data.table(
  old = c("Beta","BETA","se","SE","SE_BETA","POS","bp","position","trait","taxon"),
  new = c("beta","beta","se_beta","se_beta","se_beta","pos","pos","pos","species","species")
)
rename_if_present(beta, mapping)

# ensure required columns (fill with NA if missing), canonicalize species label
for (nm in setdiff(c("species","genus","gene"), names(beta))) beta[[nm]] <- NA_character_
for (nm in setdiff(c("beta","se_beta"), names(beta)))          beta[[nm]] <- NA_real_
beta[, species := canon(species)]

# collapse duplicates per (genus,gene,species)
dups <- beta[, .N, by=.(genus,gene,species)][N>1]
if (nrow(dups)) logf("Found duplicate rows per (genus,gene,species); collapsing with inverse-variance weighting when possible")
beta <- collapse_duplicates(beta)
logf("After collapsing: beta_cf rows=%d", nrow(beta))

## ------------------------- heritability -------------------------

logf("Loading heritability RData: %s", opt$herit_rdata)
obj_before <- ls()
load(opt$herit_rdata)
obj_after  <- setdiff(ls(), obj_before)
logf("Objects loaded: %s", paste(obj_after, collapse=","))

h2_dt <- NULL
if (exists("all_VCs_full")) {
  x <- get("all_VCs_full")
  cand_trait <- intersect(colnames(x), c("trait","species","taxon"))
  cand_h2    <- intersect(colnames(x), c("var_Ad","h2","H2"))
  if (length(cand_trait) && length(cand_h2)) {
    h2_dt <- as.data.table(x)[, .(species = canon(get(cand_trait[1])),
                                  h2      = as.numeric(get(cand_h2[1])))]
  }
}
if (is.null(h2_dt)) stop("Could not find heritability table with species and h2 in RData")
logf("Herit rows=%d, cols=%d; names=%s", nrow(h2_dt), ncol(h2_dt), paste(names(h2_dt), collapse=","))

## ------------------------- taxonomy & phylogeny -------------------------

if (is.null(opt$taxonomy) || !nzchar(opt$taxonomy) || !file.exists(opt$taxonomy)) {
  stop("No taxonomy file provided or file not found: --taxonomy")
}
cat("[load] Reading taxonomy info\n")
taxonomy <- read_tab_file(opt$taxonomy, FALSE)
colnames(taxonomy) <- c("Genome","Taxonomy")
taxonomy$Taxonomy <- gsub("; ",";", taxonomy$Taxonomy)
taxonomy$Taxonomy <- gsub(" ","_", taxonomy$Taxonomy)
taxonomy$Phylum  <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 2)
taxonomy$Class   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 3)
taxonomy$Order   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 4)
taxonomy$Family  <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 5)
taxonomy$Genus   <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 6)
taxonomy$Species <- sapply(strsplit(taxonomy$Taxonomy, ";"), "[", 7)

species2genome <- setNames(taxonomy$Genome, taxonomy$Species)
genome2species <- setNames(taxonomy$Species, taxonomy$Genome)

# emit a quick mapping preview for debugging
rep_out <- data.table(
  species = sort(unique(beta$species)),
  genome  = unname(species2genome[sort(unique(beta$species))])
)
fwrite(rep_out, file.path(opt$outdir, "genome_representatives.tsv"), sep="\t")
logf("Representative mapping preview (first 8): %s",
     paste(utils::head(sprintf("%s -> %s", rep_out$species, rep_out$genome), 8), collapse=" | "))

cat("[load] Reading phylogeny:\n      ", opt$phylotree, "\n")
phy_master <- read.tree(opt$phylotree)
phy_master <- ensure_branch_lengths(phy_master, method = "none")
logf("Master tree tips: %d", length(phy_master$tip.label))

## ------------------------- join & pre-filter (verbose harmonization) -------------------------

dt <- merge(beta, h2_dt, by="species", all.x = TRUE)
n_miss <- sum(is.na(dt$h2))
logf("Joined: rows=%d; h2 missing for %d species (written to join_misses_beta_cf.txt)", nrow(dt), n_miss)
if (n_miss) writeLines(sort(unique(dt[is.na(h2), species])),
                       file.path(opt$outdir,"join_misses_beta_cf.txt"))

# species present on the tree (via Genome mapping)
species_on_tree <- genome2species[ intersect(phy_master$tip.label, taxonomy$Genome) ]
species_on_tree <- unique(species_on_tree[!is.na(species_on_tree)])
logf("Species on tree (via taxonomy map): %d", length(species_on_tree))

dt_tree <- dt[species %in% species_on_tree]
logf("After restricting to species on tree: rows=%d (from %d); unique species=%d",
     nrow(dt_tree), nrow(dt), uniqueN(dt_tree$species))

# Check finites early and summarize
dt_finite <- dt_tree[is.finite(beta) & is.finite(h2)]
logf("Finite pairs beta/h2: rows=%d; unique species=%d", nrow(dt_finite), uniqueN(dt_finite$species))

# Write a small alignment check table
align_prev <- unique(dt_finite[, .(species, genus)])[order(genus, species)]
fwrite(head(align_prev, 40), file.path(opt$outdir, "alignment_preview.tsv"), sep="\t")

# total species per genus represented on the tree & join
n_total_by_genus <- dt_tree[!is.na(genus) & !is.na(gene),
                            .(n_total_species = uniqueN(species)), by=.(genus)]

# count valid species per (genus,gene) with finite values
valid_counts <- dt_finite[!is.na(genus) & !is.na(gene),
                          .(n_valid = uniqueN(species)), by=.(genus,gene)]
valid_counts <- merge(valid_counts, n_total_by_genus, by="genus", all.x=TRUE)

# keep only full-case groups BEFORE modeling (your request)
full_groups <- valid_counts[n_valid == n_total_species, .(genus,gene)]
logf("Pre-filter: keeping %d (genus,gene) groups with n == total species for the genus; excluding %d",
     nrow(full_groups),
     uniqueN(valid_counts[, .(genus,gene)]) - nrow(full_groups))

# dataset for modeling (only full groups & finite values)
res <- dt_finite[
  full_groups, on=.(genus,gene),
  .(species, beta, se_beta, h2), by = .(genus, gene)
]

## ---------- model fitting (verbose, PGLS-first) ----------

logf("Preparing groups …")
ug <- uniqueN(res[, .(genus,gene)])
logf("Groups: %d unique (genus,gene) after pre-filter", ug)

n_by_grp <- res[, .N, by=.(genus,gene)]
if (nrow(n_by_grp)) {
  logf("  summary N per (genus,gene): min=%d, p25=%d, median=%d, p75=%d, max=%d",
       min(n_by_grp$N),
       as.integer(quantile(n_by_grp$N, .25)),
       as.integer(median(n_by_grp$N)),
       as.integer(quantile(n_by_grp$N, .75)),
       max(n_by_grp$N))
}

fail_log <- list()
fits_list <- vector("list", max(1L, ug))
gi <- 1L

setkey(res, genus, gene)
grp_keys <- unique(res[, .(genus, gene)])

for (k in seq_len(nrow(grp_keys))) {
  g  <- grp_keys$genus[k]
  ge <- grp_keys$gene[k]
  dsub <- res[list(g, ge)]
  n <- nrow(dsub)
  
  n_total <- n_total_by_genus[genus==g]$n_total_species
  if (length(n_total) == 0L) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="no n_total_species found for genus")
    next
  }
  if (n != n_total) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason=sprintf("n=%d != n_total_species=%d", n, n_total))
    next
  }
  if (n < opt$min_species) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason=sprintf("n=%d < min_species=%d", n, opt$min_species))
    next
  }
  if (var(dsub$h2) <= 0 || var(dsub$beta) <= 0) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="zero variance in h2 or beta (pre-tree)")
    next
  }
  
  # summary stats & directions
  med_beta <- median(dsub$beta, na.rm=TRUE)
  beta_dir <- if (!is.finite(med_beta) || med_beta==0) "Zero" else if (med_beta>0) "Positive" else "Negative"
  
  ## ---- OLS (centered; slope invariant to centering) ----
  dsub$h2_c   <- as.numeric(scale(dsub$h2,   center = TRUE, scale = FALSE))
  dsub$beta_c <- as.numeric(scale(dsub$beta, center = TRUE, scale = FALSE))
  use_w_ols   <- "se_beta" %in% names(dsub) && all(is.finite(dsub$se_beta)) && all(dsub$se_beta > 0)
  ols <- if (use_w_ols) lm(beta_c ~ h2_c, data = dsub, weights = 1/(se_beta^2)) else lm(beta_c ~ h2_c, data = dsub)
  c1  <- safe_coef(ols, "h2_c")
  ols_row <- data.table(
    genus = g, gene = ge, model = "OLS",
    n = n,
    slope = c1$slope, p = c1$p,
    adjR2 = tryCatch(summary(ols)$adj.r.squared, error=function(e) NA_real_),
    lambda = NA_real_, aic = tryCatch(AIC(ols), error=function(e) NA_real_),
    Median_beta_value = med_beta,
    Beta_direction = beta_dir,
    Correlation_direction = if (!is.finite(c1$slope) || c1$slope==0) "Zero" else if (c1$slope>0) "Positive" else "Negative"
  )
  
  ## ---- PGLS (Pagel's lambda) ----
  tr <- try(make_species_tree(phy_master, dsub$species, species2genome, genome2species, "none"), silent=TRUE)
  have_tr <- !(inherits(tr, "try-error") || is.null(tr))
  if (!have_tr) {
    logf("[PGLS skip] genus=%s gene=%s : failed to build tree", g, ge)
    fits_list[[gi]] <- rbindlist(list(ols_row), use.names=TRUE, fill=TRUE); gi <- gi + 1L; next
  }
  
  sp <- tr$tip.label
  d2 <- dsub[match(sp, dsub$species)]
  # align & guard finites
  keep_fin <- is.finite(d2$h2) & is.finite(d2$beta)
  if (!all(keep_fin)) d2 <- d2[keep_fin,]
  if (nrow(d2) != n_total) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason=sprintf("post-tree n=%d != n_total_species=%d", nrow(d2), n_total))
    fits_list[[gi]] <- rbindlist(list(ols_row), use.names=TRUE, fill=TRUE); gi <- gi + 1L; next
  }
  if (var(d2$h2) <= 0 || var(d2$beta) <= 0) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="zero variance in h2 or beta (post-tree)")
    fits_list[[gi]] <- rbindlist(list(ols_row), use.names=TRUE, fill=TRUE); gi <- gi + 1L; next
  }
  
  # center on aligned set; set species factor aligned to tree; use data.frame for nlme
  d2$h2_c   <- as.numeric(scale(d2$h2,   center = TRUE, scale = FALSE))
  d2$beta_c <- as.numeric(scale(d2$beta, center = TRUE, scale = FALSE))
  rownames(d2) <- d2$species
  d2$species <- factor(d2$species, levels = tr$tip.label)
  d2 <- as.data.frame(d2)
  
  # Use weights if strictly valid; IMPORTANT: varFixed(~ se_beta^2) (NOT 1/se^2)
  use_w_gls <- "se_beta" %in% names(d2) && all(is.finite(d2$se_beta)) && all(d2$se_beta > 0)
  
  pgls_rows <- list()
  try_gls_once <- function(cor_struct, label, weighted_ok = use_w_gls) {
    ctrl <- glsControl(msMaxIter = 200, msVerbose = FALSE)
    attempt <- try(gls(beta_c ~ h2_c, data = d2,
                       correlation = cor_struct,
                       weights = if (weighted_ok) varFixed(~ I(se_beta^2)) else NULL,
                       method="ML", na.action=na.omit, control = ctrl), silent = TRUE)
    if (!inherits(attempt, "try-error")) return(list(ok=TRUE, fit=attempt, weighted=weighted_ok))
    if (weighted_ok) {
      attempt2 <- try(gls(beta_c ~ h2_c, data = d2,
                          correlation = cor_struct,
                          weights = NULL,
                          method="ML", na.action=na.omit, control = ctrl), silent = TRUE)
      if (!inherits(attempt2, "try-error")) return(list(ok=TRUE, fit=attempt2, weighted=FALSE))
    }
    list(ok=FALSE, err=conditionMessage(attr(attempt, "condition")))
  }
  
  # λ free (multiple starts for stability) and with explicit species covariate
  for (start_lam in c(0.3, 0.5, 0.8)) {
    cor_free <- try(corPagel(value = start_lam, phy = tr, fixed = FALSE, form = ~ species), silent = TRUE)
    if (inherits(cor_free, "try-error")) next
    gls_try <- try_gls_once(cor_free, sprintf("PGLS_Pagel(start=%.1f)", start_lam))
    if (isTRUE(gls_try$ok)) {
      lam <- tryCatch(as.numeric(coef(cor_free, unconstrained = FALSE)[1]), error=function(e) NA_real_)
      c2  <- safe_coef(gls_try$fit, "h2_c")
      pgls_rows[["PGLS_Pagel"]] <- data.table(
        genus = g, gene = ge, model = "PGLS_Pagel",
        n = nrow(d2), slope = c2$slope, p = c2$p,
        adjR2 = NA_real_, lambda = lam, aic = tryCatch(AIC(gls_try$fit), error=function(e) NA_real_),
        Median_beta_value = med_beta,
        Beta_direction = beta_dir,
        Correlation_direction = if (!is.finite(c2$slope) || c2$slope==0) "Zero" else if (c2$slope>0) "Positive" else "Negative"
      )
      break
    } else {
      err_msg <- gls_try$err; if (is.null(err_msg)) err_msg <- "unknown error"
      logf("[PGLS fail λ-free GLS] genus=%s gene=%s : %s (start λ=%.1f)", g, ge, err_msg, start_lam)
    }
  }
  
  # λ = 0 (no phylo correlation)
  cor0 <- try(corPagel(value = 0, phy = tr, fixed = TRUE, form = ~ species), silent = TRUE)
  if (!inherits(cor0, "try-error")) {
    gls0 <- try_gls_once(cor0, "PGLS_lambda0")
    if (isTRUE(gls0$ok)) {
      c0 <- safe_coef(gls0$fit, "h2_c")
      pgls_rows[["PGLS_lambda0"]] <- data.table(
        genus = g, gene = ge, model = "PGLS_lambda0",
        n = nrow(d2), slope = c0$slope, p = c0$p,
        adjR2 = NA_real_, lambda = 0, aic = tryCatch(AIC(gls0$fit), error=function(e) NA_real_),
        Median_beta_value = med_beta,
        Beta_direction = beta_dir,
        Correlation_direction = if (!is.finite(c0$slope) || c0$slope==0) "Zero" else if (c0$slope>0) "Positive" else "Negative"
      )
    } else {
      err_msg <- gls0$err; if (is.null(err_msg)) err_msg <- "unknown error"
      logf("[PGLS fail λ=0 GLS] genus=%s gene=%s : %s", g, ge, err_msg)
    }
  } else {
    logf("[PGLS fail λ=0] genus=%s gene=%s : could not create corPagel", g, ge)
  }
  
  # λ = 1 (Brownian motion)
  cor1 <- try(corPagel(value = 1, phy = tr, fixed = TRUE, form = ~ species), silent = TRUE)
  if (!inherits(cor1, "try-error")) {
    gls1 <- try_gls_once(cor1, "PGLS_lambda1")
    if (isTRUE(gls1$ok)) {
      c1x <- safe_coef(gls1$fit, "h2_c")
      pgls_rows[["PGLS_lambda1"]] <- data.table(
        genus = g, gene = ge, model = "PGLS_lambda1",
        n = nrow(d2), slope = c1x$slope, p = c1x$p,
        adjR2 = NA_real_, lambda = 1, aic = tryCatch(AIC(gls1$fit), error=function(e) NA_real_),
        Median_beta_value = med_beta,
        Beta_direction = beta_dir,
        Correlation_direction = if (!is.finite(c1x$slope) || c1x$slope==0) "Zero" else if (c1x$slope>0) "Positive" else "Negative"
      )
    } else {
      err_msg <- gls1$err; if (is.null(err_msg)) err_msg <- "unknown error"
      logf("[PGLS fail λ=1 GLS] genus=%s gene=%s : %s", g, ge, err_msg)
    }
  } else {
    logf("[PGLS fail λ=1] genus=%s gene=%s : could not create corPagel", g, ge)
  }
  
  rows <- c(list(ols_row), pgls_rows)
  fits_list[[gi]] <- rbindlist(rows, use.names=TRUE, fill=TRUE)
  gi <- gi + 1L
}

fits <- if (gi > 1L) rbindlist(fits_list[seq_len(gi-1L)], use.names=TRUE, fill=TRUE) else data.table()
logf("Models fit: %d rows", nrow(fits))

## ------------------------- add n_total, significance, outputs -------------------------

# Attach n_total_species per genus
fits <- merge(fits, n_total_by_genus, by="genus", all.x=TRUE)

# p-guards and corrections
fits[, p := pmin(pmax(p, .Machine$double.xmin), 1)]
fits[, bonf_p := p.adjust(p, method = "bonferroni")]
fits[, qvalue := p.adjust(p, method = "fdr")]
fits[, Significance := fifelse(qvalue < 0.10 & bonf_p < 0.05, "Significant (Bonferroni)",
                               fifelse(qvalue < 0.10, "Significant (FDR)", "Non-significant"))]

# Main output (already full-case because of prefilter)
main_out <- fits

# Exclusions (for reference): groups that failed preconditions (n != n_total, etc.)
excl_keys <- unique(valid_counts[n_valid != n_total_species, .(genus,gene)])
if (length(fail_log)) {
  flog <- rbindlist(fail_log, use.names=TRUE, fill=TRUE)
  excl_log <- merge(excl_keys, flog, by=c("genus","gene"), all=TRUE)
} else {
  excl_log <- excl_keys
}

## ------------------------- write outputs -------------------------

out_file <- file.path(opt$outdir, sprintf("correlations_by_gene__%s.tsv", opt$branch))
if (!nrow(main_out)) {
  main_out <- data.table(
    genus=character(), gene=character(), model=character(),
    n=integer(), slope=numeric(), p=numeric(),
    adjR2=numeric(), lambda=numeric(), aic=numeric(),
    Median_beta_value = numeric(),
    Beta_direction = character(),
    Correlation_direction = character(),
    n_total_species = integer(),
    bonf_p = numeric(),
    qvalue = numeric(),
    Significance = character()
  )
}
fwrite(main_out, out_file, sep="\t")
logf("Wrote: %s (rows=%d)", out_file, nrow(main_out))

ex_file <- file.path(opt$outdir, sprintf("correlations_by_gene__%s__EXCLUDED.tsv", opt$branch))
fwrite(excl_log, ex_file, sep="\t")
logf("Wrote: %s (rows=%d)", ex_file, nrow(excl_log))

if (length(fail_log)) {
  fwrite(rbindlist(fail_log, use.names=TRUE, fill=TRUE),
         file.path(opt$outdir, sprintf("pgls_failures__%s.tsv", opt$branch)), sep="\t")
}

logf("DONE")
