#!/usr/bin/env Rscript
## ============================================================================
## correlate_h2_simple.R  (updated – species-space tree + safe empty-handling)
## - Keep syntax basic and defensive
## - Build phylogeny subtrees directly in SPECIES space (no genome mapping)
## - Deduplicate species tips in the relabeled master tree
## - Print before/after counts everywhere (beta_cf, join, tree, groups)
## - Clean exit if no successful fits (prevents data.table merge error)
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ape)
  library(nlme)
  library(ggplot2)
})

## ------------------------- small utils -------------------------
read_gz  <- function(f) if (grepl("\\.gz$", f)) gzfile(f) else f
read_tab <- function(f, header=FALSE) read.delim(read_gz(f), header=header, check.names=FALSE)
canon    <- function(s) sub("^\\s+|\\s+$","", gsub("__+","__", as.character(s)))
invrank  <- function(x) qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(is.finite(x)))

ensure_branch_lengths <- function(tr, method=c("none","grafen")){
  method <- match.arg(method)
  if (is.null(tr$edge.length)) {
    if (method=="grafen") {
      tr <- compute.brlen(tr, method="Grafen")            # add Grafen lengths
    } else {
      tr$edge.length <- rep(1, nrow(tr$edge))             # unit lengths
    }
  }
  if (any(tr$edge.length <= 0)) tr$edge.length[tr$edge.length <= 0] <- 1e-8
  tr
}

## NEW: species-space subtree builder (no genome mapping)
make_species_tree <- function(master_tree, species_vec, branchlen_method="none"){
  # Work directly with species labels (master_tree$tip.label are species)
  sp <- intersect(master_tree$tip.label, unique(as.character(species_vec)))
  if (length(sp) < 3) return(NULL)
  sub <- keep.tip(master_tree, sp)
  # Deduplicate species tips if master tree still contains repeated species labels
  if (anyDuplicated(sub$tip.label)) {
    keep <- !duplicated(sub$tip.label)
    sub  <- keep.tip(sub, sub$tip.label[keep])
  }
  if (length(sub$tip.label) < 3) return(NULL)
  ensure_branch_lengths(sub, method=branchlen_method)
}

safe_sign <- function(x, eps=1e-12){
  out <- rep(NA_integer_, length(x))
  out[is.finite(x) & x >  eps] <-  1L
  out[is.finite(x) & x < -eps] <- -1L
  out[is.finite(x) & abs(x) <= eps] <- 0L
  out
}
sanitize <- function(s) gsub("[^A-Za-z0-9_.-]", "_", as.character(s))

## ------------------------- CLI -------------------------
opt_list <- list(
  make_option("--beta_cf",     type="character"),
  make_option("--herit_rdata", type="character"),
  make_option("--taxonomy",    type="character"),
  make_option("--phylotree",   type="character"),
  make_option("--tax_ranks",   type="character", default="Phylum,Class,Order,Family,Genus,Species"),
  make_option("--min_species", type="integer",   default=6),
  make_option("--branch",      type="character", default="genus_sentinels"),
  make_option("--q",           type="double",    default=0.10),
  make_option("--outdir",      type="character", default="."),
  # toggles
  make_option("--transform",   type="character", default="rankz", help="none | z | rankz"),
  make_option("--report",      type="character", default="ML",    help="ML | AIC"),
  make_option("--branchlen",   type="character", default="none",  help="none | grafen")
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$beta_cf), !is.null(opt$herit_rdata), !is.null(opt$taxonomy), !is.null(opt$phylotree))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
plots_dir <- file.path(opt$outdir,"correlation_plots"); dir.create(plots_dir, showWarnings=FALSE, recursive=TRUE)

LOG_EVERY <- as.integer(Sys.getenv("LOG_EVERY", "25"))  # progress cadence (set 0 to disable)
logf <- function(...) message(sprintf("[correlate_h2_simple] %s", sprintf(...)))

logf("START")
logf("Args: --beta_cf=%s --herit_rdata=%s --taxonomy=%s --phylotree=%s", opt$beta_cf, opt$herit_rdata, opt$taxonomy, opt$phylotree)
logf("Settings: transform=%s report=%s branchlen=%s min_species=%d outdir=%s q=%.3f",
     opt$transform, opt$report, opt$branchlen, opt$min_species, opt$outdir, opt$q)
if (LOG_EVERY > 0) logf("Progress logging every %d groups (override with env LOG_EVERY)", LOG_EVERY)

## ------------------------- load beta_cf -------------------------
logf("Reading beta_cf: %s", opt$beta_cf)
beta <- data.table::fread(opt$beta_cf)
beta_ori <- beta
logf("beta_cf columns: %s", paste(names(beta), collapse=", "))

## harmonize expected columns (keep it simple + explicit)
if ("trait" %in% names(beta)) setnames(beta, "trait", "species")
beta <- beta[!grepl("s__Bacteroides_F_", beta$species), ] # your previous filter
beta[, genus := gsub("^s__","", species)]
beta[, genus := sapply(strsplit(genus, "_", fixed=TRUE), `[`, 1)]
beta[, gene  := paste(chr, pos, sep=":")]

logf("beta_cf rows=%d | unique species=%d | unique genera=%d | unique genes=%d",
     nrow(beta), uniqueN(beta$species), uniqueN(beta$genus), uniqueN(beta$gene))

# collapse duplicates (genus,gene,species) by mean of beta_mean (simple path)
dups <- beta[, .N, by=.(genus,gene,species)][N>1]
if (nrow(dups)) logf("Found %d duplicate (genus,gene,species) entries; collapsing by mean(beta_mean)", nrow(dups))
beta <- beta[, .(beta = mean(beta_mean, na.rm=TRUE)), by=.(genus,gene,species)]
logf("After collapse: rows=%d", nrow(beta))

## ------------------------- heritability -------------------------
logf("Loading RData: %s", opt$herit_rdata)
obj_before <- ls(); load(opt$herit_rdata); obj_after <- setdiff(ls(), obj_before)
logf("Loaded objects: %s", paste(obj_after, collapse=", "))
if (!exists("all_VCs_full")) stop("Object 'all_VCs_full' not found in RData")
x <- get("all_VCs_full")

tcols <- intersect(colnames(x), c("trait","species","taxon"))
hcols <- intersect(colnames(x), c("var_Ad","h2","H2"))
if (!(length(tcols) && length(hcols))) stop("No suitable columns in all_VCs_full")
h2_dt <- as.data.table(x)[, .(species = canon(get(tcols[1])),
                              h2      = as.numeric(get(hcols[1])))]
logf("Heritability table: rows=%d, unique species=%d", nrow(h2_dt), uniqueN(h2_dt$species))
if (sum(is.finite(h2_dt$h2)) > 0) {
  hs <- summary(na.omit(h2_dt$h2))
  logf("h2 summary (finite): Min=%.3f Q1=%.3f Median=%.3f Mean=%.3f Q3=%.3f Max=%.3f",
       hs[1], hs[2], hs[3], mean(h2_dt$h2, na.rm=TRUE), hs[5], hs[6])
}

## ------------------------- taxonomy / tree -------------------------
logf("Reading taxonomy: %s", opt$taxonomy)
taxonomy <- read_tab(opt$taxonomy, FALSE)
colnames(taxonomy) <- c("Genome","Taxonomy")
taxonomy$Taxonomy <- gsub("; ",";", taxonomy$Taxonomy)
taxonomy$Taxonomy <- gsub(" ","_", taxonomy$Taxonomy)
taxonomy$Species  <- sapply(strsplit(taxonomy$Taxonomy, ";", fixed=TRUE), `[`, 7)
logf("Taxonomy loaded: rows=%d | unique genomes=%d | unique species=%d",
     nrow(taxonomy), uniqueN(taxonomy$Genome), uniqueN(taxonomy$Species))

logf("Reading phylogeny: %s", opt$phylotree)
phy_master <- read.tree(opt$phylotree)
logf("Master tree: tips=%d, nodes=%d; ensuring branch lengths (method=%s)",
     length(phy_master$tip.label), phy_master$Nnode, opt$branchlen)
phy_master <- ensure_branch_lengths(phy_master, method=opt$branchlen)

## RELABEL master tree tips from genome IDs -> species names
phy_master$tip.label <- taxonomy$Species[match(phy_master$tip.label, taxonomy$Genome)]

## RESTRICT tree to the species that appear in beta (and DEDUP species tips)
present <- unique(beta$species)
drop_these <- phy_master$tip.label[!(phy_master$tip.label %in% present)]
phy_master <- drop.tip(phy_master, drop_these)

if (anyDuplicated(phy_master$tip.label)) {
  keep <- !duplicated(phy_master$tip.label)
  phy_master <- drop.tip(phy_master, phy_master$tip.label[!keep])
}

logf("After relabel+restrict+dedup: unique species on tree=%d", length(phy_master$tip.label))

## ------------------------- join & restrict -------------------------
dt <- merge(beta, h2_dt, by="species", all.x=TRUE)
n_na_h2 <- sum(is.na(dt$h2))
if (n_na_h2) {
  misses <- sort(unique(dt[is.na(h2), species]))
  writeLines(misses, file.path(opt$outdir,"join_misses_beta_cf.txt"))
  logf("Species with missing h2 after join: %d (written join_misses_beta_cf.txt)", length(misses))
  dt <- dt[!is.na(h2)]
}
logf("After join: rows=%d; unique species=%d", nrow(dt), data.table::uniqueN(dt$species))

# full-case rule
n_total_by_genus <- dt[!is.na(genus) & !is.na(gene), .(n_total_species = uniqueN(species)), by=.(genus)]
valid_counts <- dt[is.finite(beta) & is.finite(h2) & !is.na(genus) & !is.na(gene),
                   .(n_valid = uniqueN(species)), by=.(genus,gene)]
valid_counts <- merge(valid_counts, n_total_by_genus, by="genus", all.x=TRUE)
full_groups <- valid_counts[n_valid == n_total_species, .(genus,gene)]
logf("Group counts: total groups=%d | full-case eligible=%d | excluded=%d",
     nrow(valid_counts), nrow(full_groups), nrow(valid_counts) - nrow(full_groups))

res <- dt[is.finite(beta) & is.finite(h2)][full_groups, on=.(genus,gene)]
if (!nrow(res)) {
  logf("No groups left after full-case filtering. Writing empty outputs and exiting.")
  fwrite(data.table(), file.path(opt$outdir, sprintf("correlations_by_gene__%s.tsv", opt$branch)), sep="\t")
  fwrite(valid_counts[n_valid != n_total_species, .(genus,gene)],
         file.path(opt$outdir, sprintf("correlations_by_gene__%s__EXCLUDED.tsv", opt$branch)), sep="\t")
  quit(save="no")
}

## ------------------------- transform (robust) -------------------------
transform_pair <- function(b, h, what){
  if (what=="rankz") list(beta = invrank(b), h2 = invrank(h))
  else if (what=="z") list(beta = as.numeric(scale(b)), h2 = as.numeric(scale(h)))
  else list(beta = as.numeric(b), h2 = as.numeric(h))
}
logf("Transform: %s (applied to beta and h2)", opt$transform)

## ------------------------- PGLS fits (no weights) -------------------------
fit_one <- function(tr, d2){
  # d2 has species, beta, h2
  d2 <- d2[is.finite(beta) & is.finite(h2), ]
  if (nrow(d2) < 5) return(NULL)
  sp <- intersect(tr$tip.label, d2$species)
  if (length(sp) < 5) return(NULL)
  tr2 <- keep.tip(tr, sp)
  d3  <- d2[match(tr2$tip.label, d2$species), , drop=FALSE]
  rownames(d3) <- as.character(d3$species)
  d3$species   <- factor(d3$species, levels=tr2$tip.label)
  
  coef_gls <- function(fit, term) {
    tt <- try(summary(fit)$tTable, silent=TRUE)
    if (inherits(tt, "try-error") || !(term %in% rownames(tt))) return(c(NA_real_, NA_real_))
    c(as.numeric(tt[term,"Value"]), as.numeric(tt[term,"p-value"]))
  }
  
  # ML (free λ)
  ml <- try({
    cML <- ape::corPagel(value=0.5, phy=tr2, fixed=FALSE, form=~species)
    gls(beta ~ h2, data=d3, correlation=cML, method="ML", control=glsControl(msMaxIter=200))
  }, silent=TRUE)
  beta_ml <- p_ml <- lam_ml <- aic_ml <- NA_real_
  if (!inherits(ml, "try-error")) {
    bp <- coef_gls(ml, "h2"); beta_ml <- bp[1]; p_ml <- bp[2]
    lam_ml <- tryCatch(as.numeric(coef(ml$modelStruct$corStruct, unconstrained=FALSE)[1]), error=function(e) NA_real_)
    aic_ml <- tryCatch(AIC(ml), error=function(e) NA_real_)
  }
  
  # λ=0
  l0 <- try({
    c0 <- ape::corPagel(value=0, phy=tr2, fixed=TRUE, form=~species)
    gls(beta ~ h2, data=d3, correlation=c0, method="ML")
  }, silent=TRUE)
  beta_l0 <- p_l0 <- aic_l0 <- NA_real_
  if (!inherits(l0, "try-error")) {
    bp <- coef_gls(l0, "h2"); beta_l0 <- bp[1]; p_l0 <- bp[2]
    aic_l0 <- tryCatch(AIC(l0), error=function(e) NA_real_)
  }
  
  # λ=1
  l1 <- try({
    c1 <- ape::corPagel(value=1, phy=tr2, fixed=TRUE, form=~species)
    gls(beta ~ h2, data=d3, correlation=c1, method="ML")
  }, silent=TRUE)
  beta_l1 <- p_l1 <- aic_l1 <- NA_real_
  if (!inherits(l1, "try-error")) {
    bp <- coef_gls(l1, "h2"); beta_l1 <- bp[1]; p_l1 <- bp[2]
    aic_l1 <- tryCatch(AIC(l1), error=function(e) NA_real_)
  }
  
  # AIC choice
  AICs <- c(ML=aic_ml, L0=aic_l0, L1=aic_l1)
  best <- names(which.min(AICs))
  dAIC <- AICs - min(AICs, na.rm=TRUE)
  
  # pick report
  if (toupper(opt$report)=="AIC" && is.finite(AICs[best])) {
    if (best=="ML") { beta_used <- beta_ml; p_used <- p_ml }
    if (best=="L0") { beta_used <- beta_l0; p_used <- p_l0 }
    if (best=="L1") { beta_used <- beta_l1; p_used <- p_l1 }
    model_used <- best; dAIC_used <- dAIC[best]
  } else {
    beta_used <- beta_ml; p_used <- p_ml; model_used <- "ML"; dAIC_used <- dAIC["ML"]
  }
  
  data.frame(
    pgls_beta = beta_ml, pgls_p = p_ml, lambda_ML = lam_ml, AIC_ML = aic_ml,
    beta_lam0 = beta_l0, p_lam0 = p_l0, AIC_lam0 = aic_l0,
    beta_lam1 = beta_l1, p_lam1 = p_l1, AIC_lam1 = aic_l1,
    beta_used = beta_used, p_used = p_used, model_used = model_used,
    dAIC_ML = dAIC["ML"], dAIC_lam0 = dAIC["L0"], dAIC_lam1 = dAIC["L1"], dAIC_used = dAIC_used
  )
}

## ------------------------- loop per (genus,gene) -------------------------
setkey(res, genus, gene)
groups <- unique(res[, .(genus,gene)])
out_list <- vector("list", nrow(groups))
fail_log <- list()

logf("Fitting PGLS across groups: %d (genus,gene) pairs", nrow(groups))
if (nrow(groups) > 0) {
  head_g <- head(groups, n=min(5L, nrow(groups)))
  logf("First groups: %s", paste(sprintf("(%s,%s)", head_g$genus, head_g$gene), collapse=" | "))
}

for (i in seq_len(nrow(groups))) {
  g  <- groups$genus[i]
  ge <- groups$gene[i]
  d0 <- res[list(g, ge)][, .(species, beta, h2)]
  
  # transform
  trp <- transform_pair(d0$beta, d0$h2, opt$transform)
  d0$beta <- trp$beta; d0$h2 <- trp$h2
  
  # species-space subtree (now consistent with relabeled master tree)
  tr <- try(make_species_tree(phy_master, d0$species, opt$branchlen), silent=TRUE)
  if (inherits(tr,"try-error") || is.null(tr)) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="tree_build_failed")
    if (LOG_EVERY > 0 && (i %% LOG_EVERY == 0L)) logf("(%d/%d) %s | %s -> tree_build_failed", i, nrow(groups), g, ge)
    next
  }
  
  # variance guards
  if (var(d0$beta, na.rm=TRUE) <= 0 || var(d0$h2, na.rm=TRUE) <= 0) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="zero_variance")
    if (LOG_EVERY > 0 && (i %% LOG_EVERY == 0L)) logf("(%d/%d) %s | %s -> zero_variance", i, nrow(groups), g, ge)
    next
  }
  
  fit <- fit_one(tr, d0)
  if (is.null(fit)) {
    fail_log[[length(fail_log)+1]] <- data.table(genus=g, gene=ge, reason="gls_fit_failed")
    if (LOG_EVERY > 0 && (i %% LOG_EVERY == 0L)) logf("(%d/%d) %s | %s -> gls_fit_failed", i, nrow(groups), g, ge)
    next
  }
  
  fit$genus <- g; fit$gene <- ge
  out_list[[i]] <- fit
  
  if (LOG_EVERY > 0 && (i %% LOG_EVERY == 0L)) {
    logf("(%d/%d) %s | %s -> model=%s, beta=%.4f, p=%.3g (dAIC used=%.3f)",
         i, nrow(groups), g, ge, fit$model_used, fit$beta_used, fit$p_used, fit$dAIC_used)
  }
}

out <- data.table::rbindlist(out_list, use.names=TRUE, fill=TRUE)
if (is.null(out) || !nrow(out)) {
  logf("Finished fits: OK=0 | failed=%d", length(fail_log))
  # Write empty but valid outputs and failure details, then exit safely
  main_path <- file.path(opt$outdir, sprintf("correlations_by_gene__%s.tsv", opt$branch))
  excl_path <- file.path(opt$outdir, sprintf("correlations_by_gene__%s__EXCLUDED.tsv", opt$branch))
  fwrite(data.table(), main_path, sep="\t")
  excl <- valid_counts[n_valid != n_total_species, .(genus,gene)]
  fwrite(excl, excl_path, sep="\t")
  if (length(fail_log)) {
    fail_dt <- data.table::rbindlist(fail_log, use.names=TRUE, fill=TRUE)
    fail_path <- file.path(opt$outdir, sprintf("pgls_failures__%s.tsv", opt$branch))
    fwrite(fail_dt, fail_path, sep="\t")
    tally <- fail_dt[, .N, by=reason][order(-N)]
    logf("Failure breakdown: %s", paste(sprintf("%s=%d", tally$reason, tally$N), collapse=" | "))
    logf("Wrote failure details: %s", fail_path)
  }
  logf("DONE (no successful fits)")
  quit(save="no")
}

n_fit_ok <- nrow(out)
n_fail   <- length(fail_log)
logf("Finished fits: OK=%d | failed=%d", n_fit_ok, n_fail)

# attach n_total_species and basic directions
out <- merge(out, n_total_by_genus, by="genus", all.x=TRUE)
out[, p_used := pmin(pmax(p_used, .Machine$double.xmin), 1)]
out[, logP := -log10(p_used)]
out[, bonf_p := p.adjust(p_used, method = "bonferroni")]
out[, qvalue := p.adjust(p_used, method = "fdr")]
out[, Direction := fifelse(beta_used > 0, "Positive", fifelse(beta_used < 0, "Negative", "Zero"))]
out[, Significance := fifelse(qvalue < opt$q & bonf_p < 0.05, "Significant (Bonferroni)",
                              fifelse(qvalue < opt$q, "Significant (FDR)", "Non-significant"))]

# write main/excluded/fail logs
main_out <- out[, .(genus,gene,model_used,beta_used,p_used,logP,qvalue,bonf_p,
                    Direction,Significance,
                    lambda_ML,AIC_ML,beta_lam0,p_lam0,AIC_lam0,
                    beta_lam1,p_lam1,AIC_lam1,dAIC_used,n_total_species)]

# add coordinates (pull from original beta table)
main_out[, chr := beta_ori$chr[match(gene, paste(beta_ori$chr, beta_ori$pos, sep=":"))]]
main_out[, pos := beta_ori$pos[match(gene, paste(beta_ori$chr, beta_ori$pos, sep=":"))]]

# filter & order using --q
main_out <- main_out[!is.na(qvalue) & qvalue < opt$q]
data.table::setorder(main_out, -logP)

main_path <- file.path(opt$outdir, sprintf("correlations_by_gene__%s.tsv", opt$branch))
excl_path <- file.path(opt$outdir, sprintf("correlations_by_gene__%s__EXCLUDED.tsv", opt$branch))
fwrite(main_out, main_path, sep="\t")
excl <- valid_counts[n_valid != n_total_species, .(genus,gene)]
fwrite(excl, excl_path, sep="\t")
logf("Wrote main results: %s (%d rows)", main_path, nrow(main_out))
logf("Wrote excluded groups (not full-case): %s (%d rows)", excl_path, nrow(excl))

## ------------------------- SIMPLE SCATTERPLOTS (p-value only) ----------------
if (nrow(main_out)) {
  logf("Creating scatterplots for %d significant groups (q < %.2f)", nrow(main_out), opt$q)
  idx_rows <- list()
  for (k in seq_len(nrow(main_out))) {
    g  <- main_out$genus[k]
    ge <- main_out$gene[k]
    pv <- main_out$p_used[k]
    # data for this (genus,gene)
    d0 <- res[list(g, ge)][, .(species, beta, h2)]
    trp <- transform_pair(d0$beta, d0$h2, opt$transform)
    d0$beta_t <- trp$beta
    d0$h2_t   <- trp$h2
    # assemble plot
    ttl <- sprintf("%s | %s", g, ge)
    sub <- sprintf("PGLS p = %.3g", pv)
    pl <- ggplot(d0, aes(x=h2_t, y=beta_t)) +
      geom_point(alpha=0.85, size=1.8) +
      labs(title=ttl, subtitle=sub,
           x=sprintf("Heritability (transform=%s)", opt$transform),
           y=sprintf("Variant beta (transform=%s)", opt$transform)) +
      theme_minimal()
    # save
    fn <- file.path(plots_dir,
                    sprintf("scatter__%s__%s__p%s.png",
                            sanitize(g), sanitize(ge),
                            gsub("\\+", "p", gsub("\\-","m", format(pv, digits=3)))))
    try(ggsave(filename=fn, plot=pl, width=5, height=4, dpi=300), silent=TRUE)
    idx_rows[[length(idx_rows)+1]] <- data.table(genus=g, gene=ge, p_used=pv, file=fn)
    if (LOG_EVERY > 0 && (k %% LOG_EVERY == 0L)) logf("... saved %d / %d plots", k, nrow(main_out))
  }
  idx <- rbindlist(idx_rows, fill=TRUE)
  fwrite(idx, file.path(plots_dir, "scatter_index.tsv"), sep="\t")
  logf("Scatterplots saved to: %s  (n=%d)", plots_dir, nrow(main_out))
}

if (length(fail_log)) {
  fail_dt <- data.table::rbindlist(fail_log, use.names=TRUE, fill=TRUE)
  fail_path <- file.path(opt$outdir, sprintf("pgls_failures__%s.tsv", opt$branch))
  fwrite(fail_dt, fail_path, sep="\t")
  tally <- fail_dt[, .N, by=reason][order(-N)]
  logf("Failure breakdown: %s", paste(sprintf("%s=%d", tally$reason, tally$N), collapse=" | "))
  logf("Wrote failure details: %s", fail_path)
}

logf("DONE")
