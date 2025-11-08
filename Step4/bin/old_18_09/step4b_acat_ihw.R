#!/usr/bin/env Rscript
# Step 4B: IHW (primary on ACAT) + weighted BH / BH, PLUS a sensitivity branch:
# IHW/BH directly on GENE×GENUS set tests from Phase 3B.
# Keeps gz-safe reading, species exclusions parity, and weight safeguards.

suppressPackageStartupMessages({
  req <- c("optparse","data.table","IHW")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing R packages: ", paste(miss, collapse=", "))
  library(optparse); library(data.table); library(IHW)
})

LOG <- function(...) cat("[step4B]", ..., "\n", sep="")

## -------------------- CLI --------------------
opt_list <- list(
  make_option("--acat_h2",           type="character", help="acat_genus_h2.tsv"),
  make_option("--acat_eq",           type="character", help="acat_genus_eq.tsv"),
  make_option("--assoc_species_dir", type="character", help="assoc_set_species_full dir"),
  make_option("--assoc_genus_dir",   type="character", help="assoc_set_genus_full dir (Phase 3B)"),
  make_option("--herit_rda",         type="character", help="heritability.RData"),
  make_option("--alpha",             type="double",    help="FDR target (e.g., 0.10)"),
  make_option("--out_sel_h2",        type="character", default="selected_genes_h2.tsv"),
  make_option("--out_sel_eq",        type="character", default="selected_genes_eq.tsv"),
  make_option("--make_plots",        type="character", default="TRUE",
              help="Diagnostics [TRUE/FALSE] (default TRUE)"),
  make_option("--plots_prefix",      type="character", default="diag_acat",
              help="Plot prefix (default 'diag_acat')")
)
opt <- parse_args(OptionParser(option_list = opt_list))
required <- c("acat_h2","acat_eq","assoc_species_dir","assoc_genus_dir","herit_rda","alpha")
miss <- required[!nzchar(vapply(required, function(k) as.character(opt[[k]]), ""))]
if (length(miss)) stop("Missing required options: ", paste(miss, collapse=", "))
make_plots <- toupper(opt$make_plots) %in% c("TRUE","T","1","YES","Y")
alpha <- as.numeric(opt$alpha)

## -------------------- Load ACAT tables --------------------
acat_h2 <- data.table::fread(opt$acat_h2)
acat_eq <- data.table::fread(opt$acat_eq)
stopifnot(all(c("gene","genus","p") %in% names(acat_h2)),
          all(c("gene","genus","p") %in% names(acat_eq)))
setDT(acat_h2); setDT(acat_eq)
acat_h2[, p := as.numeric(p)]; acat_eq[, p := as.numeric(p)]
acat_h2 <- acat_h2[is.finite(p) & p > 0 & p <= 1]
acat_eq <- acat_eq[is.finite(p) & p > 0 & p <= 1]

## -------------------- Heritability loader (all_VCs_full preferred) ----------
load(opt$herit_rda)
objs <- mget(ls(), inherits = TRUE)
h2_from_objs <- function(objs) {
  if ("all_VCs_full" %in% names(objs)) {
    dt <- as.data.table(objs[["all_VCs_full"]])
    if (all(c("trait","var_Ad") %in% names(dt))) {
      dt <- dt[, .(trait = as.character(trait),
                   h2    = pmax(as.numeric(var_Ad), 0))]
      return(dt[!is.na(trait) & is.finite(h2)])
    }
  }
  for (nm in names(objs)) {
    x <- objs[[nm]]
    if (!is.data.frame(x) && !data.table::is.data.table(x)) next
    x <- as.data.table(x)
    idc <- intersect(names(x), c("trait","species","name","id","taxon"))
    hc  <- intersect(names(x), c("h2","h2_add","heritability","H2","h2_obs","h2_adj"))
    if (length(idc) && length(hc)) {
      dt <- data.table(trait = as.character(x[[idc[1]]]),
                       h2    = pmax(as.numeric(x[[hc[1]]]), 0))
      dt <- dt[!is.na(trait) & is.finite(h2)]
      if (nrow(dt)) return(dt)
    }
  }
  NULL
}
h2map <- h2_from_objs(objs)
if (is.null(h2map)) stop("Could not identify a (trait,h2) table in: ", opt$herit_rda)

trait_to_genus <- function(tr){
  m <- regexec("^s__([^_]+)_", tr); mt <- regmatches(tr, m)
  if (length(mt) && length(mt[[1]])==2) mt[[1]][2] else NA_character_
}
h2map <- h2map[grepl("^s__", trait)]
h2map[, genus := vapply(trait, trait_to_genus, character(1))]
h2map <- h2map[!is.na(genus)]
h2map[, h2_sum := sum(h2), by = genus]
h2map[, w0 := ifelse(h2_sum > 0, h2 / h2_sum, 1/.N), by = genus]
h2w <- h2map[, .(trait, genus, w0)]

## -------------------- Species set files → build covariate (median n_var) ----
by_trait_species <- file.path(opt$assoc_species_dir, "by_gene_by_trait")
if (!dir.exists(by_trait_species)) stop("Missing dir: ", by_trait_species)
fls_sp <- list.files(by_trait_species, pattern="^s__.*\\.geneset\\.tsv(\\.gz)?$", full.names=TRUE)
fls_sp <- fls_sp[!grepl("/s__Bacteroides_F_", fls_sp)]
if (!length(fls_sp)) stop("No species gene-set files in: ", by_trait_species)

species_trait <- sub("\\.(geneset\\.tsv(\\.gz)?)$","", basename(fls_sp))  # s__Genus_species
genus_sp <- vapply(species_trait, function(tr) sub("^s__([^_]+).*","\\1", tr), character(1), USE.NAMES=FALSE)
keep <- !is.na(genus_sp); fls_sp <- fls_sp[keep]; species_trait <- species_trait[keep]; genus_sp <- genus_sp[keep]

.safe_read <- function(f, select=NULL){
  dt <- try(data.table::fread(f, sep="\t", header=TRUE, data.table=TRUE,
                              showProgress=FALSE, select=select), silent=TRUE)
  if (inherits(dt, "try-error") || !is.data.frame(dt)) {
    con <- gzfile(f, open="rt")
    dt2 <- try(read.table(con, header=TRUE, sep="\t", quote="", comment.char="",
                          check.names=FALSE, stringsAsFactors=FALSE), silent=TRUE)
    close(con)
    if (inherits(dt2, "try-error") || !is.data.frame(dt2)) return(NULL)
    if (!is.null(select)) dt2 <- dt2[, intersect(select, colnames(dt2)), drop=FALSE]
    data.table::setDT(dt2); return(dt2)
  }
  dt
}

get_gene_nvar <- function(f){
  dt <- .safe_read(f, select = c("GENE_ID","n_var","n.site"))
  if (is.null(dt)) return(NULL)
  # unify set-size column name
  if (!"n_var" %in% names(dt) && "n.site" %in% names(dt)) set(dt, j="n_var", value=dt[["n.site"]])
  if (!all(c("GENE_ID","n_var") %in% names(dt))) return(NULL)
  dt[, n_var := suppressWarnings(as.numeric(n_var))]
  dt <- dt[is.finite(n_var) & n_var > 0]
  if (!nrow(dt)) return(NULL)
  dt[, .(GENE_ID = as.character(GENE_ID), n_var)]
}
nv_by_file <- lapply(fls_sp, get_gene_nvar)

cov_nvar_dt <- rbindlist(lapply(unique(genus_sp), function(gn){
  idxs <- which(genus_sp == gn)
  if (!length(idxs)) return(NULL)
  nv <- rbindlist(nv_by_file[idxs], use.names = TRUE, fill = TRUE)
  if (!is.data.frame(nv) || !nrow(nv)) return(NULL)
  nv[, .(gene = GENE_ID, cov_nvar = median(n_var, na.rm = TRUE)), by = GENE_ID][, genus := gn][]
}), use.names = TRUE, fill = TRUE)

# attach covariate to ACAT tables; jitter to avoid zero/constant issues
eps <- 1e-8
acat_h2 <- merge(acat_h2, cov_nvar_dt, by = c("gene","genus"), all.x = TRUE)
acat_eq <- merge(acat_eq, cov_nvar_dt, by = c("gene","genus"), all.x = TRUE)
acat_h2[is.na(cov_nvar), cov_nvar := 0]
acat_eq[is.na(cov_nvar), cov_nvar := 0]
acat_h2[, cov_nvar_j := pmax(cov_nvar + eps, eps)]
acat_eq[, cov_nvar_j := pmax(cov_nvar + eps, eps)]

# Weighted BH weights from this covariate (mean 1 within genus)
acat_h2[, w_wbh := cov_nvar_j / mean(cov_nvar_j), by = genus]
acat_eq[, w_wbh := cov_nvar_j / mean(cov_nvar_j), by = genus]
for (DT in list(acat_h2, acat_eq)) {
  DT[!is.finite(w_wbh) | w_wbh <= 0, w_wbh := 1]
  DT[, w_wbh := w_wbh / mean(w_wbh), by = genus]
}

## -------------------- FDR helpers --------------------
ihw_select <- function(dt, cov_col){
  rbindlist(lapply(split(dt, dt$genus), function(x){
    covv <- x[[cov_col]]
    if (!any(is.finite(covv)) || isTRUE(all(abs(covv - covv[1]) < 1e-12))) {
      x$q_bh <- p.adjust(x$p, method="BH"); 
      return(x[q_bh <= alpha, .(gene, genus, p, q_ihw = q_bh, cov = covv[q_bh <= alpha])])
    }
    res <- try(ihw(x$p ~ covv, alpha = alpha), silent=TRUE)
    if (inherits(res, "try-error")) {
      x$q_bh <- p.adjust(x$p, method="BH")
      return(x[q_bh <= alpha, .(gene, genus, p, q_ihw = q_bh, cov = covv[q_bh <= alpha])])
    } else {
      q <- adj_pvalues(res)
      return(data.table(gene = x$gene, genus = x$genus, p = x$p,
                        q_ihw = as.numeric(q), cov = covv)[q_ihw <= alpha])
    }
  }), use.names=TRUE)
}
wbh_select <- function(dt){
  rbindlist(lapply(split(dt, dt$genus), function(x){
    w <- x$w_wbh; w[!is.finite(w) | w <= 0] <- 1; w <- w / mean(w)
    pprime <- x$p / w
    x$q_wbh <- p.adjust(pprime, method="BH")
    x[q_wbh <= alpha, .(gene, genus, p, q_wbh)]
  }), use.names=TRUE)
}
bh_select <- function(dt){
  rbindlist(lapply(split(dt, dt$genus), function(x){
    x$q_bh <- p.adjust(x$p, method="BH")
    x[q_bh <= alpha, .(gene, genus, p, q_bh)]
  }), use.names=TRUE)
}

## -------------------- PRIMARY: IHW on ACAT(h2/eq) ---------------------------
LOG("Running IHW (primary) on ACAT(h2)")
sel_h2 <- ihw_select(acat_h2, "cov_nvar_j")
LOG("Running IHW (primary) on ACAT(eq)")
sel_eq <- ihw_select(acat_eq, "cov_nvar_j")

# sensitivities on ACAT
LOG("Running weighted BH (sensitivity) on ACAT(h2)")
sel_h2_wbh <- wbh_select(acat_h2)
LOG("Running vanilla BH (sensitivity) on ACAT(h2)")
sel_h2_bh  <- bh_select(acat_h2)

## -------------------- NEW SENSITIVITY: IHW/BH on GENE×GENUS (Phase 3B) ------
by_trait_genus <- file.path(opt$assoc_genus_dir, "by_gene_by_trait")
if (!dir.exists(by_trait_genus)) stop("Missing dir: ", by_trait_genus)

fls_gn <- list.files(by_trait_genus, pattern="^g__.*\\.geneset\\.tsv(\\.gz)?$", full.names=TRUE)
if (length(fls_gn)) {
  # read all genus gene-set tables (expect columns: GENE_ID, p, n_var or n.site)
  read_genus <- function(f){
    dt <- .safe_read(f, select = c("GENE_ID","p","n_var","n.site"))
    if (is.null(dt)) return(NULL)
    if (!"p" %in% names(dt)) {
      # ultra-robust fallback if a p-like col appears with different name
      pcol <- intersect(tolower(names(dt)),
                        c("p","p.value","pval","pval_smmat","pval_skato","min.pval","score.pval"))
      if (length(pcol)) setnames(dt, names(dt)[match(pcol[1], tolower(names(dt)))], "p")
    }
    if (!"n_var" %in% names(dt) && "n.site" %in% names(dt)) set(dt, j="n_var", value=dt[["n.site"]])
    if (!all(c("GENE_ID","p") %in% names(dt))) return(NULL)
    dt[, p := suppressWarnings(as.numeric(p))]
    if ("n_var" %in% names(dt)) dt[, n_var := suppressWarnings(as.numeric(n_var))] else dt[, n_var := NA_real_]
    gn <- sub("^g__([^\\.]+).*","\\1", basename(f))
    dt <- dt[is.finite(p) & p > 0 & p <= 1]
    if (!nrow(dt)) return(NULL)
    dt[, .(gene = as.character(GENE_ID), genus = gn, p, n_var)]
  }
  G <- rbindlist(lapply(fls_gn, read_genus), use.names=TRUE, fill=TRUE)
  if (nrow(G)) {
    # covariate = n_var (set size); jitter & sanitize
    eps2 <- 1e-8
    G[!is.finite(n_var) | n_var <= 0, n_var := 0]
    G[, cov_nvar_j := pmax(n_var + eps2, eps2)]
    # IHW per genus
    LOG("Running IHW (sensitivity) on GENE×GENUS set tests")
    sel_genus_ihw <- ihw_select(G, "cov_nvar_j")
    # BH per genus
    LOG("Running BH (sensitivity) on GENE×GENUS set tests")
    sel_genus_bh  <- bh_select(G)
    # write
    fwrite(sel_genus_ihw, file = "selected_genes_genus_ihw.tsv", sep = "\t")
    fwrite(sel_genus_bh,  file = "selected_genes_genus_bh.tsv",  sep = "\t")
  } else {
    LOG("GENE×GENUS tables found, but no valid rows parsed; skipping genus IHW/BH outputs.")
    fwrite(data.table(), file = "selected_genes_genus_ihw.tsv", sep = "\t")
    fwrite(data.table(), file = "selected_genes_genus_bh.tsv",  sep = "\t")
  }
} else {
  LOG("No GENE×GENUS files found; skipping genus IHW/BH outputs.")
  fwrite(data.table(), file = "selected_genes_genus_ihw.tsv", sep = "\t")
  fwrite(data.table(), file = "selected_genes_genus_bh.tsv",  sep = "\t")
}

## -------------------- Write ACAT outputs --------------------
fwrite(sel_h2, file = opt$out_sel_h2, sep = "\t")
fwrite(sel_eq, file = opt$out_sel_eq, sep = "\t")
fwrite(sel_h2_wbh, file = "selected_genes_h2_wbh.tsv", sep = "\t")
fwrite(sel_h2_bh,  file = "selected_genes_h2_bh.tsv",  sep = "\t")
LOG("Selections written: ", opt$out_sel_h2, " and ", opt$out_sel_eq)

## -------------------- (plots unchanged; omitted for brevity) -----------------
if (make_plots) {
  has_gg <- requireNamespace("ggplot2", quietly = TRUE)
  if (!has_gg) LOG("ggplot2 not installed; skipping plots.")
  # ... keep your existing plotting block, or leave as-is ...
}

LOG("Done.")
