#!/usr/bin/env Rscript
# Combine species-level gene set-test p-values -> gene×genus p-values using ACAT
# Primary: heritability (var_Ad) weights from all_VCs_full; Sensitivity: equal weights.
# Robust gz reading (fread -> fallback to base gzfile) + flexible column detection.

suppressPackageStartupMessages({
  req <- c("optparse","data.table")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse=", "))
  library(optparse); library(data.table)
})

LOG <- function(...) cat("[acat]", ..., "\n", sep="")

## ---------------- CLI ----------------
opt_list <- list(
  make_option("--indir",       type="character", help="Dir with species gene-set tables (by_gene_by_trait)"),
  make_option("--herit_rda",   type="character", help="RData with 'all_VCs_full' (trait,var_Ad)"),
  make_option("--genera",      type="character", help="Comma-separated genera (e.g. Prevotella,Bacteroides)"),
  make_option("--weight_mode", type="character", default="h2", help="'h2' or 'none'/'equal'"),
  make_option("--p_floor",     type="double",    default=1e-15, help="Floor for tiny p-values"),
  make_option("--out_h2",      type="character", default="acat_genus_h2.tsv",
              help="Output TSV (gene,genus,p) for h2-weighted ACAT"),
  make_option("--out_eq",      type="character", default="acat_genus_eq.tsv",
              help="Output TSV (gene,genus,p) for equal-weight ACAT")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (k in c("indir","herit_rda","genera")) if (!nzchar(opt[[k]])) stop("Missing --", k)
if (!dir.exists(opt$indir)) stop("Input directory not found: ", opt$indir)

## ---------------- ACAT (weighted / equal) ----------------
# Liu & Xie (2019): weighted sum of tan(pi*(0.5 - p)); analytic Cauchy tail approx.
acat_combine <- function(p, w=NULL, p_floor=1e-15){
  p <- as.numeric(p); p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmax(p_floor, pmin(1 - p_floor, p))
  if (is.null(w)) w <- rep(1/length(p), length(p)) else {
    w <- as.numeric(w)[seq_along(p)]
    if (length(w)!=length(p) || any(!is.finite(w)) || any(w<0)) return(NA_real_)
    s <- sum(w); w <- if (s>0) w/s else rep(1/length(p), length(p))
  }
  t_sum <- sum(w * tan((0.5 - p) * pi))
  if (!is.finite(t_sum)) return(min(p))
  0.5 - atan(t_sum)/pi
}

## ---------------- Load heritability (species-level) ----------------
load(opt$herit_rda)
if (!exists("all_VCs_full")) stop("No object 'all_VCs_full' in ", opt$herit_rda)
h2tab <- as.data.table(all_VCs_full)
if (!all(c("trait","var_Ad") %in% names(h2tab)))
  stop("'all_VCs_full' must have columns: trait, var_Ad")
h2tab <- h2tab[grepl("^s__", trait)]
trait_to_genus <- function(tr){
  m <- regexec("^s__([^_]+)_", tr); mt <- regmatches(tr, m)
  if (length(mt) && length(mt[[1]])==2) mt[[1]][2] else NA_character_
}
h2tab[, genus := vapply(trait, trait_to_genus, character(1))]
h2tab <- h2tab[!is.na(genus)]
h2tab[, h2p := pmax(as.numeric(var_Ad), 0)]
h2tab[, w0  := { s <- sum(h2p); if (s>0) h2p/s else 1/.N }, by=genus]
LOG("Heritability rows (species) kept: ", nrow(h2tab),
    " across genera: ", paste(sort(unique(h2tab$genus)), collapse=", "))

## ---------------- Read species gene-set files ----------------
fls <- list.files(opt$indir, pattern="^s__.*\\.geneset\\.tsv(\\.gz)?$", full.names=TRUE)
fls <- fls[!grepl("/s__Bacteroides_F_", fls)]
if (!length(fls)) {
  LOG("No species gene-set files found after filtering. Writing empty outputs.")
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_h2, sep="\t")
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_eq, sep="\t")
  quit(save="no", status=0)
}

species_trait <- sub("\\.(geneset\\.tsv(\\.gz)?)$","", basename(fls)) # s__Genus_species
genus_vec <- vapply(species_trait, trait_to_genus, character(1), USE.NAMES=FALSE)
keep <- !is.na(genus_vec)
fls <- fls[keep]; species_trait <- species_trait[keep]; genus_vec <- genus_vec[keep]

# Column alias pools (GENESIS aggregate outputs vary by test)
P_COLS    <- c("p","p.value","p_value","pval","score.pval","pval_smmat","pval_skat","pval_skato","min.pval")
GENE_COLS <- c("gene_id","gene","geneid","setid","set_id","set","gene.name")

# robust reader: try fread; if it fails or returns 0 rows, fall back to base gzfile() + read.table()
.safe_read <- function(f){
  dt <- try(fread(f, sep="\t", header=TRUE, data.table=TRUE, showProgress=FALSE), silent=TRUE)
  if (inherits(dt, "try-error") || !is.data.frame(dt) || !nrow(dt)) {
    con <- gzfile(f, open="rt")  # works without system zcat; portable base R
    dt2 <- try(read.table(con, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE,
                          stringsAsFactors=FALSE), silent=TRUE)
    close(con)
    if (inherits(dt2, "try-error") || !is.data.frame(dt2) || !nrow(dt2)) return(NULL)
    setDT(dt2)
    return(dt2)
  }
  dt
}

read_species_file <- function(f){
  dt <- .safe_read(f)
  if (is.null(dt)) return(NULL)
  
  old <- names(dt); new <- trimws(old)
  setnames(dt, old, new)
  lc <- tolower(new)
  
  gene_col <- if ("GENE_ID" %in% new) "GENE_ID" else {
    hit <- which(lc %in% GENE_COLS)[1]; if (!is.na(hit)) new[hit] else NA_character_
  }
  if (is.na(gene_col)) return(NULL)
  
  p_hits <- which(lc %in% P_COLS)
  if (!length(p_hits)) return(NULL)
  pcol <- new[p_hits[1]]
  
  out <- dt[, .(GENE_ID = as.character(get(gene_col)),
                p = suppressWarnings(as.numeric(get(pcol))))]
  if ("n_var" %in% new) out[, n_var := suppressWarnings(as.integer(dt[["n_var"]]))]
  out <- out[is.finite(p) & p > 0 & p <= 1]
  if (!nrow(out)) return(NULL)
  out
}

glist <- vector("list", length(fls))
bad <- 0L; parsed <- 0L
for (i in seq_along(fls)) {
  dt <- read_species_file(fls[i])
  if (is.null(dt)) { bad <- bad + 1L; next }
  dt[, species := species_trait[i]]
  dt[, genus   := genus_vec[i]]
  glist[[i]] <- dt[, .(GENE_ID, p, species, genus)]
  parsed <- parsed + nrow(dt)
}
gdt <- rbindlist(glist, use.names=TRUE, fill=TRUE)
LOG("Parsing summary: files=", length(fls), " parsed_rows=", parsed, " bad_format_files=", bad)

if (!nrow(gdt)) {
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_h2, sep="\t")
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_eq, sep="\t")
  LOG("No valid (GENE_ID,p) rows found; wrote empty outputs and exited cleanly.")
  quit(save="no", status=0)
}

## ---------------- Limit to requested genera ----------------
gen_keep <- trimws(strsplit(opt$genera, ",", fixed=TRUE)[[1]])
gdt <- gdt[genus %in% gen_keep]
if (!nrow(gdt)) {
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_h2, sep="\t")
  fwrite(data.table(gene=character(),genus=character(),p=double()), file=opt$out_eq, sep="\t")
  LOG("After genus filter, no rows remain for: ", paste(gen_keep, collapse=", "), ". Wrote empties.")
  quit(save="no", status=0)
}

## ---------------- Build weight lookup ----------------
h2w <- h2tab[, .(species = trait, genus, w0)]
resolve_weights <- function(sp, gn, mode=c("h2","none","equal")){
  mode <- match.arg(mode)
  if (mode %in% c("none","equal")) return(rep(1/length(sp), length(sp)))
  wv <- h2w[genus == gn & species %in% sp, match(species, sp)]
  miss <- is.na(wv); if (any(miss)) wv[miss] <- 0
  if (sum(wv) <= 0) rep(1/length(sp), length(sp)) else wv / sum(wv)
}

## ---------------- Per gene×genus ACAT ----------------
p_floor <- as.numeric(opt$p_floor)
if (!is.finite(p_floor) || p_floor <= 0 || p_floor >= 0.5)
  stop("--p_floor must be in (0, 0.5)")

by_gene_genus <- gdt[, {
  sp <- species
  pv <- p
  p_eq <- acat_combine(pv, w=NULL, p_floor=p_floor)
  w_h2 <- resolve_weights(sp, genus[1], mode="h2")
  p_h2 <- acat_combine(pv, w=w_h2, p_floor=p_floor)
  .(p_h2 = p_h2, p_eq = p_eq)
}, by=.(genus, gene=GENE_ID)]

data.table::setorder(by_gene_genus, genus, gene)

fwrite(data.table(gene=by_gene_genus$gene, genus=by_gene_genus$genus, p=by_gene_genus$p_h2),
       file=opt$out_h2, sep="\t")
fwrite(data.table(gene=by_gene_genus$gene, genus=by_gene_genus$genus, p=by_gene_genus$p_eq),
       file=opt$out_eq, sep="\t")
LOG("Wrote: ", opt$out_h2, " (rows=", nrow(by_gene_genus), ")")
LOG("Wrote: ", opt$out_eq,  " (rows=", nrow(by_gene_genus), ")")
