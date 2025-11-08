#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# GENUS-LEVEL GWAS with GENESIS + LOCO + SRMs
# Robust outcome handling:
#   - Build pd_all from FULL GDS sample list (no reliance on filtered sampleData)
#   - Subset to pd_sel that matches sel_int order exactly
#   - Fit null model on pd_sel (data.frame), then run assoc on SeqVarData iterator
# References:
#   - SeqVarTools sampleData/GDS alignment & filtering behavior  [docs]        (1)
#   - GENESIS fitNullModel cov.mat (list of covariance matrices) [manual/vig]   (2)
# (1) https://rdrr.io/bioc/SeqVarTools/man/SeqVarData-class.html ; PDF manual   (turn0search1, turn0search8)
# (2) GENESIS manual/vignettes: fitNullModel, assocTestSingle                    (turn0search3, turn0search4)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(Biobase)
  library(SeqArray); library(SeqVarTools); library(GENESIS)
})

VAR_BLOCK <- as.integer(Sys.getenv("GENESIS_VARBLOCK", "2000"))

opt <- parse_args(OptionParser(option_list=list(
  make_option("--gds", type="character"),
  make_option("--scan_lookup", type="character"),
  make_option("--resid_rda", type="character"),
  make_option("--metadata", type="character"),
  make_option("--grm_rdata", type="character"),
  make_option("--srm_dam", type="character"),
  make_option("--srm_cage", type="character"),
  make_option("--genera", type="character", default="Prevotella,Bacteroides"),
  make_option("--genus_prefix", type="character", default="g__"),
  make_option("--outdir", type="character", default="assoc_genus"),
  # Variant QC
  make_option("--maf_min",     type="double",  default=NA),
  make_option("--mac_min",     type="integer", default=NA),
  make_option("--miss_max",    type="double",  default=NA),
  make_option("--hwe_p_min",   type="double",  default=NA),
  make_option("--hwe_maf_min", type="double",  default=NA),
  # LOCO
  make_option("--loco_grm_rds", type="character", default=NA)
)))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
log <- function(...) cat("[genus_screen]", ..., "\n")

# ---------- helpers ----------
read_grm_prefix <- function(prefix, sel_ids_char) {
  idp <- paste0(prefix, ".grm.id"); bin <- paste0(prefix, ".grm.bin")
  if (!file.exists(idp) || !file.exists(bin)) stop("Missing GRM files for prefix: ", prefix)
  ids <- read.table(idp, stringsAsFactors=FALSE)
  stopifnot(ncol(ids) >= 2)
  iid <- as.character(ids[[2]]); n <- length(iid)
  con <- file(bin, "rb"); on.exit(close(con), add=TRUE)
  lt  <- readBin(con, what=numeric(), n=n*(n+1)/2, size=4)
  if (length(lt) != n*(n+1)/2) stop("GRM read size mismatch for ", prefix)
  M <- matrix(0, n, n); k <- 1L
  for (i in 1:n) for (j in 1:i) { M[i,j] <- lt[k]; M[j,i] <- lt[k]; k <- k+1L }
  rownames(M) <- iid; colnames(M) <- iid
  if (!all(sel_ids_char %in% iid)) {
    miss <- sel_ids_char[!(sel_ids_char %in% iid)]
    stop("LOCO GRM ", prefix, " missing IDs; e.g. ", paste(head(miss,5), collapse=", "))
  }
  M[sel_ids_char, sel_ids_char, drop=FALSE]
}

# ---------- residuals (GENUS layer) ----------
genus <- strsplit(opt$genera, ",", fixed=TRUE)[[1]]
log("Loading residuals:", opt$resid_rda)
load(opt$resid_rda)
if (!exists("residuals_qned_counts_objs"))
  stop("Expected object 'residuals_qned_counts_objs' in ", opt$resid_rda)

residuals <- t(residuals_qned_counts_objs[[length(residuals_qned_counts_objs)-1]])
rownames(residuals) <- sapply(strsplit(rownames(residuals), "_"), `[`, 1)
log("#samples=", nrow(residuals), " #traits=", ncol(residuals))

want <- paste0(opt$genus_prefix, genus)
want <- want[!grepl("^g__Bacteroides_F", want)]
keep <- intersect(colnames(residuals), want)
if (!length(keep)) stop("No genus columns found among: ", paste(want, collapse=", "))
resid_g <- residuals[, keep, drop=FALSE]
log("Genus traits: ", paste(colnames(resid_g), collapse=", "))

# ---------- metadata / GRM / SRMs ----------
meta <- fread(opt$metadata)
id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
if (is.na(id_col)) stop("Metadata must have ID column among sample/Sample/RFID/ID/iid")
meta <- meta[!duplicated(meta[[id_col]])]
load(opt$grm_rdata); if (!exists("grm")) stop("No 'grm' in grm_rdata")
srm_dam  <- readRDS(opt$srm_dam)
srm_cage <- readRDS(opt$srm_cage)

# ---------- GDS + mappings ----------
gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
gds_all <- seqGetData(gds, "sample.id")                        # FULL sample list in GDS
lookup  <- fread(opt$scan_lookup)
stopifnot(all(c("scan_id","sample_label") %in% names(lookup)))
label2int <- setNames(lookup$scan_id,  lookup$sample_label)
int2label <- setNames(lookup$sample_label, lookup$scan_id)

sids_all <- rownames(resid_g)
lab <- Reduce(intersect, list(sids_all, meta[[id_col]], rownames(grm),
                              rownames(srm_dam), rownames(srm_cage),
                              int2label[as.character(gds_all)]))
if (!length(lab)) stop("No overlap across inputs")

sel_int  <- as.integer(label2int[lab]); if (anyNA(sel_int)) stop("Some labels didn’t map to scan_id")
sel_char <- as.character(sel_int)

# reorder phenos & covariances by label order; rename to scan IDs
resid_g  <- resid_g[lab,,drop=FALSE]
grm      <- grm[lab, lab];        rownames(grm)      <- colnames(grm)      <- sel_char
srm_dam  <- srm_dam[lab, lab];    rownames(srm_dam)  <- colnames(srm_dam)  <- sel_char
srm_cage <- srm_cage[lab, lab];   rownames(srm_cage) <- colnames(srm_cage) <- sel_char

# ---------- SeqVarData for iterators / QC only ----------
annot0 <- AnnotatedDataFrame(data.frame(sample.id=gds_all, row.names=as.character(gds_all)))
seqData <- SeqVarData(gds, sampleData=annot0)  # sampleData must match GDS sample.id (unfiltered) :contentReference[oaicite:3]{index=3}
seqSetFilter(seqData, sample.id=sel_int, action="set", verbose=FALSE)

nsamp <- length(seqGetData(seqData, "sample.id"))
nvar  <- length(seqGetData(seqData, "variant.id"))
message("# of selected samples: ", nsamp)
message("# of selected variants: ", format(nvar, big.mark=","))

# ---------- Variant QC ----------
apply_qc <- !(is.na(opt$maf_min) & is.na(opt$mac_min) & is.na(opt$miss_max) &
                is.na(opt$hwe_p_min) & is.na(opt$hwe_maf_min))
keep_var <- rep(TRUE, nvar)
if (apply_qc) {
  miss  <- missingGenotypeRate(seqData, margin="by.variant")
  refAF <- alleleFrequency(seqData, n=0)
  maf   <- pmin(refAF, 1 - refAF)
  n_called <- round((1 - miss) * nsamp)
  mac      <- round(maf * 2L * n_called)
  miss[is.na(miss)] <- 1; maf[is.na(maf)] <- 0; mac[is.na(mac)] <- 0L
  
  hwe_p <- rep(1, length(maf)); n_hwe_tested <- 0L; n_hwe_fail <- 0L
  if (!is.na(opt$hwe_p_min)) {
    hw <- try(hwe(seqData, permute=FALSE), silent=TRUE)
    if (!inherits(hw, "try-error") && is.data.frame(hw) && "p" %in% names(hw)) {
      hwe_p <- hw$p; hwe_p[is.na(hwe_p)] <- 0
      maf_cut <- ifelse(is.na(opt$hwe_maf_min), 0.05, opt$hwe_maf_min)
      idx_hwe <- which(maf >= maf_cut)
      n_hwe_tested <- length(idx_hwe); n_hwe_fail <- sum(hwe_p[idx_hwe] < opt$hwe_p_min)
    } else { log("WARNING: HWE computation failed; skipping HWE filter"); opt$hwe_p_min <- NA }
  }
  
  keep_var <- rep(TRUE, length(miss))
  if (!is.na(opt$miss_max)) keep_var <- keep_var & (miss <= opt$miss_max)
  if (!is.na(opt$maf_min))  keep_var <- keep_var & (maf  >= opt$maf_min)
  if (!is.na(opt$mac_min))  keep_var <- keep_var & (mac  >= opt$mac_min)
  if (!is.na(opt$hwe_p_min)) {
    maf_cut <- ifelse(is.na(opt$hwe_maf_min), 0.05, opt$hwe_maf_min)
    keep_var <- keep_var & ( (maf < maf_cut) | (hwe_p >= opt$hwe_p_min) )
  }
  
  qc_dt <- data.table(
    metric=c("nvar_total","nvar_keep","prop_keep","miss_max","maf_min","mac_min",
             "hwe_p_min","hwe_maf_min","n_hwe_tested","n_hwe_fail","n_samples"),
    value = c(nvar, sum(keep_var), ifelse(nvar>0, sum(keep_var)/nvar, NA),
              opt$miss_max, opt$maf_min, opt$mac_min,
              opt$hwe_p_min, opt$hwe_maf_min, n_hwe_tested, n_hwe_fail, nsamp)
  )
  fwrite(qc_dt, file.path(opt$outdir, "genesis_variant_qc.tsv"), sep="\t")
}

# ---------- LOCO GRM preload ----------
use_loco <- !is.na(opt$loco_grm_rds) && nzchar(opt$loco_grm_rds)
loco_mats <- NULL
if (use_loco) {
  loco_map <- readRDS(opt$loco_grm_rds) # named vector: chr -> prefix
  if (is.null(names(loco_map)) || !length(loco_map)) stop("loco_grm_rds empty/malformed")
  loco_mats <- lapply(names(loco_map), function(chr) read_grm_prefix(loco_map[[chr]], sel_char))
  names(loco_mats) <- names(loco_map)
}

# ---------- Chromosomes present ----------
chr_all    <- seqGetData(seqData, "chromosome")
chr_levels <- sort(unique(chr_all[keep_var]), na.last=NA)

# ---------- Association loop ----------
for (tr in colnames(resid_g)) {
  y <- as.numeric(resid_g[, tr])
  if (all(is.na(y)) || length(unique(na.omit(y))) < 3) { log("Skipping ", tr, " (too few unique)"); next }
  
  # Build outcome on FULL GDS list, then subset to cohort (avoid filtered sampleData trap) :contentReference[oaicite:4]{index=4}
  pd_all <- data.frame(sample.id = gds_all, row.names = as.character(gds_all), stringsAsFactors = FALSE)
  pd_all$y <- NA_real_
  stopifnot(all(sel_char %in% rownames(pd_all)))
  pd_all[sel_char, "y"] <- y
  
  # Subset EXACTLY to cohort order
  pd_sel <- pd_all[sel_char, , drop=FALSE]
  message(sprintf("[DEBUG:%s] n_all=%d n_pd_sel=%d n_y_nonNA=%d",
                  tr, nrow(pd_all), nrow(pd_sel), sum(!is.na(pd_sel$y))))
  stopifnot(identical(rownames(pd_sel), sel_char))
  
  res_parts <- vector("list", length(chr_levels)); names(res_parts) <- as.character(chr_levels); part_i <- 1L
  
  for (chr in chr_levels) {
    sel_chr <- (chr_all == chr) & keep_var
    if (!any(sel_chr)) next
    seqSetFilter(seqData, variant.sel = sel_chr, action="set", verbose=FALSE)
    
    covlist <- list(dam=srm_dam, cage=srm_cage,
                    grm = if (use_loco) loco_mats[[as.character(chr)]] else grm)
    
    message(sprintf("[COV DEBUG chr=%s] dims: grm=%dx%d dam=%dx%d cage=%dx%d",
                    as.character(chr),
                    nrow(covlist$grm),  ncol(covlist$grm),
                    nrow(covlist$dam),  ncol(covlist$dam),
                    nrow(covlist$cage), ncol(covlist$cage)))
    
    # Fit null on pd_sel (data.frame) with matching order to covlist :contentReference[oaicite:5]{index=5}
    nullmod <- fitNullModel(pd_sel,
                            outcome="y",
                            covars=NULL,          # intercept-only fixed effects
                            cov.mat=covlist,      # multiple random effects are allowed
                            family="gaussian",
                            verbose=FALSE)
    
    it <- SeqVarBlockIterator(seqData, variantBlock=VAR_BLOCK, verbose=FALSE)
    assoc <- assocTestSingle(it, nullmod, test="Score", verbose=FALSE)  # :contentReference[oaicite:6]{index=6}
    if (!is.null(assoc) && nrow(assoc) > 0) res_parts[[part_i]] <- as.data.frame(assoc)
    rm(nullmod, it, assoc); gc()
    part_i <- part_i + 1L
  }
  
  res_all <- do.call(rbind, Filter(Negate(is.null), res_parts))
  fn <- file.path(opt$outdir, paste0(tr, "_assoc.tsv.gz"))
  if (!is.null(res_all) && nrow(res_all) > 0) { fwrite(res_all, fn, sep="\t"); log("Wrote ", fn, " [", nrow(res_all), " rows]") }
  else { con <- gzfile(fn, "w"); close(con); log("No associations (empty): ", fn) }
}

# Sidecars for MAGMA
fwrite(data.table(IID=sel_char), file.path(opt$outdir,"genesis_final_scan_ids.txt"), col.names=FALSE)
fwrite(data.table(label=lab),    file.path(opt$outdir,"genesis_final_label_ids.txt"), col.names=FALSE)
log("DONE")
