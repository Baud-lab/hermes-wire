#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(SeqArray); library(GENESIS); library(Biobase)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--gds",         type="character"),
  make_option("--scan_lookup", type="character"),
  make_option("--trait",       type="character"),
  make_option("--resid_rda",   type="character"),
  make_option("--metadata",    type="character"),
  make_option("--grm_rdata",   type="character"),
  make_option("--srm_dam",     type="character"),
  make_option("--srm_cage",    type="character"),
  make_option("--outdir",      type="character", default="cohort_ids")
)))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

cat("[cohort_ids] loading residuals (for ID universe)\n")
load(opt$resid_rda)
if (opt$trait=="Taxa"){
  if (!exists("residuals_qned_counts_objs"))
    stop("Expected 'residuals_qned_counts_objs' in --resid_rda")
  genus_layer_idx <- length(residuals_qned_counts_objs) - 1L
  residuals <- t(residuals_qned_counts_objs[[genus_layer_idx]])
  
} else if (opt$trait=="Guild") {
  residuals <- t(residuals_qned_counts_clusters_objs[[1]])
  
} else if (opt$trait=="Beta") {
  residuals <- t(residuals_qned_counts_beta_objs[[1]])
  
} else if (opt$trait=="Alpha") {
  residuals <- t(residuals_qned_counts_alpha_objs[[1]])
}
rownames(residuals) <- sapply(strsplit(rownames(residuals), "_"), `[`, 1)
sids <- rownames(residuals)

cat("[cohort_ids] loading metadata\n")
meta <- fread(opt$metadata)
id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
if (is.na(id_col)) stop("Metadata needs an ID column among sample/Sample/RFID/ID/iid")
meta <- meta[!duplicated(meta[[id_col]])]

cat("[cohort_ids] loading GRM & SRMs (for ID universe)\n")
load(opt$grm_rdata); if (!exists("grm")) stop("No 'grm' in grm_rdata")
srm_dam  <- readRDS(opt$srm_dam)
srm_cage <- readRDS(opt$srm_cage)

cat("[cohort_ids] opening GDS & scan lookup\n")
gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
lookup <- fread(opt$scan_lookup)
stopifnot(all(c("scan_id","sample_label") %in% names(lookup)))
label2int <- setNames(lookup$scan_id,  lookup$sample_label)
int2label <- setNames(lookup$sample_label, lookup$scan_id)
gds_int_all <- seqGetData(gds, "sample.id")
gds_lab_all <- int2label[as.character(gds_int_all)]

cat("[cohort_ids] intersecting ID sets across inputs\n")
lab <- Reduce(intersect, list(sids, meta[[id_col]], rownames(grm),
                              rownames(srm_dam), rownames(srm_cage), gds_lab_all))
if (!length(lab)) stop("No overlap across inputs.")
sel_int <- as.integer(label2int[lab])
if (any(is.na(sel_int))) stop("Label->scan_id mapping failed.")

fwrite(data.table(IID=as.character(sel_int)),
       file.path(opt$outdir,"genesis_final_scan_ids.txt"), col.names=FALSE)
fwrite(data.table(label=lab),
       file.path(opt$outdir,"genesis_final_label_ids.txt"), col.names=FALSE)
cat("[cohort_ids] N =", length(sel_int), "\n")
