#!/usr/bin/env Rscript
# Simple synthesis: correlate chr:pos to GDS, align samples via scan_ids/labels, save RData
# Inputs:
#   --corr        correlations_by_gene__Bacteroides.tsv (must have chr,pos; p/q optional)
#   --genes_xlsx  putative_genes_final.xlsx (Gene, Chromosome, Start, End, Bioprocess)
#   --gds         genotypes_candidates.seq.gds
#   --scan_lookup genotypes_candidates_scan_lookup.tsv (scan_id,sample_label)
#   --scan_ids    genesis_final_scan_ids.txt
#   --label_ids   genesis_final_label_ids.txt
#
# Output:  selected_genes.RData  with:
#   correlations : correlations table annotated (Gene_name, Bioprocess, Category, variant.id)
#   genotypes    : data.frame sample + one column per Gene_name with 0/1/2 ALT dosages
#                  (rows ordered as label_ids; rows with all-NA genotypes removed)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(readxl)
  library(SeqArray)
  library(SeqVarTools)
})

logf <- function(...) cat("[synthesize]", format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--corr",        type="character"),
  make_option("--genes_xlsx",  type="character"),
  make_option("--gds",         type="character"),
  make_option("--scan_lookup", type="character"),
  make_option("--scan_ids",    type="character"),
  make_option("--label_ids",   type="character")
)))

## --------- basic guards ----------
stopifnot(!is.null(opt$corr),        file.exists(opt$corr))
stopifnot(!is.null(opt$genes_xlsx),  file.exists(opt$genes_xlsx))
stopifnot(!is.null(opt$gds),         file.exists(opt$gds))
stopifnot(!is.null(opt$scan_lookup), file.exists(opt$scan_lookup))
stopifnot(!is.null(opt$scan_ids),    file.exists(opt$scan_ids))
stopifnot(!is.null(opt$label_ids),   file.exists(opt$label_ids))

## --------- 1) Read correlations + genes ----------
logf("Reading correlations:", basename(opt$corr))
C <- fread(opt$corr)
if (!nrow(C)) stop("Empty correlations file.")
names(C) <- tolower(names(C))
if (!all(c("chr","pos") %in% names(C))) stop("Correlations must contain 'chr' and 'pos'.")

# normalize chr/pos
C[, chr := as.character(chr)]
C[, chr := tolower(trimws(chr))]
C[, chr := sub("^chr", "", chr)]
C[, pos := as.integer(pos)]

logf("Reading genes excel:", basename(opt$genes_xlsx))
GX <- as.data.table(read_excel(opt$genes_xlsx))
need <- c("Gene","Chromosome","Start","End","Bioprocess")
miss <- setdiff(need, names(GX))
if (length(miss)) stop("genes_xlsx missing columns: ", paste(miss, collapse = ", "))
setnames(GX, c("Gene","Chromosome","Start","End","Bioprocess"),
         c("Gene","Chromosome","Start","End","Bioprocess"))
GX[, Chromosome := as.character(Chromosome)]
GX[, Chromosome := tolower(trimws(Chromosome))]
GX[, Chromosome := sub("^chr", "", Chromosome)]
GX[, Start := as.numeric(Start)]
GX[, End   := as.numeric(End)]

## --------- 2) Annotate Gene_name / Bioprocess / Category ----------
logf("Annotating Gene_name via interval overlap (pos in [Start, End] & chr match)")
C[, Gene_name := NA_character_]
C[, Bioprocess := NA_character_]
C[, Category := NA_character_]

chrs <- sort(unique(C$chr))
for (ch in chrs) {
  idx_c <- which(C$chr == ch & !is.na(C$pos))
  if (!length(idx_c)) next
  gxsub <- GX[Chromosome == ch & !is.na(Start) & !is.na(End)]
  if (!nrow(gxsub)) next
  for (i in idx_c) {
    p <- C$pos[i]
    hit <- gxsub[Start <= p & End >= p]
    if (nrow(hit) > 0) {
      C$Gene_name[i] <- as.character(hit$Gene[1])
      bp <- as.character(hit$Bioprocess[1])
      if (!is.na(bp)) {
        C$Bioprocess[i] <- sub("\\s*\\(.*$","", bp)
        C$Category[i]   <- ifelse(grepl("\\(", bp), sub(".*\\(([^\\)]+)\\).*","\\1", bp), NA_character_)
      }
    }
  }
}
C <- C[!is.na(Gene_name)]
if (!nrow(C)) stop("No (chr,pos) entries fell within any gene interval.")

## --------- 3) Pick one sentinel per Gene_name (by q then p then first) ----------
logf("Selecting one sentinel (chr:pos) per Gene_name (by q, then p, else first).")
if ("q" %in% names(C)) {
  setorder(C, is.na(q), q, is.na(p), p)
} else if ("p" %in% names(C)) {
  setorder(C, is.na(p), p)
}
C1 <- C[!duplicated(Gene_name)]
logf("Genes with a chosen sentinel:", nrow(C1))

## --------- 4) Align samples exactly like run_genesis_genus.R ----------
logf("Opening GDS and aligning cohort samples.")
gds <- seqOpen(opt$gds); on.exit(seqClose(gds), add=TRUE)
lookup <- fread(opt$scan_lookup)
names(lookup) <- tolower(names(lookup))
if (!all(c("scan_id","sample_label") %in% names(lookup)))
  stop("scan_lookup must have 'scan_id' and 'sample_label'.")

lookup[, scan_id := as.integer(scan_id)]
lookup[, sample_label := as.character(sample_label)]

scan_ids  <- scan(opt$scan_ids,  what=integer(),  quiet=TRUE)
label_ids <- scan(opt$label_ids, what=character(), quiet=TRUE)
if (length(scan_ids) != length(label_ids))
  stop("scan_ids and label_ids length differ.")

# authoritative labels via scan_ids -> lookup
map_lab <- as.character(lookup$sample_label[match(scan_ids, lookup$scan_id)])
if (any(is.na(map_lab))) stop("Some scan_ids are missing in scan_lookup.")
if (!all(label_ids == map_lab)) {
  warning("label_ids do not exactly match lookup; using labels derived from scan_ids.")
  label_ids <- map_lab
}

# Ensure all scan_ids present in GDS, then filter to exactly those samples
g_samp <- as.integer(seqGetData(gds, "sample.id"))
missing_in_gds <- setdiff(scan_ids, g_samp)
if (length(missing_in_gds)) {
  stop("Some scan_ids are not present in GDS sample.id, e.g.: ",
       paste(utils::head(missing_in_gds, 10), collapse=", "))
}
seqResetFilter(gds, verbose=FALSE)
seqSetFilter(gds, sample.sel = g_samp %in% scan_ids, verbose=FALSE)
cur_ids <- as.integer(seqGetData(gds, "sample.id"))
# row reordering index to match scan_ids order exactly
row_idx <- match(scan_ids, cur_ids)
if (any(is.na(row_idx))) stop("Internal sample mapping error after setting sample.sel.")

## --------- 5) Map sentinels (chr:pos) to GDS variant.id ----------
logf("Extracting 0/1/2 alt-allele dosages for selected sentinels from GDS.")
g_chr <- seqGetData(gds, "chromosome")
g_pos <- seqGetData(gds, "position")
g_vid <- seqGetData(gds, "variant.id")

# normalize chr/pos to same type/format as correlations
VMAP <- data.table(chr = as.character(g_chr),
                   pos = as.integer(g_pos),
                   variant.id = as.integer(g_vid))
VMAP[, chr := tolower(trimws(chr))]
VMAP[, chr := sub("^chr", "", chr)]

# keep first occurrence per (chr,pos)
setkey(VMAP, chr, pos)
VMAP <- VMAP[!duplicated(paste(chr, pos))]
# match
want <- C1[, .(chr, pos, Gene_name)]
setkey(want, chr, pos)
MM <- VMAP[want, nomatch=0L]  # inner join on chr,pos
# MM now has: chr, pos, variant.id, Gene_name (from 'want' through i.* columns)
# data.table keeps only VMAP cols unless we merge differently; rebuild clean table:
MM <- merge(want, VMAP, by=c("chr","pos"), all.x=FALSE, all.y=FALSE)
MM <- unique(MM, by=c("chr","pos","variant.id","Gene_name"))

if (!nrow(MM)) {
  logf("WARN: No sentinel columns matched by chr:pos in GDS.")
  correlations <- C1
  genotypes <- data.frame(sample = label_ids, stringsAsFactors = FALSE)
  # drop rows with all NA genotype columns (none exist yet, so nothing to drop)
  save(correlations, genotypes, file = "selected_genes.RData")
  logf("Saved: selected_genes.RData")
  logf("DONE")
  quit(save="no", status=0)
}

# keep one variant per gene (if same gene has multiple exact matches, choose first)
MM <- MM[!duplicated(Gene_name)]

# Set variant filter and pull dosages (samples × variants)
seqResetFilter(gds, verbose=FALSE)
seqSetFilter(gds,
             sample.sel  = g_samp %in% scan_ids,
             variant.id  = MM$variant.id,
             verbose=FALSE)

D <- altDosage(gds)   # rows = current sample.id order; cols = in same order as MM$variant.id
# reorder rows to match scan_ids and rename to label_ids
D <- D[row_idx, , drop=FALSE]
rownames(D) <- label_ids

## --------- 6) Assemble 'genotypes' and final 'correlations' ----------
# unique, safe column names from Gene_name
gene_cols <- as.character(MM$Gene_name)
dup <- duplicated(gene_cols)
if (any(dup)) {
  gene_cols[dup] <- paste0(gene_cols[dup], "_dup", seq_len(sum(dup)))
}

genotypes <- data.frame(sample = label_ids, stringsAsFactors = FALSE)
for (j in seq_len(ncol(D))) {
  genotypes[[ gene_cols[j] ]] <- suppressWarnings(as.integer(D[, j]))
}

# drop rows where ALL genotype columns are NA (your 2 no-genotype samples)
if (ncol(genotypes) > 1) {
  has_any <- rowSums(!is.na(genotypes[ , -1, drop=FALSE])) > 0
  dropped <- sum(!has_any)
  if (dropped > 0) logf(sprintf("Dropping %d samples with no genotype across all sentinels.", dropped))
  genotypes <- genotypes[has_any, , drop=FALSE]
}

# correlations table with chosen variant.id
correlations <- merge(C1, MM[,.(chr,pos,Gene_name,variant.id)], by=c("chr","pos","Gene_name"), all.x=FALSE)
# keep tidy order
keep_corr <- intersect(c("genus","gene","model_used","beta_used","p","q",
                         "chr","pos","Gene_name","Bioprocess","Category","variant.id"),
                       names(correlations))
correlations <- correlations[, ..keep_corr]

## --------- 7) Save ----------
save(correlations, genotypes, file = "selected_genes.RData")
logf("Saved: selected_genes.RData")
logf("DONE")
