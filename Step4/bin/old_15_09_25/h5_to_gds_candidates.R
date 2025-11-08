#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(data.table)
  library(rhdf5)
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
})

# ================================================================
# CLI
# ================================================================
option_list <- list(
  make_option("--h5",          type="character", help="Input HDF5 path"),
  make_option("--genes",       type="character", help="XLSX with Chromosome, Start, End, Gene"),
  make_option("--vers",        type="character", help="'candidate' or 'full'"),
  make_option("--out",         type="character", default="genotypes_candidates.gds",
              help="Output SNPRelate GDS [default: %default]"),
  make_option("--geno-mode",   type="character", default="round",
              help="'round' (default) or 'missing'"),
  make_option("--window_up",   type="integer",   default=0,
              help="Upstream window in bp"),
  make_option("--window_down", type="integer",   default=0,
              help="Downstream window in bp")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- sanity ----
if (!file.exists(opt$h5)) stop("[FATAL] HDF5 not found: ", opt$h5)
if (!file.exists(opt$genes)) stop("[FATAL] genes XLSX not found: ", opt$genes)
if (!(opt$`geno-mode` %in% c("round","missing")))
  stop("[FATAL] --geno-mode must be 'round' or 'missing' (got: ", opt$`geno-mode`, ")")
if (opt$window_up < 0 || opt$window_down < 0)
  stop("[FATAL] windows must be >=0")

cat("[INFO] =============================================\n")
cat("[INFO] h5        : ", opt$h5, "\n", sep="")
cat("[INFO] genes     : ", opt$genes, "\n", sep="")
cat("[INFO] vers      : ", opt$vers, "\n", sep="")
cat("[INFO] out       : ", opt$out, "\n", sep="")
cat("[INFO] geno-mode : ", opt$`geno-mode`, "\n", sep="")
cat("[INFO] win_up/down (bp): ", opt$window_up, " / ", opt$window_down, "\n", sep="")
cat("[INFO] =============================================\n")

# ================================================================
# Read genes and normalize columns
# Expect: "Bioprocess","Gene","Chromosome","Start","End","Number of SNPs","Reference"
# ================================================================
cat("[INFO] Reading genes XLSX ...\n")
putative <- as.data.table(readxl::read_excel(opt$genes))
req <- c("Gene","Chromosome","Start","End")
if (!all(req %in% names(putative))) {
  stop("[FATAL] XLSX must contain columns: ", paste(req, collapse=", "))
}

# Standardize names and types
putative <- putative %>%
  rename(GENE = Gene, CHR = Chromosome, START = Start, END = End) %>%
  mutate(CHR = as.character(CHR),
         START = as.integer(START),
         END   = as.integer(END)) %>%
  select(GENE, CHR, START, END) %>%
  as.data.table()

# Rename for downstream consistency
setnames(putative, c("GENE","Chromosome","Start0","End0"))

# Filter malformed / empty
putative <- putative[!is.na(Start0) & !is.na(End0) & Start0 <= End0]
if (nrow(putative) == 0) stop("[FATAL] No valid gene intervals after filtering.")

# Expanded windows; keep original Start0/End0 for 'in-gene'
putative[, `:=`(
  Start  = pmax(1L, Start0 - as.integer(opt$window_up)),
  End    = End0 + as.integer(opt$window_down),
  GeneID = paste0(Chromosome, ":", Start0, "-", End0)
)]
cat("[INFO] Genes loaded: ", nrow(putative), " across ",
    length(unique(putative$Chromosome)), " chromosomes\n", sep="")

# ================================================================
# Helpers
# ================================================================
mk_intervals <- function(dt, win_up, win_down) {
  body <- dt[, .(Chromosome, start=Start0, end=End0, GeneID)]
  up   <- dt[, .(Chromosome, start=pmax(1L, Start0 - win_up), end=Start0-1L, GeneID)]
  down <- dt[, .(Chromosome, start=End0+1L, end=End0 + win_down, GeneID)]
  body <- body[end >= start]; up <- up[end >= start]; down <- down[end >= start]
  list(body=body, upstream=up, downstream=down)
}
pad2 <- function(x) sprintf("%02d", as.integer(x))
chr_node <- function(h5f, ch) {
  n1 <- paste0("/direct_unpruned_Rn7/chr", as.character(ch))
  n2 <- paste0("/direct_unpruned_Rn7/chr", pad2(ch))
  if (H5Lexists(h5f, n1)) return(n1)
  if (H5Lexists(h5f, n2)) return(n2)
  stop("Chromosome node not found in HDF5 for chr ", ch)
}

# ================================================================
# Open H5, read sample IDs
# ================================================================
cat("[INFO] Opening HDF5 ...\n")
h5f <- H5Fopen(opt$h5)

first_ch <- putative$Chromosome[1]
cat("[INFO] Reading sample IDs from first chromosome: ", first_ch, "\n", sep="")
rfids <- h5read(h5f, paste0(chr_node(h5f, first_ch), "/row_header/sample_ID"))
rfids <- as.character(rfids)
n_scan <- length(rfids)
scan_id <- seq_len(n_scan)
lookup_scan <- data.table(scan_id = scan_id, sample_label = rfids)
cat("[INFO] #samples: ", n_scan, "\n", sep="")

# ================================================================
# Main loop: per-chromosome mapping + GDS parts
# ================================================================
intervals <- mk_intervals(putative, opt$window_up, opt$window_down)
sum_rows <- list()
per_gene <- list()
gds_parts <- character(0)

for (ch in unique(putative$Chromosome)) {
  cat("[INFO] ---- Chromosome ", ch, " ----\n", sep="")
  node <- chr_node(h5f, ch)
  
  # positions
  pos <- as.integer(h5read(h5f, paste0(node, "/col_header/pos")))
  if (!length(pos)) { cat("[WARN] No SNP positions on chr ", ch, "\n", sep=""); next }
  
  # point intervals for foverlaps (start=end=pos)
  dt_pos <- data.table(Chromosome=ch, start=pos, end=pos, idx=seq_along(pos))
  setkey(dt_pos, Chromosome, start, end)
  
  # expanded window table & keys
  sub_exp <- putative[Chromosome == ch, .(Chromosome, start=Start, end=End)]
  if (!nrow(sub_exp)) { cat("[WARN] No gene windows on chr ", ch, "\n", sep=""); next }
  setkey(sub_exp, Chromosome, start, end)
  
  # SNPs in ANY expanded window
  hits_exp <- foverlaps(dt_pos, sub_exp, nomatch=0L)$idx
  hits_exp <- unique(hits_exp)
  if (!length(hits_exp)) { cat("[WARN] No SNPs in expanded windows on chr ", ch, "\n", sep=""); next }
  
  # ---- classify with upstream priority (body -> upstream -> downstream) ----
  body_tbl <- intervals$body[Chromosome == ch];     setkey(body_tbl, Chromosome, start, end)
  up_tbl   <- intervals$upstream[Chromosome == ch]; setkey(up_tbl,   Chromosome, start, end)
  dn_tbl   <- intervals$downstream[Chromosome == ch]; setkey(dn_tbl, Chromosome, start, end)
  
  classify_one <- function(p) {
    hit <- body_tbl[start <= p & end >= p]
    if (nrow(hit)) return(list(gene=hit$GeneID[1], type="in_gene"))
    hit <- up_tbl[start <= p & end >= p]
    if (nrow(hit)) return(list(gene=hit$GeneID[1], type="upstream"))
    hit <- dn_tbl[start <= p & end >= p]
    if (nrow(hit)) return(list(gene=hit$GeneID[1], type="downstream"))
    return(NULL)
  }
  
  cls <- lapply(hits_exp, function(i) classify_one(pos[i]))
  keep <- !vapply(cls, is.null, logical(1))
  if (!any(keep)) { cat("[WARN] Mapped none on chr ", ch, "\n", sep=""); next }
  cls_dt <- rbindlist(cls[keep], fill=TRUE)  # <-- ensures atomic columns
  cdt <- cbind(data.table(idx = hits_exp[keep]), cls_dt)
  setDT(cdt)
  cdt[, gene := as.character(gene)]
  cdt[, type := as.character(type)]
  
  # ---- summaries ----
  n_body <- sum(cdt$type=="in_gene"); n_up <- sum(cdt$type=="upstream"); n_dn <- sum(cdt$type=="downstream")
  cat(sprintf("[INFO] Counts chr %s: in_gene=%d upstream=%d downstream=%d total=%d\n",
              ch, n_body, n_up, n_dn, nrow(cdt)))
  
  sum_rows[[length(sum_rows)+1]] <- data.table(
    Chromosome=ch, n_in_gene=n_body, n_upstream=n_up, n_downstream=n_dn, n_total=nrow(cdt)
  )
  
  # per-gene
  pg <- cdt[, .(n_in_gene=sum(type=="in_gene"),
                n_upstream=sum(type=="upstream"),
                n_downstream=sum(type=="downstream"),
                n_total=.N),
            by=.(GeneID = gene)]
  pg <- merge(pg, putative[, .(GeneID, Chromosome, Start0, End0, Start, End)], by="GeneID", all.x=TRUE)
  per_gene[[length(per_gene)+1]] <- pg
  
  # ---- extract genotypes for selected SNPs and write a per-chr GDS part ----
  cat("[INFO] Reading genotype submatrix for chr ", ch, " ...\n", sep="")
  G <- h5read(h5f, paste0(node, "/matrix"), index = list(1:n_scan, cdt$idx))
  
  if (opt$`geno-mode` == "round") {
    G <- round(G); G[is.na(G)] <- 3; G[G < 0] <- 0; G[G > 2] <- 2
  } else {
    G[!(G %in% c(0,1,2))] <- 3
  }
  storage.mode(G) <- "integer"
  GENO <- t(G)
  
  allele <- tryCatch({
    ref <- h5read(h5f, paste0(node, "/col_header/ref"))[cdt$idx]
    alt <- h5read(h5f, paste0(node, "/col_header/alt"))[cdt$idx]
    paste(ref, alt, sep="/")
  }, error=function(e) rep("A/B", length(cdt$idx)))
  
  snp_pos <- pos[cdt$idx]
  snp_id  <- seq_len(nrow(GENO))
  part    <- sprintf("cand_chr%s.gds", pad2(ch))
  if (file.exists(part)) file.remove(part)
  
  cat("[INFO] Writing GDS part: ", part, " (", length(snp_id), " SNPs)\n", sep="")
  snpgdsCreateGeno(
    gds.fn         = part,
    genmat         = GENO,
    sample.id      = scan_id,
    snp.id         = snp_id,
    snp.chromosome = rep(as.integer(ch), length(snp_id)),
    snp.position   = as.integer(snp_pos),
    snp.allele     = allele,
    snpfirstdim    = TRUE
  )
  g <- openfn.gds(part, readonly=FALSE)
  add.gdsn(g, "sample.label", rfids)
  add.gdsn(g, "snp.label", paste0("chr", ch, ":", snp_pos))
  closefn.gds(g)
  
  gds_parts <- c(gds_parts, part)
}

# close H5
H5Fclose(h5f)

# ================================================================
# Write summaries (TSVs)
# ================================================================
cat("[INFO] Writing mapping summaries ...\n")
summary_dt <- rbindlist(sum_rows, fill=TRUE)
fwrite(summary_dt, "candidate_mapping_summary.tsv")

totals <- summary_dt[, .(
  n_in_gene   = sum(n_in_gene  , na.rm=TRUE),
  n_upstream  = sum(n_upstream , na.rm=TRUE),
  n_downstream= sum(n_downstream, na.rm=TRUE),
  n_total     = sum(n_total    , na.rm=TRUE)
)]
fwrite(totals, "candidate_mapping_total.tsv")

pg_all <- rbindlist(per_gene, fill=TRUE)
fwrite(pg_all, "candidate_mapping_per_gene.tsv")
cat("[INFO] Wrote TSVs: candidate_mapping_summary.tsv, candidate_mapping_total.tsv, candidate_mapping_per_gene.tsv\n")

# ================================================================
# Combine per-chr GDS parts -> final GDS, then SeqArray
# ================================================================
if (!length(gds_parts)) stop("[FATAL] No GDS parts were created; nothing to combine.")

if (file.exists(opt$out)) file.remove(opt$out)
cat("[INFO] Combining ", length(gds_parts), " parts into ", opt$out, " ...\n", sep="")
snpgdsCombineGeno(
  gds.fn              = gds_parts,
  out.fn              = opt$out,
  method              = "position",
  compress.annotation = "ZIP_RA.max",
  compress.geno       = "",
  same.strand         = FALSE,
  snpfirstdim         = TRUE,
  verbose             = TRUE
)
g <- openfn.gds(opt$out, readonly=FALSE)
add.gdsn(g, "sample.label", rfids)
closefn.gds(g)
cat("[INFO] Final candidate GDS: ", opt$out, "\n", sep="")

# scan lookup
lookup_path <- sub("\\.gds$", "_scan_lookup.tsv", opt$out)
fwrite(lookup_scan, lookup_path)
cat("[INFO] Wrote scan lookup: ", lookup_path, "\n", sep="")

# SeqArray conversion
out_seq <- sub("\\.gds$", ".seq.gds", opt$out)
if (file.exists(out_seq)) file.remove(out_seq)
cat("[INFO] Converting to SeqArray: ", out_seq, " ...\n", sep="")
seqSNP2GDS(opt$out, out_seq, verbose=TRUE)
cat("[INFO] Done. SeqArray GDS: ", out_seq, "\n", sep="")

cat("[INFO] All tasks completed successfully.\n")
