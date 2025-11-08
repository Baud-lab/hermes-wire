#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(data.table)
})

# -------- CLI --------
option_list <- list(
  make_option("--bim",   type="character", help="Input PLINK .bim file"),
  make_option("--genes", type="character", help="Putative genes XLSX with Chromosome/Start/End"),
  make_option("--magma_up",   type="integer", default=50000, help="MAGMA upstream window [default: %default]"),
  make_option("--magma_down", type="integer", default=50000, help="MAGMA downstream window [default: %default]"),
  make_option("--outdir", type="character", default=".", help="Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

bim_path   <- opt$bim
genes_xlsx <- opt$genes
up         <- as.integer(opt$magma_up)
down       <- as.integer(opt$magma_down)
outdir     <- normalizePath(opt$outdir, mustWork=FALSE)

cat("[make_gene_annotation] ===== START =====\n")
cat("[make_gene_annotation] BIM:", bim_path, "\n")
cat("[make_gene_annotation] Genes XLSX:", genes_xlsx, "\n")
cat("[make_gene_annotation] MAGMA window (bp):", up, "up /", down, "down\n")

# -------- Load BIM --------
bim <- fread(bim_path, header=FALSE, col.names=c("chr","snp","cm","pos","a1","a2"))
# normalize BIM chr to character (don't lose leading X etc.)
bim[, chr := as.character(chr)]
bim[, chr := sub("^chr", "", chr, ignore.case=TRUE)]

# -------- Load genes --------
genes <- as.data.frame(readxl::read_excel(genes_xlsx))
setDT(genes)
req <- c("Chromosome","Start","End")
if (!all(req %in% names(genes))) stop("Genes file must have columns: ", paste(req, collapse=", "))

genes[, `:=`(
  Chromosome = as.character(Chromosome),
  Start      = as.integer(Start),
  End        = as.integer(End)
)]
# swap if Start > End
bad <- which(genes$Start > genes$End)
if (length(bad)) {
  tmp <- genes$Start[bad]; genes$Start[bad] <- genes$End[bad]; genes$End[bad] <- tmp
  message("[make_gene_annotation] Swapped Start/End for ", length(bad), " gene(s) with Start>End")
}

# strip any 'chr' prefix
genes[, Chromosome := sub("^chr", "", Chromosome, ignore.case=TRUE)]

# -------- Harmonize chromosome coding (X/Y/MT) to match BIM ----------
bim_chr_is_numeric <- suppressWarnings(all(!is.na(as.integer(bim$chr))))
if (bim_chr_is_numeric) {
  # PLINK commonly uses 23=X, 24=Y, 26=MT
  genes[Chromosome %in% c("X","x"),  Chromosome := "23"]
  genes[Chromosome %in% c("Y","y"),  Chromosome := "24"]
  genes[Chromosome %in% c("MT","M","Mt","m","MtDNA","mtDNA"), Chromosome := "26"]
} else {
  # ensure both are letters if BIM is not numeric
  genes[Chromosome == "23", Chromosome := "X"]
  genes[Chromosome == "24", Chromosome := "Y"]
  genes[Chromosome == "26", Chromosome := "MT"]
}

# -------- Write MAGMA gene.loc (no window) --------
genes[, GENE_ID := sprintf("chr%s:%s-%s", Chromosome, Start, End)]
gene_loc <- file.path(outdir, "magma.genes.loc")
fwrite(genes[,.(GENE_ID, Chromosome, Start, End)], gene_loc, sep="\t", col.names=FALSE)
cat("[make_gene_annotation] Wrote MAGMA gene.loc:", gene_loc, " nGenes=", nrow(genes), "\n")

# -------- SNP partition for GREML with window (crypt/rest) ----------
# (use data.table non-equi joins for speed)
bim4join <- bim[, .(chr, start=pos, end=pos, snp)]
# widen gene intervals for GREML partition (cap at 1 bp)
genes_w <- copy(genes)[, .(Chromosome,
                           wstart = pmax(1L, Start - up),
                           wend   = End + down,
                           GENE_ID)]
setkey(bim4join, chr, start, end)
setkey(genes_w, Chromosome, wstart, wend)

ov <- foverlaps(bim4join, genes_w,
                by.x=c("chr","start","end"),
                by.y=c("Chromosome","wstart","wend"),
                nomatch=0L)

crypt_snps <- unique(ov$snp)
rest_snps  <- setdiff(bim$snp, crypt_snps)

crypt_file <- file.path(outdir, "crypt_snps.txt")
rest_file  <- file.path(outdir, "rest_snps.txt")
fwrite(data.table(snp=crypt_snps), crypt_file, col.names=FALSE, quote=FALSE)
fwrite(data.table(snp=rest_snps),  rest_file,  col.names=FALSE, quote=FALSE)

# -------- Light summary & per-gene counts ----------
per_gene_counts <- ov[, .N, by=GENE_ID]
summary_file <- file.path(outdir, "gene_annotation_summary.tsv")
summary <- data.table(
  nGenes = nrow(genes),
  nSNP   = nrow(bim),
  nCrypt = length(crypt_snps),
  nRest  = length(rest_snps),
  magma_up_bp = up,
  magma_down_bp = down,
  nGenes_with_SNPs = sum(genes$GENE_ID %in% per_gene_counts$GENE_ID)
)
fwrite(summary, summary_file, sep="\t")

counts_file <- file.path(outdir, "per_gene_bim_overlap.tsv")
fwrite(per_gene_counts[order(-N)], counts_file, sep="\t")

cat("[make_gene_annotation] crypt_snps:", length(crypt_snps),
    " rest_snps:", length(rest_snps), "\n")
cat("[make_gene_annotation] Per-gene overlap: ", counts_file, "\n")
cat("[make_gene_annotation] Summary: ", summary_file, "\n")
cat("[make_gene_annotation] ===== DONE =====\n")
