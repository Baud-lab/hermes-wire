#!/usr/bin/env Rscript
# Build GRangesList of genes that intersect the cohort PLINK panel
# Input:
#   --gene_loc magma.genes.loc      (tab: GENE_ID, CHR, START, END; >=4 columns)
#   --bim      candidates_kept.bim  (cohort-subset)
#   --up/--down window padding (bp) already applied upstream usually; we re-apply if desired
# Output:
#   genes.grl.rds  (GRangesList named by GENE_ID)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

opt_list <- list(
  make_option("--gene_loc", type="character", help="magma.genes.loc"),
  make_option("--bim",      type="character", help="PLINK .bim for cohort subset"),
  make_option("--up",       type="integer", default=0),
  make_option("--down",     type="integer", default=0),
  make_option("--out",      type="character", default="genes.grl.rds")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(file.exists(opt$gene_loc), file.exists(opt$bim))

cat("[make_gene_granges] Reading gene_loc…\n")
gl <- fread(opt$gene_loc, header = FALSE, data.table = FALSE)
if (ncol(gl) < 4) stop("gene_loc must have >=4 columns: GENE_ID CHR START END")
colnames(gl)[1:4] <- c("GENE_ID","CHR","START","END")
gl <- gl[!is.na(gl$GENE_ID) & !is.na(gl$CHR) & !is.na(gl$START) & !is.na(gl$END), ]

# normalize chr names (keep as character)
gl$CHR <- as.character(gl$CHR)

# pad windows if requested
if (opt$up > 0L || opt$down > 0L) {
  gl$START <- pmax(1L, as.integer(gl$START) - as.integer(opt$up))
  gl$END   <- as.integer(gl$END) + as.integer(opt$down)
}

cat("[make_gene_granges] Reading BIM…\n")
bim <- fread(opt$bim, data.table = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))
bim$CHR <- as.character(bim$CHR)
bim$BP  <- as.integer(bim$BP)

# keep only chromosomes present in BIM (prevents empty ranges on non-existent contigs)
keep_chr <- intersect(unique(gl$CHR), unique(bim$CHR))
gl <- gl[gl$CHR %in% keep_chr, , drop = FALSE]

# create GRanges per gene
cat("[make_gene_granges] Building GRangesList…\n")
gr <- GRanges(seqnames = gl$CHR,
              ranges   = IRanges(start = as.integer(gl$START), end = as.integer(gl$END)),
              GENE_ID  = gl$GENE_ID)
# Split by GENE_ID
split_index <- split(seq_along(gr), mcols(gr)$GENE_ID)
grl <- GRangesList(lapply(split_index, function(ix) gr[ix]))
names(grl) <- names(split_index)

# optional light pruning: drop genes with zero overlap to BIM positions
cat("[make_gene_granges] Pruning genes with no SNPs in BIM…\n")
bim_gr <- GRanges(seqnames = bim$CHR, ranges = IRanges(bim$BP, bim$BP))
has_overlap <- vapply(grl, function(g) suppressWarnings(length(findOverlaps(g, bim_gr, select="first")) > 0), logical(1))
grl <- grl[has_overlap]

cat("[make_gene_granges] Genes retained:", length(grl), "\n")
saveRDS(grl, file = opt$out)
cat("[make_gene_granges] Wrote:", opt$out, "\n")
