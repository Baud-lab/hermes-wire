#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--selected_genes_rdata", type="character"),
  make_option("--bim", type="character"),
  make_option("--out_prefix", type="character", default="snp_extract")
)))
load(opt$selected_genes_rdata)   # expects 'correlations' and/or 'genotypes'
if (!exists("correlations")) stop("selected_genes_rdata must include 'correlations'")
# correlations must have chr,pos or 'gene' (chr:pos). Build set of sentinel or per-gene SNPs
DT <- as.data.table(correlations)
if (!"gene" %in% names(DT) && all(c("chr","pos") %in% names(DT)))
  DT[, gene := paste(chr, pos, sep=":")]
DT <- unique(DT[!is.na(gene), .(gene)])
BIM <- fread(opt$bim, header=FALSE)
setnames(BIM, c("CHR","SNP","CM","BP","A1","A2"))
BIM[, tag := paste(CHR, BP, sep=":")]
focus <- unique(BIM$SNP[BIM$tag %in% DT$gene])
fwrite(data.table(SNP=focus), file=paste0(opt$out_prefix,"_focus.txt"),
       col.names=FALSE, quote=FALSE)

bg <- setdiff(BIM$SNP, focus)
fwrite(data.table(SNP=bg), file=paste0(opt$out_prefix,"_bg.txt"),
       col.names=FALSE, quote=FALSE)
