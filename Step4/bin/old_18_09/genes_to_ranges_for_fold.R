#!/usr/bin/env Rscript
# genes_to_ranges_for_fold.R
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--genes",   type="character", help="fullmed_genes__GENUS__foldXX.txt (one gene per line)"),
  make_option("--gene_loc",type="character", help="magma.genes.loc (columns: GENE_ID, CHR, START, END)"),
  make_option("--out",     type="character", help="Output: ranges__GENUS__foldXX.txt (CHR START END GENE_ID)"),
  make_option("--verbose", action="store_true", default=TRUE)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$genes), !is.null(opt$gene_loc), !is.null(opt$out))

msg <- function(...) if (isTRUE(opt$verbose)) cat(sprintf(...), "\n")

G <- fread(opt$genes, header=FALSE)$V1
if (length(G) == 0L) {
  msg("[ranges] Empty gene list; writing empty file: %s", opt$out)
  fwrite(data.table(), file=opt$out, sep="\t", col.names=FALSE)
  quit(save="no")
}

L <- fread(opt$gene_loc, sep="\t", header=TRUE)
setnames(L, tolower(names(L)))  # expect gene_id, chr, start, end
stopifnot(all(c("gene_id","chr","start","end") %in% names(L)))

W <- L[tolower(gene_id) %in% tolower(G), .(chr, start, end, gene_id)]
if (nrow(W) == 0L) {
  msg("[ranges] No coordinates found for the requested genes; writing empty")
  fwrite(data.table(), file=opt$out, sep="\t", col.names=FALSE)
  quit(save="no")
}

# PLINK range format: CHR START END NAME (1-based, inclusive)
fwrite(W, file=opt$out, sep="\t", col.names=FALSE)
msg("[ranges] Wrote: %s (rows=%d)", opt$out, nrow(W))
