#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt_list <- list(
make_option("--bh_kept", type="character", help="sentinels_bh__*.tsv"),
make_option("--pgls_keep", type="character", help="pgls_keep__*.tsv"),
make_option("--out", type="character", default="sentinels_for_mediation.tsv")
)
opt <- parse_args(OptionParser(option_list=opt_list))
B <- fread(opt$bh_kept); P <- fread(opt$pgls_keep)
setnames(B, tolower(names(B))); setnames(P, tolower(names(P)))
B[, key := if ("sentinel_snp"%in%names(B)) sentinel_snp else sentinel]
P[, key := if ("variant"%in%names(P)) variant else paste(chr,pos,sep=":")]
K <- merge(B[, .(key)], P[, .(key)], by="key")
if (nrow(K)==0L) stop("Intersection empty")
setnames(K, "key","snp")
fwrite(unique(K), file=opt$out, sep="\t")
message("[gate] wrote ", opt$out, " rows=", nrow(unique(K)))
