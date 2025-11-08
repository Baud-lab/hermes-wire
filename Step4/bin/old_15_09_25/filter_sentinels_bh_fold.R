#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(optparse); library(R.utils) })
opt_list <- list(
make_option("--sentinels", type="character", help="sentinels_genus__{genus}__foldX.tsv"),
make_option("--assoc_dir", type="character", help="fold_X/assoc_genus_full__{genus}/by_trait or parent"),
make_option("--bh_q", type="double", default=0.10),
make_option("--out", type="character", default="sentinels_bh.tsv")
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(opt$sentinels), dir.exists(opt$assoc_dir))
S <- fread(opt$sentinels)
setDT(S); setnames(S, tolower(names(S)))
if (!"sentinel_snp"%in%names(S)) {
if ("snp"%in%names(S)) setnames(S, "snp","sentinel_snp") else stop("sentinel file needs column sentinel_snp or snp")
}
# discover assoc file
cand <- list.files(opt$assoc_dir, pattern="assoc\\.tsv\\.gz$", full.names=TRUE, recursive=TRUE)
if (length(cand)==0L) stop("No assoc .tsv.gz under ", opt$assoc_dir)
A <- rbindlist(lapply(cand, function(f) {
tryCatch({
x <- fread(f); setDT(x); setnames(x, tolower(names(x)))
if (!"snp"%in%names(x) && "variant_id"%in%names(x)) setnames(x, "variant_id","snp")
if (!"score.pval"%in%names(x)) {
# allow p columns with different names
pcol <- grep("pval$|p_value$|p$", names(x), value=TRUE)
if (length(pcol)==0L) stop("No p-value column in ", f)
setnames(x, pcol[1], "score.pval")
}
x[, .(snp, score.pval)]
}, error=function(e) data.table())
}))
if (nrow(A)==0L) stop("No assoc rows parsed")
M <- merge(unique(S[, .(sentinel_snp)]), unique(A), by.x="sentinel_snp", by.y="snp", all.x=TRUE)
M[, q := p.adjust(score.pval, method="BH")]
KEEP <- M[!is.na(q) & q < opt$bh_q]
if (nrow(KEEP)==0L) stop("No sentinels pass BH q<", opt$bh_q)
fwrite(KEEP, file=opt$out, sep="\t")
message("[bh] wrote ", opt$out, " rows=", nrow(KEEP))
