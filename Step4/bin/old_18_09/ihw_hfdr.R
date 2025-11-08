#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(optparse);library(data.table);library(IHW)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--acat_table"),  # gene, p_acat_h2
  make_option("--cov_table"),   # gene, NSNPS, meanMAF
  make_option("--alpha", type="double", default=0.10),
  make_option("--out")
)))
P <- fread(opt$acat_table); C <- fread(opt$cov_table)
D <- merge(P,C,by="gene",all.x=TRUE)
cov <- pmin(pmax(D$NSNPS,1), 1000)  # example covariate; replace/extend with meanMAF etc.
ihw_res <- ihw(D$p_acat_h2, covariate=cov, alpha=opt$alpha)
D$q_ihw <- adj_pvalues(ihw_res)
fwrite(D, opt$out, sep="\t")

