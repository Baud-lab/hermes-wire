#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
args <- commandArgs(trailingOnly=TRUE)
kv <- function(k) sub(paste0(".*",k,"="),"",grep(paste0("^",k,"="),args,val=TRUE))
inp <- kv("--mr_input"); pref <- kv("--out_prefix")

H <- fread(inp)
setDT(H)
if (!nrow(H)) quit(status=0)

# split by gene/outcome; run standard MR set
keys <- unique(H[, .(gene, outcome)])
for (i in seq_len(nrow(keys))) {
  g <- keys$gene[i]; o <- keys$outcome[i]
  D <- H[gene==g & outcome==o]
  if (nrow(D) < 1) next
  res_snp <- mr_singlesnp(D)
  res_mr  <- mr(D, method_list=c("mr_ivw_fe","mr_egger_regression","mr_weighted_median"))
  pleio   <- mr_pleiotropy_test(D)
  het     <- mr_heterogeneity(D)

  main <- res_mr[, .(exposure=id.exposure, outcome=id.outcome, method, b, se, pval, nsnp)]
  main[, `:=`(gene=g, outcome=o)]
  qc <- merge(
    pleio[, .(egger_intercept, se, pval)], 
    het[, .(Q, Q_df, Q_pval= Q_pval %||% Q_p)], all=TRUE
  )
  fn_pref <- paste0(pref, "__gene_", g, "__out_", o)
  fwrite(main, paste0("mr_results__", fn_pref, ".tsv"), sep="\t")
  fwrite(qc,   paste0("mr_qc__",      fn_pref, ".tsv"), sep="\t")
}
