#!/usr/bin/env Rscript
# meta_mr_across_folds.R
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--inputs_list", type="character", help="Text file with mr_fold__*.tsv, one per line"),
  make_option("--out",         type="character", help="Output: mr_meta__GENUS.tsv"),
  make_option("--verbose",     action="store_true", default=TRUE)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$inputs_list), !is.null(opt$out))

msg <- function(...) if (isTRUE(opt$verbose)) cat(sprintf(...), "\n")

files <- readLines(opt$inputs_list)
if (length(files) == 0L) {
  msg("[meta] No inputs; writing empty")
  fwrite(data.table(), file=opt$out, sep="\t")
  quit(save="no")
}

L <- lapply(files, function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
  x <- tryCatch(fread(f), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  req <- c("trait","ivw_beta","ivw_se","ivw_p","ivw_q")
  if (!all(req %in% names(x))) return(NULL)
  x[, ..req]
})
DT <- rbindlist(L, fill=TRUE)
if (nrow(DT) == 0L) {
  msg("[meta] All inputs empty; writing empty")
  fwrite(data.table(), file=opt$out, sep="\t")
  quit(save="no")
}

META <- DT[is.finite(ivw_beta) & is.finite(ivw_se),
           .( n_folds = .N,
              beta = { w <- 1/(ivw_se^2); sum(w*ivw_beta)/sum(w) },
              se   = sqrt(1/sum(1/(ivw_se^2))) ),
           by=trait]
if (nrow(META) > 0L) {
  META[, z := beta/se ]
  META[, p := 2*pnorm(-abs(z)) ]
  META[, q := p.adjust(p, method="BH") ]
}
setcolorder(META, c("trait","n_folds","beta","se","z","p","q"))
fwrite(META, file=opt$out, sep="\t")
msg("[meta] Wrote: %s", opt$out)
