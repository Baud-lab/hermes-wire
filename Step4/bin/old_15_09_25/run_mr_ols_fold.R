#!/usr/bin/env Rscript
# run_mr_ols_fold.R
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--df",       type="character", help="DF_fullmed__GENUS__foldXX.tsv[.gz] (id + genes + last5)"),
  make_option("--raw",      type="character", help="geno__GENUS__foldXX.raw (PLINK --recode A)"),
  make_option("--out_fold", type="character", help="Output: mr_fold__GENUS__foldXX.tsv"),
  make_option("--out_snp",  type="character", help="Output: mr_per_snp__GENUS__foldXX.tsv.gz"),
  make_option("--verbose",  action="store_true", default=TRUE)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$df), !is.null(opt$raw), !is.null(opt$out_fold), !is.null(opt$out_snp))

msg <- function(...) if (isTRUE(opt$verbose)) cat(sprintf(...), "\n")

read_any <- function(path) {
  if (grepl("\\.gz$", path)) fread(cmd=paste("gzip -dc", shQuote(path))) else fread(path)
}

# Load DF (id + genes + last 5 outcomes)
DF <- tryCatch(read_any(opt$df), error=function(e) data.table())
if (nrow(DF) == 0L) {
  msg("[mr] DF is empty; writing empty outputs")
  fwrite(data.table(), file=opt$out_fold, sep="\t")
  fwrite(data.table(), file=opt$out_snp,  sep="\t")
  quit(save="no")
}
stopifnot("id" %in% names(DF))
k <- ncol(DF); stopifnot(k >= 7)
last5 <- (k-4):k
setnames(DF, last5, c("Mediator","OtherGenus","BetaPCo1","Guild3","Glucose"))

# Load PLINK .raw (FID IID PAT MAT SEX PHENOTYPE + SNPs)
RAW <- fread(opt$raw, sep=" ")
if (nrow(RAW) == 0L || ncol(RAW) < 7L) {
  msg("[mr] .raw is empty or malformed; writing empty outputs")
  fwrite(data.table(), file=opt$out_fold, sep="\t")
  fwrite(data.table(), file=opt$out_snp,  sep="\t")
  quit(save="no")
}
setnames(RAW, 1:2, c("FID","IID"))

# Merge DF$id ↔ RAW$IID
setkey(DF, id); setkey(RAW, IID)
M <- RAW[DF, nomatch=0]
if (nrow(M) == 0L) stop("[mr] No IID match between DF$id and .raw IID")

# SNP columns from RAW: everything after the first 6 metadata cols
raw_snp_cols <- setdiff(names(RAW), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
if (length(raw_snp_cols) == 0L) {
  msg("[mr] No SNP columns in .raw; writing empty outputs")
  fwrite(data.table(), file=opt$out_fold, sep="\t")
  fwrite(data.table(), file=opt$out_snp,  sep="\t")
  quit(save="no")
}

# Helpers
ols_slope <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 8) return(c(beta=NA_real_, se=NA_real_, F=NA_real_))
  xm <- mean(x); ym <- mean(y)
  Sxx <- sum((x - xm)^2)
  if (Sxx <= 0) return(c(beta=NA_real_, se=NA_real_, F=NA_real_))
  beta <- sum((x - xm)*(y - ym)) / Sxx
  resid <- y - (beta*(x - xm) + ym)
  s2 <- sum(resid^2) / (n - 2)
  se <- sqrt(s2 / Sxx)
  F  <- (beta*beta) / (se*se)
  c(beta=beta, se=se, F=F)
}

traits <- c("OtherGenus","BetaPCo1","Guild3","Glucose")

# Loop SNPs (column-wise to avoid huge copies)
res_list <- vector("list", length(raw_snp_cols))
for (i in seq_along(raw_snp_cols)) {
  sn <- raw_snp_cols[i]
  g  <- M[[sn]]
  bx <- ols_slope(g, M$Mediator)
  if (is.na(bx["beta"])) {
    res_list[[i]] <- data.table(SNP=sn, bx=NA_real_, se_bx=NA_real_, F_first=NA_real_,
                                trait=NA_character_, by=NA_real_, se_by=NA_real_,
                                ratio=NA_real_, se_ratio=NA_real_)
    next
  }
  tmp <- rbindlist(lapply(traits, function(tr) {
    by <- ols_slope(g, M[[tr]])
    if (is.na(by["beta"]) || bx["beta"] == 0) {
      return(data.table(SNP=sn, bx=bx["beta"], se_bx=bx["se"], F_first=bx["F"],
                        trait=tr, by=NA_real_, se_by=NA_real_,
                        ratio=NA_real_, se_ratio=NA_real_))
    }
    ratio <- by["beta"] / bx["beta"]
    se_r  <- sqrt((by["se"]^2)/(bx["beta"]^2) + ((by["beta"]^2)*(bx["se"]^2))/(bx["beta"]^4))
    data.table(SNP=sn, bx=bx["beta"], se_bx=bx["se"], F_first=bx["F"],
               trait=tr, by=by["beta"], se_by=by["se"], ratio=ratio, se_ratio=se_r)
  }), fill=TRUE)
  res_list[[i]] <- tmp
}
RES <- rbindlist(res_list, fill=TRUE)

# Write per-SNP results (gz)
fwrite(RES, file=opt$out_snp, sep="\t")

# IVW per trait
IVW <- RES[is.finite(ratio) & is.finite(se_ratio),
           .( n_snp = .N,
              ivw_beta = { w <- 1/(se_ratio^2); sum(w*ratio)/sum(w) },
              ivw_se   = sqrt(1/sum(1/(se_ratio^2))),
              F_first_min = suppressWarnings(min(F_first, na.rm=TRUE)),
              F_first_med = suppressWarnings(median(F_first, na.rm=TRUE)),
              F_first_mean= suppressWarnings(mean(F_first, na.rm=TRUE))
            ),
           by=trait]
if (nrow(IVW) > 0L) {
  IVW[, ivw_z := ivw_beta/ivw_se ]
  IVW[, ivw_p := 2*pnorm(-abs(ivw_z)) ]
  IVW[, ivw_q := p.adjust(ivw_p, method="BH") ]
}
fwrite(IVW, file=opt$out_fold, sep="\t")
msg("[mr] Wrote: %s and %s", opt$out_fold, opt$out_snp)
