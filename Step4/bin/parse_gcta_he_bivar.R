#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(optparse); library(data.table); library(stringr) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--glob", type="character", help="glob like 'he_bivar/*.HEreg'"),
  make_option("--out",  type="character", default="he_bivar_summary.tsv")
)))
files <- Sys.glob(opt$glob)
if (!length(files)) { file.create(opt$out); quit(save="no") }

parse_he <- function(f) {
  fn <- basename(f)
  m  <- str_match(fn, "__(.+)__VS__(.+)\\.HEreg$")
  trait1 <- if (is.na(m[1,2])) NA_character_ else gsub("_"," ",m[1,2])
  trait2 <- if (is.na(m[1,3])) NA_character_ else gsub("_"," ",m[1,3])
  x <- readLines(f, warn=FALSE)
  # GCTA HEreg reports cov_g (genetic covariance) and rg (at the end). We try to parse 'rg'.
  rg_line <- grep("^rg\\s", x, value=TRUE)
  est <- se <- p <- NA_real_
  if (length(rg_line)) {
    parts <- strsplit(rg_line, "\\s+")[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 4) {
      est <- as.numeric(parts[2]); se <- as.numeric(parts[3]); p <- as.numeric(parts[4])
    }
  }
  data.table(file=f, trait1=trait1, trait2=trait2, rg=est, se=se, p=p)
}
DT <- rbindlist(lapply(files, parse_he), fill=TRUE)
fwrite(DT, file=opt$out, sep="\t")
