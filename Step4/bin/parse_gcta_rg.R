#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(optparse); library(data.table); library(stringr) })
opt <- parse_args(OptionParser(option_list=list(
  make_option("--glob", type="character", help="glob like 'bivar_all/*.hsq'"),
  make_option("--out",  type="character", default="rg_summary.tsv")
)))
files <- Sys.glob(opt$glob)
if (!length(files)) { file.create(opt$out); quit(save="no") }

parse_hsq <- function(f) {
  x <- readLines(f, warn=FALSE)
  # trait names are embedded in filename: rg_bivar_all__trait1__VS__trait2.hsq
  fn <- basename(f)
  m  <- str_match(fn, "__(.+)__VS__(.+)\\.hsq$")
  trait1 <- if (is.na(m[1,2])) NA_character_ else gsub("_"," ",m[1,2])
  trait2 <- if (is.na(m[1,3])) NA_character_ else gsub("_"," ",m[1,3])
  # pull rg line: "rg  <est> <se> <p>"
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
DT <- rbindlist(lapply(files, parse_hsq), fill=TRUE)
fwrite(DT, file=opt$out, sep="\t")
