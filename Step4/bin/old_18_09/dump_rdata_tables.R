#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE); if (!length(args)) stop("usage: dump_rdata_tables.R heritability.RData")
env <- new.env(); load(args[1], envir=env)
rec <- function(obj, name="") {
  if (is.data.frame(obj) || data.table::is.data.table(obj)) {
    cat("TABLE:", name, "  nrow=", nrow(obj), "  cols=", paste(names(obj), collapse=","), "\n", sep="")
    print(utils::head(obj, 3))
  } else if (is.list(obj)) {
    nms <- names(obj); nms[is.null(nms)] <- as.character(seq_along(obj))
    for (i in seq_along(obj)) rec(obj[[i]], paste0(name, if(nchar(name))"$","", nms[i]))
  }
}
for (nm in ls(env)) rec(get(nm, env), nm)

