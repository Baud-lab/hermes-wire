#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--scan_ids",  type="character"),
  make_option("--label_ids", type="character"),
  make_option("--metadata",  type="character"),
  make_option("--seed",      type="integer", default=17),
  make_option("--outdir",    type="character", default=".")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

scan_ids  <- scan(opt$scan_ids,  what=integer(),  quiet=TRUE)
label_ids <- scan(opt$label_ids, what=character(), quiet=TRUE)
stopifnot(length(scan_ids) == length(label_ids))
dt <- data.table(scan=as.integer(scan_ids), label=as.character(label_ids))

meta <- fread(opt$metadata)
id_col <- intersect(c("sample","Sample","RFID","ID","iid"), names(meta))[1]
if (is.na(id_col)) stop("Metadata needs an ID column among sample/Sample/RFID/ID/iid")
meta <- meta[!duplicated(meta[[id_col]])]
setnames(meta, id_col, "label")

# build grouping unit (dam|cage); fallback to single available or to label
if (!("dam" %in% names(meta)) && !("cage" %in% names(meta))) {
  meta[, grp := label]
} else if ("dam" %in% names(meta) && "cage" %in% names(meta)) {
  meta[, grp := paste0(as.character(dam), "|", as.character(cage))]
} else if ("dam" %in% names(meta)) {
  meta[, grp := as.character(dam)]
} else {
  meta[, grp := as.character(cage)]
}

dt <- merge(dt, meta[, .(label, grp)], by="label", all.x=TRUE)
if (any(is.na(dt$grp))) {
  missing <- unique(dt$label[is.na(dt$grp)])
  stop("Missing grouping info for ", length(missing), " samples; first: ", paste(head(missing,3), collapse=", "))
}

set.seed(opt$seed)
ug <- unique(dt$grp)
perm <- sample(ug, length(ug))
k <- ceiling(length(perm)/2)
grpA <- perm[seq_len(k)]
grpB <- setdiff(ug, grpA)

splitA <- dt[grp %in% grpA]
splitB <- dt[grp %in% grpB]

writeLines(as.character(splitA$scan),  file.path(opt$outdir,"splitA_scan_ids.txt"))
writeLines(as.character(splitA$label), file.path(opt$outdir,"splitA_label_ids.txt"))
writeLines(as.character(splitB$scan),  file.path(opt$outdir,"splitB_scan_ids.txt"))
writeLines(as.character(splitB$label), file.path(opt$outdir,"splitB_label_ids.txt"))

summ <- data.table(
  n_total = nrow(dt),
  nA = nrow(splitA), nB = nrow(splitB),
  n_grp = length(ug), n_grpA = length(grpA), n_grpB = length(grpB)
)
fwrite(summ, file.path(opt$outdir,"two_split_summary.tsv"), sep="\t")
cat("[two_splits] n_total=", nrow(dt), " nA=", nrow(splitA), " nB=", nrow(splitB),
    " groups=", length(ug), "\n", sep="")

