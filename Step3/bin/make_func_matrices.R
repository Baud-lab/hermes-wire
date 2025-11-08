#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
argv <- setNames(as.list(rep(NA_character_, length(args))), NULL)
# Tiny key=val parser: --key value
getopt <- function(key, default=NULL) {
  if (!paste0("--",key) %in% args) return(default)
  args[which(args==paste0("--",key))+1]
}

ann_dir   <- getopt("ann_dir")
functions <- getopt("functions", "eggNOG_OGs,CAZy,EC,KEGG_ko,KEGG_Module,KEGG_Pathway,PFAMs,COG_category,COG_pathway")
cog_def   <- getopt("cog_def", "")
out_pref  <- getopt("out_prefix", ".")

if (is.null(ann_dir) || !dir.exists(ann_dir)) stop("Provide --ann_dir with .annotations.tsv files")
fun_list <- strsplit(functions, ",")[[1]] |> trimws()

# Load COG defs if needed
cog_map <- NULL
if ("COG_pathway" %in% fun_list) {
  if (!nzchar(cog_def) || !file.exists(cog_def)) stop("COG_pathway requested but --cog_def not found")
  cog_dt <- fread(cog_def, header = FALSE)
  # Expected columns in your file:
  # V1=COG, V2=Type, V3=Name, V4=Gene, V5=Pathway
  setnames(cog_dt, c("COG","Type","Name","Gene","Pathway","col6","col7")[seq_len(ncol(cog_dt))])
  cog_map <- cog_dt[, .(COG, Pathway)]
}

ann_files <- list.files(ann_dir, pattern = "\\.annotations\\.tsv$", full.names = TRUE)
if (!length(ann_files)) stop("No *.annotations.tsv in ", ann_dir)

message("[make_func_matrices] N annotation files: ", length(ann_files))

# helper: split a comma-separated column to counts
count_split <- function(vec, split_char = ",") {
  vec <- vec[!is.na(vec) & nzchar(vec)]
  if (!length(vec)) return(integer())
  items <- unlist(strsplit(vec, split_char, fixed = TRUE))
  items <- trimws(items)
  table(items[nzchar(items)])
}

# For COG_category, split by single letters (no commas)
count_cog_cat <- function(vec) {
  vec <- vec[!is.na(vec) & nzchar(vec)]
  if (!length(vec)) return(integer())
  items <- unlist(strsplit(vec, ""))
  items <- items[nzchar(items)]
  table(items)
}

# Extract COG/ENOG id before '@'
get_cog_enog <- function(eggNOG_OGs) {
  x <- strsplit(eggNOG_OGs, "@", fixed=TRUE)
  sapply(x, function(z) if (length(z)) z[[1]] else NA_character_)
}

# Build matrices as named list: name -> matrix (rows=function ids; cols=genomes)
mats <- vector("list", length(fun_list)); names(mats) <- fun_list
genome_ids <- character(0)

for (f in ann_files) {
  bn <- basename(f)
  stem <- sub("\\.annotations\\.tsv$", "", bn)
  genome_ids <- c(genome_ids, stem)

  DT <- fread(f, sep="\t", header = TRUE, quote = "", na.strings = c("", "-"))
  cols <- names(DT)

  for (fn in fun_list) {
    cnts <- integer()
    if (fn == "eggNOG_OGs" && "eggNOG_OGs" %in% cols) {
      ids <- get_cog_enog(DT$eggNOG_OGs)
      ids <- ids[!is.na(ids) & nzchar(ids)]
      cnts <- table(ids)
    } else if (fn == "CAZy" && "CAZy" %in% cols) {
      cnts <- count_split(DT$CAZy)
    } else if (fn == "EC" && "EC" %in% cols) {
      cnts <- count_split(DT$EC)
    } else if (fn == "KEGG_ko" && "KEGG_ko" %in% cols) {
      cnts <- count_split(DT$KEGG_ko)
    } else if (fn == "KEGG_Module" && "KEGG_Module" %in% cols) {
      cnts <- count_split(DT$KEGG_Module)
    } else if (fn == "KEGG_Pathway" && "KEGG_Pathway" %in% cols) {
      cnts <- count_split(DT$KEGG_Pathway)
      # keep only map pathways
      if (length(cnts)) cnts <- cnts[grepl("^map", names(cnts))]
    } else if (fn == "PFAMs" && "PFAMs" %in% cols) {
      cnts <- count_split(DT$PFAMs)
    } else if (fn == "COG_category" && "COG_category" %in% cols) {
      cnts <- count_cog_cat(DT$COG_category)
    } else if (fn == "COG_pathway") {
      if (!is.null(cog_map) && "eggNOG_OGs" %in% cols) {
        ids <- get_cog_enog(DT$eggNOG_OGs)
        ids <- ids[!is.na(ids) & nzchar(ids)]
        ids <- ids[grepl("^COG", ids)]
        if (length(ids)) {
          pw <- cog_map$Pathway[match(ids, cog_map$COG)]
          pw <- pw[!is.na(pw) & nzchar(pw)]
          cnts <- table(pw)
        }
      }
    } else {
      next
    }

    # Accumulate into a sparse-like list; we’ll densify at the end
    if (!length(cnts)) next
    if (is.null(mats[[fn]])) {
      mats[[fn]] <- list(ids=names(cnts), vals=as.integer(cnts), cols=c(stem), n=1L)
    } else {
      mats[[fn]]$ids  <- c(mats[[fn]]$ids,  names(cnts))
      mats[[fn]]$vals <- c(mats[[fn]]$vals, as.integer(cnts))
      mats[[fn]]$cols <- c(mats[[fn]]$cols, rep(stem, length(cnts)))
    }
  }
}

genome_ids <- unique(genome_ids)

# Densify into matrices
final <- list()
qc <- list()

for (fn in fun_list) {
  acc <- mats[[fn]]
  if (is.null(acc)) next
  fun_ids <- unique(acc$ids)
  M <- matrix(0L, nrow = length(fun_ids), ncol = length(genome_ids),
              dimnames = list(fun_ids, genome_ids))
  # fill
  idx_fun <- match(acc$ids, rownames(M))
  idx_gen <- match(acc$cols, colnames(M))
  for (i in seq_along(acc$vals)) {
    M[idx_fun[i], idx_gen[i]] <- M[idx_fun[i], idx_gen[i]] + acc$vals[i]
  }
  final[[fn]] <- M
  qc[[fn]] <- data.table(FunctionType = fn,
                         n_functions = nrow(M),
                         n_genomes = ncol(M),
                         nnz = sum(M > 0))
  assign(paste0("func_matrix_", fn), M, inherits = TRUE)
}

qc_dt <- rbindlist(qc)
fwrite(qc_dt, file = file.path(out_pref, "qc_functions_summary.tsv"), sep = "\t")

# Save: combined + eggNOG only (for backward compatibility with your analysis)
save(list = ls(pattern = "^func_matrix_"), file = file.path(out_pref, "func_matrices.RData"))
if ("eggNOG_OGs" %in% names(final)) {
  save(func_matrix_eggNOG_OGs, file = file.path(out_pref, "func_eggNOG_OGs.RData"))
}

message("[make_func_matrices] Done. Wrote func_matrices.RData / func_eggNOG_OGs.RData and QC.")
