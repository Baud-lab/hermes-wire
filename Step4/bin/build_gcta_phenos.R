#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(jsonlite)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--fam",          type="character"),
  make_option("--scan_lookup",  type="character"),
  make_option("--resid_rda",    type="character"),
  make_option("--beta_rda",     type="character"),
  make_option("--alpha_rda",    type="character"),
  make_option("--clusters_rda", type="character"),
  make_option("--phen_rda",             type="character"),
  make_option("--phen_glucose_col",     type="character", default="Glucose"),
  make_option("--phen_bmi_col",         type="character", default="BMI"),
  make_option("--genus",        type="character"),
  make_option("--rg_traits",    type="character",
              default="beta__PD_PC1,alpha__PD_q2,cluster__3"),
  make_option("--rg_pairs",     type="character",
              default="genus~beta__PD_PC1,genus~alpha__PD_q2,genus~cluster__3,genus~other_genus")
)))

# ---------- helpers ----------
norm_id <- function(x) sub("_.*$", "", as.character(x))

map_token_to_trait <- function(tok, opt) {
  t <- trimws(as.character(tok))
  tl <- tolower(t)
  if (tl %in% c("Glucose","glucose","glu")) return(opt$phen_glucose_col)
  if (tl %in% c("BMI","bmi","body_mass_index","bodymassindex")) return(opt$phen_bmi_col)
  t
}

split_csv <- function(x) unique(trimws(unlist(strsplit(as.character(x), ",", fixed=TRUE))))

# ---------- IDs: FAM -> lookup labels ----------
FAM <- fread(opt$fam, header = FALSE)
setnames(FAM, c("FID","IID","PAT","MAT","SEX","PHENO"))
FAM[, IID := as.character(IID)]
PLINK_IDS <- FAM[, .(FID, IID)]

LU <- fread(opt$scan_lookup)
stopifnot(all(c("scan_id","sample_label") %in% names(LU)))
LU[, scan_id := as.character(scan_id)]
LU[, label   := norm_id(sample_label)]
IDS <- merge(PLINK_IDS, LU[, .(scan_id, label)], by.x="IID", by.y="scan_id", all.x=TRUE)

# ---------- Load residual layers ----------
load(opt$resid_rda)
R_taxa <- NULL
if (exists("residuals_qned_counts_objs")) {
  idx <- max(1L, length(residuals_qned_counts_objs) - 1L)
  R_taxa <- residuals_qned_counts_objs[[idx]]
}
load(opt$beta_rda);    R_beta    <- residuals_qned_counts_beta_objs[[1]]
load(opt$alpha_rda);   R_alpha   <- residuals_qned_counts_alpha_objs[[1]]
load(opt$clusters_rda);R_cluster <- residuals_qned_counts_clusters_objs[[1]]

# ---------- Host phenotypes ----------
load(opt$phen_rda)
if (!exists("df") || !("samples" %in% names(df)))
  stop("--phen_rda must contain object 'df' with a 'samples' column")
need_cols <- c(opt$phen_glucose_col, opt$phen_bmi_col)
miss <- setdiff(need_cols, names(df))
if (length(miss)) stop("--phen_rda$df missing columns: ", paste(miss, collapse=", "))
glucose_vec <- setNames(as.numeric(df[[opt$phen_glucose_col]]), norm_id(df$samples))
bmi_vec     <- setNames(as.numeric(df[[opt$phen_bmi_col]]),     norm_id(df$samples))

# ---------- projectors ----------
project_vec <- function(M, rowkey, ids = IDS) {
  if (is.null(M)) return(rep(NA_real_, nrow(ids)))
  rn <- rownames(M)
  hit <- which(rn == rowkey); if (!length(hit)) hit <- which(startsWith(rn, rowkey))
  if (!length(hit)) return(rep(NA_real_, nrow(ids)))
  v <- as.numeric(M[hit[1], , drop=TRUE])
  cn <- colnames(M)
  labmap <- setNames(seq_along(cn), norm_id(cn))
  out <- rep(NA_real_, nrow(ids))
  ok  <- !is.na(ids$label) & ids$label %in% names(labmap)
  out[ok] <- v[ labmap[ ids$label[ok] ] ]
  out
}
project_host <- function(v, ids = IDS) {
  out <- rep(NA_real_, nrow(ids))
  ok  <- !is.na(ids$label) & ids$label %in% names(v)
  out[ok] <- as.numeric(v[ ids$label[ok] ])
  out
}

# ---------- genus names ----------
g_this  <- paste0("g__", opt$genus)
g_other <- if (grepl("Bacteroides", opt$genus, ignore.case=TRUE)) "g__Prevotella" else "g__Bacteroides"

# ---------- build phenotype table ----------
PH <- IDS[, .(FID, IID)]
vecs <- list()
if (!is.null(R_taxa)) {
  vecs[[g_this]] <- project_vec(R_taxa, g_this)
  vo <- project_vec(R_taxa, g_other); if (any(!is.na(vo))) vecs[[g_other]] <- vo
}
vecs[["beta__PD_PC1"]] <- project_vec(R_beta,    "beta__PD_PC1")
vecs[["alpha__PD_q2"]] <- project_vec(R_alpha,   "alpha__PD_q2")
vecs[["cluster__3"]]   <- project_vec(R_cluster, "cluster__3")
vecs[[opt$phen_glucose_col]] <- project_host(glucose_vec)
vecs[[opt$phen_bmi_col]]     <- project_host(bmi_vec)

# ---------- decide which traits to include ----------
# 1) from --rg_traits
want_from_traits <- vapply(split_csv(opt$rg_traits), map_token_to_trait, character(1), opt=opt)

# 2) also include anything referenced in --rg_pairs (both sides of each "~")
pairs_tokens <- split_csv(opt$rg_pairs)
want_from_pairs <- unlist(lapply(pairs_tokens, function(s) {
  ab <- unlist(strsplit(s, "~", fixed=TRUE))
  if (length(ab) != 2) return(character(0))
  c(ab[1], ab[2])
}), use.names = FALSE)
want_from_pairs <- vapply(want_from_pairs, function(tok) {
  t <- trimws(tok)
  if (t == "genus") return(g_this)
  if (t == "other_genus") return(g_other)
  map_token_to_trait(t, opt)
}, character(1))

# union → these are the columns we’ll try to add
want_all <- unique(c(want_from_traits, want_from_pairs))

always <- names(vecs)[grepl("^g__", names(vecs))]           # keep genus columns when present
keep   <- unique(c(always, intersect(want_all, names(vecs))))  # add all referenced traits
for (nm in keep) PH[, (nm) := vecs[[nm]]]

# drop all-NA traits
ph_cols <- setdiff(names(PH), c("FID","IID"))
allna <- ph_cols[sapply(ph_cols, function(nm) all(is.na(PH[[nm]])))]
if (length(allna)) PH[, (allna) := NULL]
ph_cols <- setdiff(names(PH), c("FID","IID"))
stopifnot(length(ph_cols) > 0)

cat("Included phenotypes:\n"); print(ph_cols)
cat("Non-missing per phenotype:\n"); print(sapply(ph_cols, function(nm) sum(!is.na(PH[[nm]]))))

# ---------- write outputs ----------
PH_OUT <- copy(PH)
for (nm in ph_cols) set(PH_OUT, which(is.na(PH_OUT[[nm]])), nm, -9)
fwrite(PH_OUT, file="gcta.pheno", sep="\t", col.names=FALSE, quote=FALSE)

trait_map <- data.table(trait = ph_cols, mpheno = seq_along(ph_cols))
fwrite(trait_map, file="trait_map.tsv", sep="\t")

resolve_token <- function(tok) {
  t <- trimws(tok)
  if (t == "genus")       return(g_this)
  if (t == "other_genus") return(if (g_other %in% ph_cols) g_other else NA_character_)
  map_token_to_trait(t, opt)
}

pairs_tbl <- rbindlist(lapply(pairs_tokens, function(s) {
  bits <- unlist(strsplit(s, "~", fixed=TRUE))
  if (length(bits)!=2) return(NULL)
  a <- resolve_token(bits[1]); b <- resolve_token(bits[2])
  if (!all(c(a,b) %in% ph_cols)) return(NULL)
  data.table(
    trait1 = a, trait2 = b,
    col1   = trait_map$mpheno[match(a, trait_map$trait)],
    col2   = trait_map$mpheno[match(b, trait_map$trait)]
  )
}), fill=TRUE)
pairs_tbl <- unique(pairs_tbl[complete.cases(pairs_tbl)])
fwrite(pairs_tbl, file="pairs.tsv", sep="\t")
