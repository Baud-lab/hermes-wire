#!/usr/bin/env bash
# Fast mixed model GWAS with GCTA --fastGWA using a sparse GRM.
# - Harmonises IDs across residuals (host_subject_id) ↔ metadata$RFID ↔ PLINK .fam
# - Starts from the residuals sample set, then:
#     1) removes samples without genotypes (.fam),
#     2) removes samples missing dam or cage in metadata,
#   and enforces identical order everywhere (RFID order).
# - Writes a single qcovar matrix (dummy-coded: sex, cohort, batch, dam, cage)
#   and a .pheno per trait (species), excluding 's__Bacteroides_F' for Bacteroides.
# - Creates a sparse GRM if none is supplied.
# - Runs fastGWA-mlm per trait; collates .fastGWA to per-trait .tsv.gz.
#
# Notes & references:
# - fastGWA requires a sparse GRM via --grm-sparse (create with --make-bK-sparse). :contentReference[oaicite:1]{index=1}
# - PLINK BED/BIM/FAM is required (--bfile). :contentReference[oaicite:2]{index=2}
# - We dummy-code categorical covariates into a numeric design matrix for --qcovar.
# - If you already created a sparse GRM, pass --sgrm-prefix to skip building it.
# - This script assumes *continuous* traits (fastGWA-mlm). For binary, adapt to --fastGWA-mlm-binary.
#
# Exit on error/undefined pipe failures; print commands.
set -euo pipefail

say() { printf "[fastgwa] %s\n" "$*"; }

usage() {
  cat <<EOF
Usage: $(basename "$0") \
  --plink-base /path/to/plink_base_without_ext \
  --resid-rda   residuals_qned_counts.RData \
  --metadata    metadata_Shallow.txt \
  --outdir      results/phase2/fastgwa \
  [--sgrm-prefix /path/to/sparse_grm_prefix] \
  [--threads 4] \
  [--trait-prefix s__] \
  [--genera Prevotella,Bacteroides]

Required:
  --plink-base     PLINK prefix (expects .bed/.bim/.fam)
  --resid-rda      RData with residual matrix (species x samples *or* samples x species)
  --metadata       TSV with columns including RFID and host_subject_id (+ sex, cohort, batch, dam, cage)
  --outdir         Output directory

Optional:
  --sgrm-prefix    Prefix of an existing sparse GRM (expects .grm.sp/.grm.id); if omitted, we build one.
  --threads        Threads for GCTA (default 1)
  --trait-prefix   Trait name prefix (default 's__')
  --genera         Comma-separated genera to analyse (default 'Prevotella,Bacteroides')
EOF
}

# ---- Parse args
PLINK=""
RESID=""
META=""
OUT=""
SGRM=""
THREADS=1
GENERA="Prevotella,Bacteroides"
TPREF="s__"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plink-base)   PLINK="$2"; shift 2;;
    --resid-rda)    RESID="$2"; shift 2;;
    --metadata)     META="$2"; shift 2;;
    --outdir)       OUT="$2"; shift 2;;
    --sgrm-prefix)  SGRM="$2"; shift 2;;
    --threads)      THREADS="${2}"; shift 2;;
    --genera)       GENERA="${2}"; shift 2;;
    --trait-prefix) TPREF="${2}"; shift 2;;
    -h|--help)      usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "$PLINK" || -z "$RESID" || -z "$META" || -z "$OUT" ]] && { usage; exit 1; }
[[ -f "${PLINK}.bed" && -f "${PLINK}.bim" && -f "${PLINK}.fam" ]] || { echo "Missing PLINK files for prefix: $PLINK" >&2; exit 1; }

mkdir -p "$OUT"/{prep,logs,per_trait}

say "===== START fastGWA ====="
say "PLINK base         : $PLINK"
say "Residual RData     : $RESID"
say "Metadata           : $META"
say "Output dir         : $OUT"
say "Sparse GRM prefix  : ${SGRM:-'(will build)'}"
say "Threads            : $THREADS"
say "Genera             : $GENERA"
say "Trait prefix       : $TPREF"

# ---- R: harmonise samples, build covariates & per-trait phenos
Rscript - <<RSCRIPT
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

say <- function(...) cat("[fastgwa][prep]", sprintf(...), "\n")

plink_fam <- fread("${PLINK}.fam", header=FALSE)
setnames(plink_fam, c("FID","IID","PID","MID","SEX","PHENO"))
say("Read PLINK FAM: n=", nrow(plink_fam))

# metadata (robust column naming)
meta <- fread("${META}")
colnames(meta)[colnames(meta) == "RFID"] <- "sample"
setnames(meta, tolower(names(meta)))
# try to infer expected fields
nm <- names(meta)
pick <- function(cands) { x <- intersect(cands, nm); if (length(x)) x[1] else NA_character_ }
rfid_col  <- pick(c("sample","rfid","sample_id","id","rf_id"))
hsid_col  <- pick(c("host_subject_id","host_id","subject_id"))
sex_col   <- pick(c("sex","gender"))
coh_col   <- pick(c("cohort","batch_cohort","cohort_id"))
bat_col   <- pick(c("sequencing_batch","seq_batch","batch","plate","library_batch"))
dam_col   <- pick(c("dam","dam_id","mother","maternal_id"))
cage_col  <- pick(c("cage","cage_id","housing","cohousing"))

req <- c(rfid_col, hsid_col, sex_col, coh_col, bat_col, dam_col, cage_col)
if (any(is.na(req))) {
  stop(sprintf("Missing required metadata columns; detected -> RFID:%s HSID:%s sex:%s cohort:%s batch:%s dam:%s cage:%s",
               rfid_col, hsid_col, sex_col, coh_col, bat_col, dam_col, cage_col))
}

# Construct helper columns and standardize
meta[, RFID := as.character(get(rfid_col))]
meta[, host_subject_id := as.character(get(hsid_col))]
meta[, sex := as.factor(get(sex_col))]
meta[, cohort := as.factor(get(coh_col))]
meta[, sequencing_batch := as.factor(get(bat_col))]
meta[, dam := as.factor(get(dam_col))]
meta[, cage := as.factor(get(cage_col))]

say("Metadata rows: %d  unique RFID: %d", nrow(meta), uniqueN(meta$RFID))

# residuals: load from RData; find the matrix/data.frame object
env <- new.env()
load("${RESID}", envir=env)
residuals=residuals_qned_counts_objs[[length(residuals_qned_counts_objs)]]
resid <- as.matrix(residuals)

say("Residuals dims: %d x %d (rows x cols)", nrow(resid), ncol(resid))

# We expect *columns* to be sample IDs equal to metadata$host_subject_id (per user note).
# If rownames match host_subject_id better than colnames, transpose.
hsid <- meta$host_subject_id
match_cols <- sum(colnames(resid) %in% hsid)
match_rows <- sum(rownames(resid) %in% hsid)
if (is.na(match_cols)) match_cols <- 0
if (is.na(match_rows)) match_rows <- 0
if (match_rows > match_cols) {
  resid <- t(resid)
  say("Transposed residuals to make *columns* = samples.")
}

if (is.null(colnames(resid))) stop("Residuals lack column names (sample IDs).")
# Map host_subject_id -> RFID (prefix before first underscore)
rfid_from_hsid <- sub("_.*$", "", colnames(resid))
colnames(resid) <- rfid_from_hsid

# Determine traits (species) from rownames
if (is.null(rownames(resid))) stop("Residuals lack row names (traits/species).")
traits_all <- rownames(resid)

# Filter traits by genera and prefix; exclude Bacteroides_F when genus is Bacteroides
tpf <- "${TPREF}"
genera <- strsplit("${GENERA}", ",", fixed=TRUE)[[1]]
traits_keep <- character(0)
sel <- traits_all[grepl(paste0(tpf,g), traits_all)]
  if (g == "Bacteroides") {
    sel <- sel[!grepl("s__Bacteroides_F", sel)]
  }
  traits_keep <- c(traits_keep, sel)
}
traits_keep <- unique(traits_keep)
traits_keep=traits_keep[!grepl("s__Bacteroides_F",traits_keep)]
if (!length(traits_keep)) stop("No traits matched the genera/prefix filters.")
say("Traits kept: %d (first 6): %s", length(traits_keep), paste(head(traits_keep), collapse=", "))

# ------------- Sample harmonisation (hard rules) -------------
# Start from samples present in residuals matrix (columns)
rf_from_resid <- colnames(resid)
say("Start samples from residuals: %d", length(rf_from_resid))

# Keep only RFIDs that exist in PLINK fam
fam_ids <- unique(plink_fam$IID)
keep1 <- intersect(rf_from_resid, fam_ids)
drop1 <- setdiff(rf_from_resid, keep1)
say("Drop no-genotype samples: %d", length(drop1))

# Keep only RFIDs with non-missing dam & cage
meta_sub <- meta[RFID %in% keep1]
meta_sub <- meta_sub[!is.na(dam) & dam != "" & !is.na(cage) & cage != ""]
keep2 <- intersect(keep1, meta_sub$RFID)
drop2 <- setdiff(keep1, keep2)
say("Drop missing dam/cage: %d", length(drop2))

# Final intersection, ordered by PLINK .fam order (makes GCTA happy)
final_rfids <- plink_fam$IID[plink_fam$IID %in% keep2]
say("Final sample size: %d", length(final_rfids))
stopifnot(length(final_rfids) > 10) # sanity for model fitting

# Subset residuals to final samples and transpose to samples x traits for writing
resid2 <- resid[traits_keep, final_rfids, drop=FALSE]
resid2 <- t(resid2)  # rows = samples, cols = traits
stopifnot(nrow(resid2) == length(final_rfids))

# Build qcovar design matrix: dummy-code all categorical covariates (no intercept)
meta_final <- meta_sub[match(final_rfids, meta_sub$RFID)]
stopifnot(all(meta_final$RFID == final_rfids))
mm <- model.matrix(~ 0 + sex + cohort + sequencing_batch + dam + cage, data=meta_final)
say("qcovar columns: %d", ncol(mm))
qcovar <- data.table(FID=final_rfids, IID=final_rfids, mm)
fwrite(qcovar, file = file.path("${OUT}","prep","covar.qcovar"), sep="\t")

# Write per-trait .pheno (FID IID PHE)
ndone <- 0L
for (tr in colnames(resid2)) {
  phe <- data.table(FID=final_rfids, IID=final_rfids, PHE = as.numeric(resid2[, tr]))
  # If all NA or (nearly) constant, skip
  if (all(is.na(phe$PHE)) || length(unique(na.omit(phe$PHE))) < 3) {
    next
  }
  fn <- file.path("${OUT}", "prep", paste0(tr, ".pheno"))
  fwrite(phe, fn, sep="\t", col.names=TRUE)
  ndone <- ndone + 1L
}
say("Wrote phenotype files: %d", ndone)
if (ndone == 0L) stop("No valid phenotype files after QC.")

# Also write the list of traits we actually created
traits_written <- list.files(file.path("${OUT}","prep"), pattern="\\.pheno$", full.names=FALSE)
traits_written <- sub("\\.pheno$", "", traits_written)
fwrite(data.table(trait=traits_written), file.path("${OUT}","prep","traits.tsv"), sep="\t")
RSCRIPT

# ---- Build a sparse GRM (if needed)
if [[ -z "${SGRM}" ]]; then
  say "Sparse GRM not provided; building from PLINK using GCTA."
  # Step 1: dense binary GRM (small N=~856 — fine)
  gcta64 --bfile "$PLINK" --make-grm-bin --out "$OUT/prep/full_grm" >"$OUT/logs/make_grm.log" 2>&1
  # Step 2: threshold to sparse with cutoff 0.05 (typical for close relatives) :contentReference[oaicite:3]{index=3}
  gcta64 --grm "$OUT/prep/full_grm" --make-bK-sparse 0.05 --out "$OUT/prep/sparse_005" >"$OUT/logs/make_sparse.log" 2>&1
  SGRM="$OUT/prep/sparse_005"
else
  say "Using existing sparse GRM prefix: $SGRM"
fi
[[ -f "${SGRM}.grm.sp" && -f "${SGRM}.grm.id" ]] || { echo "Sparse GRM files not found for prefix: $SGRM" >&2; exit 1; }

# ---- Count traits
NTRAITS=$(wc -l < "$OUT/prep/traits.tsv" | tr -d '[:space:]')
say "Traits to analyse with fastGWA: $NTRAITS"
[[ "$NTRAITS" -gt 0 ]] || { echo "No traits." >&2; exit 1; }

# ---- Run fastGWA per trait
i=0
while read -r TR; do
  i=$((i+1))
  PHE="$OUT/prep/${TR}.pheno"
  say "[$i/$NTRAITS] fastGWA for trait: ${TR}"

  # GCTA fastGWA-mlm (continuous); covariates as qcovar dummy matrix
  # Output: *.fastGWA (tab-delimited)
  gcta64 \
    --fastGWA-mlm \
    --bfile "$PLINK" \
    --grm-sparse "$SGRM" \
    --pheno "$PHE" \
    --qcovar "$OUT/prep/covar.qcovar" \
    --threads "$THREADS" \
    --out "$OUT/per_trait/${TR}" \
    >"$OUT/logs/${TR}.log" 2>&1

  # Normalise column names and compress to .tsv.gz for downstream MAGMA
  if [[ -f "$OUT/per_trait/${TR}.fastGWA" ]]; then
    awk 'BEGIN{OFS="\t"} NR==1{print "SNP","CHR","BP","A1","A2","N","BETA","SE","P"; next}
                 {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' \
      "$OUT/per_trait/${TR}.fastGWA" | bgzip -c > "$OUT/per_trait/${TR}.tsv.gz"
    tabix -s 2 -b 3 -e 3 "$OUT/per_trait/${TR}.tsv.gz" || true
  else
    say "WARNING: No .fastGWA produced for ${TR} (see $OUT/logs/${TR}.log)"
  fi
done < <(cut -f1 "$OUT/prep/traits.tsv")

say "===== DONE fastGWA ====="
