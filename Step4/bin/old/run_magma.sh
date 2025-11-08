#!/usr/bin/env bash
# MAGMA gene-level tests + (optionally) genus-level family selection
# - Reads per-trait association .tsv.gz from --assoc_dir
# - Converts to MAGMA .pval using a robust (CHR,BP)→BIM.SNP join
# - De-duplicates rows per SNP (keeps min P, max N) to be LOCO-safe
# - Annotates with windows and runs gene analysis (with optional per-genus auto-tune)
# - Writes helper tables for downstream (magma.snps.loc, gene_name_map.tsv)
#
# Exit on error, treat unset vars as error, fail on pipe errors
set -euo pipefail
shopt -s nullglob

echo "[magma] START"

# -------------------- CLI --------------------
PLINK=""; ADIR=""; GLOC=""; GENES_XLSX=""; Q="0.10"; OUT="."
WUP="0"; WDN="0"; FDR_SCOPE="genus"
KEEP_IDS=""; GENUS_COMBINER="empirical_brown"; GENE_MODEL="multi"; AUTO_TUNE="false"

usage() {
  cat >&2 <<USAGE
Usage: $0 \\
  --plink_base candidates \\
  --assoc_dir  <assoc dir with *.tsv.gz> \\
  --gene_loc   magma.genes.loc \\
  --genes_xlsx <panel.xlsx> \\
  --q          0.10 \\
  --outdir     . \\
  --window_up  50000 --window_down 50000 \\
  --fdr_scope  genus \\
  --gene_model multi \\
  --auto_tune  false \\
  [--keep_ids  <IID list to restrict LD panel>]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --plink_base)      PLINK=$2; shift 2;;
    --assoc_dir)       ADIR=$2; shift 2;;
    --gene_loc)        GLOC=$2; shift 2;;
    --genes_xlsx)      GENES_XLSX=$2; shift 2;;
    --q)               Q=$2; shift 2;;
    --outdir)          OUT=$2; shift 2;;
    --window_up)       WUP=$2; shift 2;;
    --window_down)     WDN=$2; shift 2;;
    --fdr_scope)       FDR_SCOPE=$2; shift 2;;
    --keep_ids)        KEEP_IDS=$2; shift 2;;
    --genus_combiner)  GENUS_COMBINER=$2; shift 2;;
    --gene_model)      GENE_MODEL=$2; shift 2;;
    --auto_tune)       AUTO_TUNE=$2; shift 2;;
    *) echo "Unknown arg $1" >&2; usage;;
  esac
done

[[ -n "$PLINK" && -n "$ADIR" && -n "$GLOC" && -n "$OUT" ]] || usage
mkdir -p "$OUT" "$OUT/family_fdr_by_genus" "$OUT/species_calls_by_genus"

normalize_gene_model() {
  # Map friendly names or auto-tune outputs to MAGMA's exact CLI modifiers.
  # Valid returns for --pval mode: snp-wise=mean | snp-wise=top | snp-wise=top,1 | multi
  local m="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
  case "$m" in
    top)      echo "snp-wise=top" ;;
    top1)     echo "snp-wise=top,1" ;;
    mean)     echo "snp-wise=mean" ;;
    multi|multi=snp-wise|multi=all) echo "multi" ;;
    snp-wise=*|snpwise=*) echo "$m" ;;
    *)        echo "snp-wise=mean" ;;
  esac
}

# -------------------- Tool detection --------------------
PLINK_BIN="$(command -v plink || true)"; [[ -z "$PLINK_BIN" ]] && PLINK_BIN="$(command -v plink2 || true)"
MAGMA_BIN="$(command -v magma || true)"
[[ -z "$PLINK_BIN" ]] && { echo "[magma] ERROR: plink/plink2 not found in PATH" >&2; exit 1; }
[[ -z "$MAGMA_BIN" ]] && { echo "[magma] ERROR: magma not found in PATH" >&2; exit 1; }

echo "[magma] LD reference base: $PLINK"
echo "[magma] Assoc dir: $ADIR"
echo "[magma] Gene loc: $GLOC"
[[ -n "$GENES_XLSX" ]] && echo "[magma] Gene panel: $GENES_XLSX"
echo "[magma] Window(bp): upstream=${WUP} downstream=${WDN} (MAGMA uses kb; will pass: $((WUP/1000)),$((WDN/1000)))"
echo "[magma] FDR scope: ${FDR_SCOPE}"
echo "[magma] Gene model default: ${GENE_MODEL}"
echo "[magma] Auto-tune per genus: ${AUTO_TUNE}"
[[ -n "$KEEP_IDS" ]] && echo "[magma] Keep IDs file: $(basename "$KEEP_IDS")"

# -------------------- 0) gene_loc sanity + de-dup by GENE_ID ----------
if [[ ! -s "$GLOC" ]]; then
  echo "[magma] ERROR: gene_loc file missing/empty: $GLOC" >&2; exit 1
fi
first_fields=$(awk 'NF>0 && $1 !~ /^#/ { print NF; exit }' "$GLOC" || echo 0)
if [[ "${first_fields:-0}" -lt 4 ]]; then
  echo "[magma] ERROR: gene_loc has ${first_fields:-0} columns; need >=4 (GENE_ID CHR START END)"; exit 1
fi
echo "[magma] Checking duplicate GENE_IDs in gene_loc ..."
dup_cnt=$(awk '++seen[$1]==2{c++} END{print c+0}' "$GLOC")
if [[ "$dup_cnt" -gt 0 ]]; then
  echo "[magma] WARNING: ${dup_cnt} duplicate GENE_ID rows; keeping first occurrence per ID"
  awk 'BEGIN{FS=OFS="\t"} !seen[$1]++' "$GLOC" > "$OUT/.genes.loc.nodup"
  mv "$OUT/.genes.loc.nodup" "$GLOC"
fi
uniq_cnt=$(awk '{print $1}' "$GLOC" | sort -u | wc -l | awk '{print $1}')
tot_cnt=$(wc -l < "$GLOC")
echo "[magma] gene_loc rows: ${tot_cnt} unique IDs: ${uniq_cnt}"

# -------------------- 1) Normalize PLINK → variant-major --------------------
echo "[magma] Normalizing PLINK bed to variant-major with: $PLINK_BIN"
"$PLINK_BIN" --bfile "$PLINK" --make-bed --out "$OUT/magma_vm" >/dev/null 2>&1
BFILE="$OUT/magma_vm"

# -------------------- 1b) Restrict LD panel strictly to --keep_ids --------
if [[ -n "$KEEP_IDS" ]]; then
  tr -d '\r' < "${BFILE}.fam" > "$OUT/.fam_raw"
  tr -d '\r' < "$KEEP_IDS" > "$OUT/.keep_raw"
  ncol_keep=$(awk 'NF>0{print NF; exit}' "$OUT/.keep_raw")
  n_in=$(awk 'NF>0{c++} END{print c+0}' "$OUT/.keep_raw")
  echo "[magma] keep_ids lines: ${n_in} (detected ${ncol_keep}-column format)"
  if [[ "${ncol_keep:-1}" -ge 2 ]]; then
    awk 'NR==FNR{key[$1"\t"$2]=1; next} ( ($1"\t"$2) in key ){print $1,$2}' \
      "$OUT/.keep_raw" "$OUT/.fam_raw" > "$OUT/.keep_ids_2col.txt"
  else
    awk 'NR==FNR{ids[$1]=1; next} ($2 in ids){print $1,$2}' \
      "$OUT/.keep_raw" "$OUT/.fam_raw" > "$OUT/.keep_ids_2col.txt"
  fi
  NKEEP=$(awk 'NF>0{c++} END{print c+0}' "$OUT/.keep_ids_2col.txt")
  echo "[magma] matched to FAM: ${NKEEP}"
  [[ "$NKEEP" -gt 0 ]] || { echo "[magma] ERROR: 0 matches between keep_ids and .fam"; exit 1; }
  "$PLINK_BIN" --bfile "$BFILE" --keep "$OUT/.keep_ids_2col.txt" --make-bed --out "$OUT/magma_vm_kept" >/dev/null 2>&1
  BFILE="$OUT/magma_vm_kept"
fi

# -------------------- 2) SNP location table for MAGMA ----------------------
awk 'BEGIN{OFS="\t"} {print $2, $1, $4}' "${BFILE}.bim" > "$OUT/magma.snps.loc"
echo "[magma] Wrote SNP loc: $OUT/magma.snps.loc nSNPs=$(wc -l < "$OUT/magma.snps.loc")"

# -------------------- 3) Build gene_name_map.tsv (optional Excel) ---------
if [[ -n "$GENES_XLSX" && -s "$GENES_XLSX" ]]; then
  Rscript -e "
    suppressPackageStartupMessages({ library(readxl); library(data.table) })
    gl <- fread('${GLOC}', header=FALSE, fill=TRUE, data.table=TRUE)
    setnames(gl, 1:4, c('GENE_ID','CHR','START','END'))
    xl <- as.data.table(read_excel('${GENES_XLSX}'))
    req <- c('Bioprocess','Gene','Chromosome','Start','End')
    if (!all(req %in% names(xl))) stop('Expected Excel headers not found.')
    setnames(xl, c('Gene','Chromosome','Start','End'), c('GENE','CHR','START','END'))
    xl <- unique(xl[, .(CHR,START,END,GENE)])
    m <- merge(gl[, .(GENE_ID, CHR, START, END)],
               xl, by=c('CHR','START','END'), all.x=TRUE, allow.cartesian=FALSE)
    m[is.na(GENE) | trimws(GENE)=='', GENE := GENE_ID]
    data.table::fwrite(unique(m[, .(GENE_ID, GENE)]), file='${OUT}/gene_name_map.tsv', sep='\t')
  "
else
  awk 'BEGIN{OFS="\t"}{print $1,$1}' "$GLOC" | sed '1iGENE_ID\tGENE' > "${OUT}/gene_name_map.tsv"
fi

# -------------------- 4) Build MAGMA .pval files per trait -----------------
# Important: use a single-quoted heredoc so Bash does not expand $ inside the R code
BIM_PATH="${BFILE}.bim"
for f in "$ADIR"/*.tsv.gz; do
  t=$(basename "$f" .tsv.gz)
  echo "[magma] Trait: $t -> .pval via (CHR,BP) ⟶ BIM.SNP and N from assoc (LOCO-safe de-dup)"
  ASSOC="$f" TRAIT="$t" OUTDIR="$OUT" BIM="$BIM_PATH" Rscript - <<'RS'
suppressPackageStartupMessages({ library(data.table) })
bim <- fread(Sys.getenv("BIM"), col.names=c("CHR","SNP","CM","BP","A1","A2"))
df  <- suppressWarnings(fread(Sys.getenv("ASSOC")))
setnames(df, tolower(names(df)))

# detect columns (case-insensitive)
chr_col <- intersect(c('chr','chrom','chromosome'), names(df))[1]
pos_col <- intersect(c('pos','position','bp'), names(df))[1]
pcol    <- intersect(c('score.pval','pval','p','p.value','pvalue'), names(df))[1]
ncol    <- intersect(c('n.obs','nobs','n','neff','nca'), names(df))[1]
if (is.na(chr_col) || is.na(pos_col) || is.na(pcol)) stop('Missing chr/pos/p columns in ', Sys.getenv("ASSOC"))

# normalize columns
df[, CHR := as.character(get(chr_col))]
df[, BP  := as.integer(get(pos_col))]
df[, P   := pmax(as.numeric(get(pcol)), .Machine$double.xmin)]
if (!is.na(ncol)) {
  df[, N := as.integer(round(get(ncol)))]
} else {
  df[, N := NA_integer_]
}

# sanitize
df <- df[is.finite(P) & P <= 1 & is.finite(BP) & !is.na(CHR)]
# join with BIM on (CHR,BP)
bim[, CHR := as.character(CHR) ]
m <- merge(df[, .(CHR,BP,P,N)], bim[, .(CHR,BP,SNP)], by=c('CHR','BP'), all=FALSE)

# LOCO/dup safety: collapse by SNP (min P, max N)
if (nrow(m) == 0L) stop('After join, no SNPs matched between assoc and BIM (CHR/BP).')
m <- m[, .(P = min(P, na.rm=TRUE), N = max(N, na.rm=TRUE)), by=.(SNP)]

# If all N are NA (some outputs omit N), drop N and let MAGMA default
if (all(is.na(m$N))) {
  fwrite(unique(m[, .(SNP,P)]),
         file=file.path(Sys.getenv("OUTDIR"), paste0(Sys.getenv("TRAIT"), ".pval")),
         sep='\t', quote=FALSE, col.names=TRUE)
} else {
  fwrite(unique(m[, .(SNP,P,N)]),
         file=file.path(Sys.getenv("OUTDIR"), paste0(Sys.getenv("TRAIT"), ".pval")),
         sep='\t', quote=FALSE, col.names=TRUE)
}
RS
done

# -------------------- 5) Auto-tune per genus (optional) --------------------
AUTO_FILE="$OUT/auto_choices.tsv"
if [[ "${AUTO_TUNE}" == "true" ]]; then
  echo "[auto] Pre-analyzing per-genus structure to pick gene_model & combiner ..."
  python3 - "$OUT" "$GENE_MODEL" "$GENUS_COMBINER" "$GLOC" <<'PY'
import sys, os, glob, re, numpy as np, pandas as pd
outdir, def_model, def_combiner, gloc = sys.argv[1], sys.argv[2].lower(), sys.argv[3].lower(), sys.argv[4]
pval_paths = sorted(glob.glob(os.path.join(outdir, "*.pval")))
def trait_to_genus(t):
    t = re.sub(r"_assoc$", "", t)
    m = re.search(r"__([^_]+)", t)
    return m.group(1) if m else t.split("_")[0]
rows=[]
for p in pval_paths:
    tr = os.path.basename(p).replace(".pval","")
    try:
        df = pd.read_csv(p, sep="\t")
        if {"SNP","P"}.issubset(df.columns) and len(df)>0:
            df = df[np.isfinite(df["P"])]
            df["trait"] = tr
            df["genus"] = trait_to_genus(tr)
            rows.append(df[["SNP","P","trait","genus"]])
    except Exception:
        pass
if not rows:
    open(os.path.join(outdir, "auto_choices.tsv"), "w").write("genus\tgene_model\tcombiner\n")
    sys.exit(0)
D = pd.concat(rows, ignore_index=True)
choices=[]
for gn, sub in D.groupby("genus", sort=False):
    traits = sorted(sub["trait"].unique().tolist())
    ntraits = len(traits)
    wide = sub.pivot_table(index="SNP", columns="trait", values="P", aggfunc="min")
    sim = []
    W = -np.log10(wide.replace(0, np.nextafter(0,1)))
    W = W.replace([np.inf,-np.inf], np.nan).dropna(axis=0, how="any")
    mean_r = 0.0
    if W.shape[1] >= 2:
        for i in range(W.shape[1]):
            for j in range(i+1, W.shape[1]):
                a, b = W.iloc[:,i], W.iloc[:,j]
                r = a.corr(b, method="spearman")
                if np.isfinite(r): sim.append(r)
        mean_r = float(np.nanmean(sim)) if sim else 0.0
    # simple heuristic
    gene_model = "multi" if (ntraits >= 6 and mean_r >= 0.3) else ("top" if mean_r < 0.2 else def_model)
    combiner = "simes" if (ntraits < 2 or mean_r < 0.05) else ("empirical_brown" if def_combiner!="simes" else "simes")
    print(f"[auto] Genus {gn}: gene_model={gene_model} combiner={combiner} (traits={ntraits}, mean_r={mean_r:.3f})")
    choices.append({"genus":gn, "gene_model":gene_model, "combiner":combiner})
pd.DataFrame(choices).to_csv(os.path.join(outdir, "auto_choices.tsv"), sep="\t", index=False)
PY
else
  echo -e "genus\tgene_model\tcombiner" > "$AUTO_FILE"
fi

# -------------------- 6) MAGMA annotate + gene analysis --------------------
win_kb_u=$((WUP/1000))
win_kb_d=$((WDN/1000))

# Read auto choices into assoc arrays
declare -A CHOICE_MODEL
declare -A CHOICE_COMB
if [[ -s "$AUTO_FILE" ]]; then
  while IFS=$'\t' read -r g m c rest; do
    [[ "$g" == "genus" ]] && continue
    CHOICE_MODEL["$g"]="$m"
    CHOICE_COMB["$g"]="$c"
  done < "$AUTO_FILE"
fi

trait_to_genus() {
  local t="$1"
  t="${t%_assoc}"
  if [[ "$t" =~ __([^_]+) ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "${t%%_*}"
  fi
}

for t in "$OUT"/*.pval; do
  trait=$(basename "$t" .pval)
  gn=$(trait_to_genus "$trait")
  gm="${CHOICE_MODEL[$gn]:-$GENE_MODEL}"
  comb="${CHOICE_COMB[$gn]:-$GENUS_COMBINER}"

  echo "[magma] MAGMA annotate (nonhuman, window=${win_kb_u},${win_kb_d} kb) for ${trait} (genus=${gn})"
  "$MAGMA_BIN" --annotate window="${win_kb_u},${win_kb_d}" nonhuman \
    --snp-loc "$OUT/magma.snps.loc" \
    --gene-loc "$GLOC" \
    --out "$OUT/${trait}"

  GM_CLI="$(normalize_gene_model "$gm")"
  echo "[magma] MAGMA gene analysis for ${trait} [gene_model=${gm} -> ${GM_CLI}]"
  if grep -q -i $'\tN$' "$t" 2>/dev/null; then
    "$MAGMA_BIN" --pval "$t" ncol=N \
      --gene-annot "$OUT/${trait}.genes.annot" \
      --bfile "$BFILE" \
      --gene-model "$GM_CLI" \
      --out "$OUT/${trait}.genes"
  else
    "$MAGMA_BIN" --pval "$t" \
      --gene-annot "$OUT/${trait}.genes.annot" \
      --bfile "$BFILE" \
      --gene-model "$GM_CLI" \
      --out "$OUT/${trait}.genes"
  fi
done

echo "[magma] DONE"
