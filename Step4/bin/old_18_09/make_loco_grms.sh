#!/usr/bin/env bash
# make_loco_grms.sh — LOCO GRM builder that also works when gcta64 is an AppImage and FUSE is unavailable.
# It will auto-switch to AppImage extract-and-run, or pre-extract once, then run the extracted binary.
#
# References:
# - AppImage + FUSE requirement and extract options: https://docs.appimage.org/user-guide/run-appimages.html
# - Workarounds w/o FUSE (--appimage-extract / --appimage-extract-and-run): AskUbuntu threads
# - GCTA GRM docs: --make-grm-bin produces .grm.bin/.grm.N.bin/.grm.id
#
set -euo pipefail

echo "=== /etc/os-release ==="; cat /etc/os-release || true
echo "=== which gcta64 ==="; which gcta64 || true
echo "=== file \$(which gcta64) ==="; file "$(which gcta64 || echo /usr/local/bin/gcta64)" || true

# ---------------------------- CLI defaults ----------------------------
BFILE=""              # PLINK prefix (bed/bim/fam)
KEEP_IDS=""           # optional keep file (1- or 2-column)
OUTDIR="."
GCTA_BIN="${GCTA_BIN:-gcta64}"  # overridable via env or --gcta_exec

usage() {
  cat <<'EOF'
[loco] make_loco_grms.sh
  --bfile     <prefix>   PLINK trio prefix (e.g., candidates_kept)
  --keep_ids  <path>     (optional) 1-col (IID) or 2-col (FID IID) keep file
  --outdir    <dir>      Output directory (default: .)
  --gcta_exec <path>     gcta64 executable (default: gcta64)
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile)     BFILE="$2"; shift 2;;
    --keep_ids)  KEEP_IDS="$2"; shift 2;;
    --outdir)    OUTDIR="$2"; shift 2;;
    --gcta_exec) GCTA_BIN="$2"; shift 2;;
    -h|--help)   usage;;
    *) echo "[loco] Unknown arg: $1"; usage;;
  esac
done

echo "[loco] START"
echo "[loco] bfile   : ${BFILE}"
echo "[loco] keep_ids: ${KEEP_IDS:-<none>}"
echo "[loco] outdir  : ${OUTDIR}"
echo "[loco] gcta bin: ${GCTA_BIN}"

# ---------------------------- Sanity checks ----------------------------
[[ -n "$BFILE" && -s "${BFILE}.bed" && -s "${BFILE}.bim" && -s "${BFILE}.fam" ]] \
  || { echo "[loco] ERROR: missing PLINK files for --bfile '${BFILE}'"; exit 2; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ---------------------------- Resolve working GCTA runner ----------------------------
# We define a RUN_GCTA() helper that always runs a working gcta64, regardless of AppImage/FUSE.
# Strategy:
#   1) Try plain binary: "$GCTA_BIN" --help (fast path)
#   2) If we see AppImage/FUSE errors, prefer: "$GCTA_BIN --appimage-extract-and-run ..."
#   3) (Optional optimization) If extract-and-run is slow, pre-extract once:
#        "$GCTA_BIN --appimage-extract"  -> "squashfs-root/AppRun" or ".../usr/bin/gcta64"
#
GCTA_PATH="$(command -v "$GCTA_BIN" || true)"
[[ -n "$GCTA_PATH" ]] || { echo "[loco] ERROR: cannot find '$GCTA_BIN' in PATH"; exit 127; }

echo "[loco] Probing gcta64: ${GCTA_PATH}"
if "$GCTA_PATH" --help >/dev/null 2>&1; then
  echo "[loco] gcta64 appears runnable directly."
  RUN_GCTA() { "$GCTA_PATH" "$@"; }
else
  echo "[loco] Direct run failed; checking for AppImage options (FUSE-less path)..."
  # Does it respond to --appimage-help?
  if "$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
    echo "[loco] AppImage detected. Using --appimage-extract-and-run to bypass FUSE."
    # Fast path: extract-and-run per invocation
    RUN_GCTA() { "$GCTA_PATH" --appimage-extract-and-run "$@"; }

    # If you prefer to extract ONCE and reuse (uncomment the block below):
    : <<'OPTIONAL_PREEXTRACT'
    echo "[loco] Pre-extracting AppImage once (creates ./squashfs-root)..."
    "$GCTA_PATH" --appimage-extract >/dev/null 2>&1 || true
    if [[ -x squashfs-root/AppRun ]]; then
      GCTA_REAL="$(pwd)/squashfs-root/AppRun"
    elif [[ -x squashfs-root/usr/bin/gcta64 ]]; then
      GCTA_REAL="$(pwd)/squashfs-root/usr/bin/gcta64"
    else
      echo "[loco] WARNING: could not locate extracted gcta binary; falling back to --appimage-extract-and-run."
      GCTA_REAL=""
    fi
    if [[ -n "${GCTA_REAL}" ]]; then
      echo "[loco] Using extracted runner: ${GCTA_REAL}"
      RUN_GCTA() { "$GCTA_REAL" "$@"; }
    fi
OPTIONAL_PREEXTRACT

  else
    echo "[loco] ERROR: gcta64 is not runnable and does not expose AppImage helpers."
    echo "[loco] Please install libfuse (host) OR use a static gcta64 binary (recommended)."
    exit 127
  fi
fi

# ---------------------------- Build keep file if provided ----------------------------
if [[ -n "${KEEP_IDS}" ]]; then
  echo "[loco] Building two-column keep file against ${BFILE}.fam"
  # Detect 1-col vs 2-col input
  ncols=$(awk '{print NF; exit}' "${KEEP_IDS}")
  if [[ "${ncols:-0}" -ge 2 ]]; then
    awk 'NF>=2{print $1, $2}' "${KEEP_IDS}" > .keep.2col.txt
  else
    awk 'NR==FNR{ids[$1]=1; next} ($2 in ids){print $1, $2}' "${KEEP_IDS}" "${BFILE}.fam" > .keep.2col.txt
  fi
  nkeep=$(awk 'NF>0{c++} END{print c+0}' .keep.2col.txt)
  echo "[loco] keep matched to FAM: ${nkeep}"
  [[ "$nkeep" -gt 0 ]] || { echo "[loco] ERROR: 0 matches between keep_ids and FAM"; exit 3; }
fi

# ---------------------------- Chromosome set (numeric only) ----------------------------
awk '$1 ~ /^[0-9]+$/{print $1}' "${BFILE}.bim" | sort -n | uniq > .chr.list
echo "[loco] chromosomes found:"; cat .chr.list

# ---------------------------- LOCO GRMs ----------------------------
while read -r CHR; do
  [[ -n "$CHR" ]] || continue
  echo "[loco] chr=${CHR} → compute GRM on all SNPs except this chromosome"
  awk -v C="$CHR" '$1 ~ /^[0-9]+$/ && $1 != C {print $2}' "${BFILE}.bim" > ".include_chr${CHR}.snps"
  nsnps=$(wc -l < ".include_chr${CHR}.snps")
  if [[ "$nsnps" -eq 0 ]]; then
    echo "[loco] WARNING: include set empty for chr ${CHR}; skipping"
    continue
  fi

  if [[ -s .keep.2col.txt ]]; then
    RUN_GCTA --bfile "${BFILE}" --keep .keep.2col.txt --extract ".include_chr${CHR}.snps" --make-grm-bin --out "loco_chr${CHR}"
  else
    RUN_GCTA --bfile "${BFILE}" --extract ".include_chr${CHR}.snps" --make-grm-bin --out "loco_chr${CHR}"
  fi
done < .chr.list

# ---------------------------- Index + RDS helper ----------------------------
printf "chr\tprefix\n" > grm_loco_index.tsv
while read -r CHR; do
  if [[ -s "loco_chr${CHR}.grm.bin" ]]; then
    printf "%s\tloco_chr%s\n" "$CHR" "$CHR" >> grm_loco_index.tsv
  fi
done < .chr.list

if command -v Rscript >/dev/null 2>&1; then
  Rscript - <<'RS'
idx <- read.table("grm_loco_index.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!all(c("chr","prefix") %in% names(idx))) stop("grm_loco_index.tsv must have chr and prefix columns")
vec <- setNames(idx$prefix, idx$chr)
saveRDS(vec, file="grm_loco.rds")
RS
  echo "[loco] wrote grm_loco.rds"
else
  echo "[loco] NOTE: Rscript not found; skipping grm_loco.rds (use grm_loco_index.tsv instead)."
fi

echo "[loco] DONE"
