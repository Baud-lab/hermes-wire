#!/usr/bin/env bash
set -euo pipefail

usage(){ echo "Usage: $0 --grm_prefix <pfx> --pheno gcta.pheno --trait_map trait_map.tsv --threads N --out_prefix cv__"; exit 1; }
GRM= PHENO= TMAP= THREADS=1 OUT=cv__
while [[ $# -gt 0 ]]; do
  case "$1" in
    --grm_prefix) GRM="$2"; shift 2;;
    --pheno)      PHENO="$2"; shift 2;;
    --trait_map)  TMAP="$2"; shift 2;;
    --threads)    THREADS="$2"; shift 2;;
    --out_prefix) OUT="$2"; shift 2;;
    --gcta_exec) GCTA_BIN="$2"; shift 2;;
    *) echo "Unknown arg $1"; usage;;
  esac
done
[[ -f "${PHENO}" && -f "${TMAP}" ]] || usage

# Build GRM argument (prefix vs list-of-prefixes)
if [[ -f "$GRM" && "${GRM##*.}" == "txt" ]]; then
  GRM_ARG=(--mgrm-gz "$GRM")
else
  GRM_ARG=(--grm "$GRM")
fi

GCTA_PATH="$(command -v ${GCTA_BIN} || true)"
  [[ -n "$GCTA_PATH" ]] || { echo "[mgrm] ERROR: cannot find ${GCTA_BIN}"; exit 127; }
  if "$GCTA_PATH" --help >/dev/null 2>&1; then
    RUN_GCTA() { "$GCTA_PATH" "$@"; }
  elif "$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
    RUN_GCTA() { "$GCTA_PATH" --appimage-extract-and-run "$@"; }
  else
    echo "[mgrm] ERROR: gcta64 not runnable"; exit 127
  fi


# For each phenotype column (mpheno), run cvBLUP
tail -n +2 "${TMAP}" | while IFS=$'\t' read -r trait mpheno; do
  tag="${trait//[^A-Za-z0-9_]/_}"
  RUN_GCTA "${GRM_ARG[@]}" \
         --pheno "${PHENO}" \
         --mpheno "${mpheno}" \
         --reml --cvblup \
         --reml-no-constrain \
         --thread-num "${THREADS}" \
         --out "${OUT}${tag}" || true
  # GCTA writes *.indi.blp with per-individual cvBLUPs (and residuals)
done
