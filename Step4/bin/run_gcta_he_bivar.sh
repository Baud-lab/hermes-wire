#!/usr/bin/env bash
set -euo pipefail


usage(){ echo "Usage: $0 --grm_prefix <pfx> --pheno gcta.pheno --pairs pairs.tsv --threads N --out_prefix out__"; exit 1; }
GRM= PHENO= PAIRS= THREADS=1 OUT=he_bivar__

while [[ $# -gt 0 ]]; do
  case "$1" in
    --grm_prefix) GRM="$2"; shift 2;;
    --pheno)      PHENO="$2"; shift 2;;
    --pairs)      PAIRS="$2"; shift 2;;
    --threads)    THREADS="$2"; shift 2;;
    --out_prefix) OUT="$2"; shift 2;;
    --gcta_exec) GCTA_BIN="$2"; shift 2;;
    *) echo "Unknown arg $1"; usage;;
  esac
done

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

[[ -f "${PHENO}" && -f "${PAIRS}" ]] || usage
while IFS=$'\t' read -r t1 t2 col1 col2; do
  tag="${t1}__VS__${t2}"
  tag="${tag//[^A-Za-z0-9_]/_}"
  RUN_GCTA "${GRM_ARG[@]}" \
         --pheno "${PHENO}" \
         --HEreg-bivar "${col1}" "${col2}" \
         --thread-num "${THREADS}" \
         --out "${OUT}${tag}" || true
done < <(tail -n +2 "${PAIRS}")
