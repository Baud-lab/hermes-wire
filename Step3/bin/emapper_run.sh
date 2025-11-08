#!/usr/bin/env bash
set -euo pipefail

# Simple, verbose eggNOG-mapper runner
# Makes <basename>.annotations.tsv (cleaned) in CWD

usage() {
  cat <<EOF
Usage: emapper_run.sh --fasta FILE --mode {diamond|mmseqs} --threads N [--extra '...']
Outputs: <basename>.annotations.tsv in the current directory
EOF
}

FASTA=""
MODE="diamond"
THREADS=4
EXTRA=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --extra) EXTRA="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[emapper_run] Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -f "$FASTA" ]] || { echo "[emapper_run] FASTA not found: $FASTA" >&2; exit 2; }
BASENAME="$(basename "$FASTA")"
STEM="${BASENAME%.*}"

echo "[emapper_run] FASTA=$FASTA"
echo "[emapper_run] MODE=$MODE  THREADS=$THREADS"
echo "[emapper_run] EXTRA=$EXTRA"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

emapper.py \
  -i "$FASTA" \
  -m "$MODE" \
  --cpu "$THREADS" \
  --override \
  --output "$STEM" \
  --output_dir "$TMPDIR" \
  $EXTRA

# Clean header (drop '##', remove leading '#')
ANN="$TMPDIR/${STEM}.emapper.annotations"
[[ -f "$ANN" ]] || { echo "[emapper_run] Missing emapper output: $ANN" >&2; exit 3; }

grep -v '^##' "$ANN" | sed 's/^#//' > "${STEM}.annotations.tsv"
echo "[emapper_run] Wrote ${STEM}.annotations.tsv"
