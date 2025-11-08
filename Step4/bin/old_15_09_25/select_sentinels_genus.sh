#!/usr/bin/env bash
set -euo pipefail

ASSOC_DIR="${1}"
BIM="${2}"
GENE_LOC="${3}"
GENUS_PREFIX="${4}"
WUP="${5:-0}"
WDN="${6:-0}"

log(){ echo "[sentinels_genus] $*"; }

OUT="sentinels_genus.tsv"
# NOTE extra third column: variant_id
printf "gene\tgenus\tvariant_id\tSNP\tP\tCHR\tBP\n" > "${OUT}"

log "Indexing BIM → (KEY=CHR:BP → variant_id,SNP,CHR,BP)"
# bmap.tsv columns: 1=variant_id (line #), 2=KEY, 3=SNP, 4=CHR, 5=BP
awk 'BEGIN{OFS="\t"}{print NR, $1":"$4, $2, $1, $4}' "${BIM}" \
  | sort -k2,2 > bmap.tsv
log "BIM rows: $(wc -l < "${BIM}") | map keys: $(wc -l < bmap.tsv)"

log "Loading genes from: ${GENE_LOC}"
awk -v U="${WUP}" -v D="${WDN}" 'BEGIN{OFS="\t"}
{
  gene=$1; chr=$2; gsub(/^chr/,"",chr);
  start=$3+0; end=$4+0; if (start>end){t=start; start=end; end=t}
  wstart=start-U; if (wstart<1) wstart=1; wend=end+D;
  print gene, chr, start, end, wstart, wend
}' "${GENE_LOC}" > genes.tsv
log "Genes loaded: $(wc -l < genes.tsv) (from ${GENE_LOC})"

ls -1 "${ASSOC_DIR}/by_trait/${GENUS_PREFIX}"*.tsv.gz 2>/dev/null \
  | sed -E "s#.*/${GENUS_PREFIX}([^.]*)\\.tsv.*#\\1#" \
  | sort -u > genera.list || true
ngen=$(wc -l < genera.list 2>/dev/null || echo 0)
log "Genera with assoc files: ${ngen}"
[[ "${ngen}" -gt 0 ]] || { log "No assoc files found; wrote header only"; exit 0; }

# Assoc columns: (1)variant.id (2)chr (3)pos (4)Score.pval ...
while IFS=$'\t' read -r GENE CHR START END WSTART WEND; do
  while read -r GENUS; do
    f="${ASSOC_DIR}/by_trait/${GENUS_PREFIX}${GENUS}.tsv.gz"
    [[ -s "$f" ]] || continue

    tmp_in=$(mktemp)
    gzip -cd "$f" \
      | awk -F'\t' -v c="${CHR}" -v s="${WSTART}" -v e="${WEND}" '
          NR==1 { next }  # skip header
          $2==c && $3>=s && $3<=e && $4==$4 { print $2 ":" $3 "\t" $4 }' \
      > "${tmp_in}" || true

    if [[ -s "${tmp_in}" ]]; then
      # Read BIM map (by KEY=CHR:BP), keep the best P in-window, and print with variant_id
      awk -F $'\t' -v G="${GENE}" -v GU="${GENUS}" '
        NR==FNR { m[$2] = $1 "\t" $3 "\t" $4 "\t" $5; next }  # key → "variant_id\tSNP\tCHR\tBP"
        {
          k = $1; p = $2 + 0
          if (k in m) {
            if (!has || p < minp) { has=1; minp = p; best = m[k] }
          }
        }
        END {
          if (has) {
            n = split(best, a, "\t")  # a[1]=variant_id, a[2]=SNP, a[3]=CHR, a[4]=BP
            printf "%s\t%s\t%s\t%s\t%.16g\t%s\t%s\n", G, GU, a[1], a[2], minp, a[3], a[4]
          }
        }
      ' bmap.tsv "${tmp_in}" >> "${OUT}"
    fi
    rm -f "${tmp_in}" || true
  done < genera.list
done < genes.tsv

log "Wrote: ${OUT}  rows: $(($(wc -l < "${OUT}")-1))"
