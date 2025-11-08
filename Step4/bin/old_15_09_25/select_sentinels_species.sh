#!/usr/bin/env bash
# Pick one sentinel per (gene, genus) by scanning all species files of that genus
# Inputs:
#   $1 = ASSOC_DIR         (assoc_species_full)
#   $2 = BIM               (candidates_kept.bim)
#   $3 = GENE_LOC          (magma.genes.loc; GENE_ID,CHR,START,END; no header)
#   $4 = SPECIES_PREFIX    (e.g., s__)
#   $5 = WINDOW_UP_BP
#   $6 = WINDOW_DOWN_BP
set -euo pipefail

ASSOC_DIR="${1}"
BIM="${2}"
GENE_LOC="${3}"
SPECIES_PREFIX="${4}"
WUP="${5:-0}"
WDN="${6:-0}"

log(){ echo "[sentinels_species] $*"; }

OUT="sentinels_species.tsv"
# NOTE now includes variant_id to mirror genus version
printf "gene\tgenus\tvariant_id\tSNP\tP\tCHR\tBP\tspecies_minp\n" > "${OUT}"

log "Indexing BIM → (KEY=CHR:BP → variant_id,SNP,CHR,BP)"
# bmap.tsv: 1=variant_id (line #), 2=KEY, 3=SNP, 4=CHR, 5=BP
awk 'BEGIN{OFS="\t"}{print NR, $1":"$4, $2, $1, $4}' "${BIM}" \
  | sort -k2,2 > bmap.tsv
log "BIM rows: $(wc -l < "${BIM}") | map keys: $(wc -l < bmap.tsv)"

log "Loading genes from: ${GENE_LOC}"
# Normalize gene intervals; strip leading 'chr' on CHR; apply optional window
awk -v U="${WUP}" -v D="${WDN}" 'BEGIN{OFS="\t"}
{
  gene=$1; chr=$2; gsub(/^chr/,"",chr);
  start=$3+0; end=$4+0; if (start>end){t=start; start=end; end=t}
  wstart=start-U; if (wstart<1) wstart=1; wend=end+D;
  print gene, chr, start, end, wstart, wend
}' "${GENE_LOC}" > genes.tsv
log "Genes loaded: $(wc -l < genes.tsv)"

# Determine genera present at species level
ls -1 "${ASSOC_DIR}/by_trait/${SPECIES_PREFIX}"*.tsv.gz 2>/dev/null \
  | sed -E "s#.*/${SPECIES_PREFIX}([^_]+)_.*#\\1#" \
  | sort -u > genera.list || true
ngen=$(wc -l < genera.list 2>/dev/null || echo 0)
log "Genera with species assoc files: ${ngen}"
[[ "${ngen}" -gt 0 ]] || { log "No species assoc files found; wrote header only"; exit 0; }

# For each gene × genus: scan all species files of that genus and take min p
# Species assoc columns (after prefixing species in col 1):
#   $1=species $2=variant.id $3=chr $4=pos $5=Score.pval ...
while IFS=$'\t' read -r GENE CHR START END WSTART WEND; do
  while read -r GENUS; do
    shopt -s nullglob
    # Stream rows from all species files for this genus, prefixing species in $1
    if compgen -G "${ASSOC_DIR}/by_trait/${SPECIES_PREFIX}${GENUS}"_*.tsv.gz > /dev/null; then
      awk -F $'\t' -v G="${GENE}" -v GU="${GENUS}" -v C="${CHR}" -v S="${WSTART}" -v E="${WEND}" '
        BEGIN { OFS="\t" }
        NR == FNR { m[$2] = $1 "\t" $3 "\t" $4 "\t" $5; next }  # key → "variant_id\tSNP\tCHR\tBP"
        {
          # $1=species, $2=variant.id, $3=chr, $4=pos, $5=Score.pval (header already stripped)
          if ($3==C && $4>=S && $4<=E && $5==$5) {
            key = $3 ":" $4
            p = $5 + 0
            if (key in m) {
              if (!has || p < minp) {
                has   = 1
                minp  = p
                best  = m[key]
                bestsp= $1
              }
            }
          }
        }
        END {
          if (has) {
            n = split(best, a, "\t")  # a[1]=variant_id, a[2]=SNP, a[3]=CHR, a[4]=BP
            printf "%s\t%s\t%s\t%s\t%.16g\t%s\t%s\t%s\n", G, GU, a[1], a[2], minp, a[3], a[4], bestsp
          }
        }
      ' bmap.tsv \
        <( for f in "${ASSOC_DIR}/by_trait/${SPECIES_PREFIX}${GENUS}"_*.tsv.gz; do
             sp="$(basename "$f" .tsv.gz)"; sp="${sp%_assoc}"
             gzip -cd "$f" | awk -v SP="$sp" 'BEGIN{FS=OFS="\t"} NR>1 { print SP, $0 }'
           done ) >> "${OUT}" || true
    fi
  done < genera.list
done < genes.tsv

log "Wrote: ${OUT}  rows: $(($(wc -l < "${OUT}")-1))"
