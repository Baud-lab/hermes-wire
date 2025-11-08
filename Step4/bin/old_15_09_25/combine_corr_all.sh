#!/usr/bin/env bash
set -euo pipefail
out="corr_by_gene.tsv"
shopt -s nullglob
tsps=( *.tsv )
if [[ ${#tsps[@]} -eq 0 ]]; then
  printf "genus\tgene\tn_species_used\tspearman_rho\tspearman_p\tq_spearman\tpgls_slope\tpgls_p\tq_pgls\n" > "${out}"
  exit 0
fi
head -n1 "${tsps[0]}" > "${out}"
for f in "${tsps[@]}"; do
  tail -n +2 "${f}" >> "${out}"
done
