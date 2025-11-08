#!/usr/bin/env bash
# build_mgrm_with_srm.sh
#
# Build all-SNP GRM (text.gz) from PLINK and harmonize it with DAM/CAGE SRMs
# so all three share the *same* samples (names+order), using the cohort from
# MAKE_GENESIS_COHORT_IDS (scan IDs). Produces an --mgrm-gz list.
#
# Assumptions
#   - scan_lookup has columns akin to: scan_id, sample_label (names may vary)
#   - SRM .grm.id have label-like IIDs; keep substring before first "_"
#   - PLINK IIDs are scan IDs (candidates_kept.* already restricted to cohort)
#
# Outputs (CWD)
#   grm_all_txt.grm.gz, grm_all_txt.grm.id
#   grm_all_harmo.grm.gz, grm_all_harmo.grm.id
#   srm_dam_harmo.grm.gz, srm_dam_harmo.grm.id
#   srm_cage_harmo.grm.gz, srm_cage_harmo.grm.id
#   mgrm_with_srm.txt
set -euo pipefail

# ----------------------- CLI -----------------------
GCTA_EXEC="gcta64"
THREADS=4
BFILE=""
DAM_GZ="" ; DAM_ID=""
CAGE_GZ=""; CAGE_ID=""
SCAN_LOOKUP=""
SCAN_IDS=""
MGRM_OUT="mgrm_with_srm.txt"

usage() {
  cat <<EOF
Usage: $0 --bfile <prefix> --dam_gz <file> --dam_id <file> --cage_gz <file> --cage_id <file> \\
         --scan_lookup <tsv/csv> --scan_ids <txt> [--gcta_exec PATH] [--threads N] [--mgrm_out FILE]
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gcta_exec)   GCTA_EXEC="$2"; shift 2;;
    --threads)     THREADS="$2"; shift 2;;
    --bfile)       BFILE="$2"; shift 2;;
    --dam_gz)      DAM_GZ="$2"; shift 2;;
    --dam_id)      DAM_ID="$2"; shift 2;;
    --cage_gz)     CAGE_GZ="$2"; shift 2;;
    --cage_id)     CAGE_ID="$2"; shift 2;;
    --scan_lookup) SCAN_LOOKUP="$2"; shift 2;;
    --scan_ids)    SCAN_IDS="$2"; shift 2;;
    --mgrm_out)    MGRM_OUT="$2"; shift 2;;
    -h|--help)     usage;;
    *) echo "[mgrm] Unknown arg: $1"; usage;;
  esac
done

[[ -n "$BFILE" && -s "${BFILE}.bed" && -s "${BFILE}.bim" && -s "${BFILE}.fam" ]] \
  || { echo "[mgrm] ERROR: missing PLINK files for --bfile '$BFILE'"; exit 2; }
[[ -s "$DAM_GZ"  && -s "$DAM_ID"  ]] || { echo "[mgrm] ERROR: missing DAM SRM files"; exit 2; }
[[ -s "$CAGE_GZ" && -s "$CAGE_ID" ]] || { echo "[mgrm] ERROR: missing CAGE SRM files"; exit 2; }
[[ -s "$SCAN_LOOKUP" ]] || { echo "[mgrm] ERROR: missing --scan_lookup"; exit 2; }
[[ -s "$SCAN_IDS"    ]] || { echo "[mgrm] ERROR: missing --scan_ids"; exit 2; }

# ----------------------- Resolve GCTA runner -----------------------
GCTA_PATH="$(command -v "$GCTA_EXEC" || true)"
[[ -n "$GCTA_PATH" ]] || { echo "[mgrm] ERROR: cannot find '$GCTA_EXEC' in PATH"; exit 127; }
if "$GCTA_PATH" --help >/dev/null 2>&1; then
  RUN_GCTA() { "$GCTA_PATH" "$@"; }
elif "$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
  RUN_GCTA() { "$GCTA_PATH" --appimage-extract-and-run "$@"; }
else
  echo "[mgrm] ERROR: gcta64 not runnable"; exit 127
fi

# ----------------------- GRM (text.gz) from PLINK -------------------------
echo "[mgrm] Building all-SNP GRM (text.gz) from ${BFILE}.*"
RUN_GCTA --bfile "$BFILE" --autosome --make-grm-gz --thread-num "$THREADS" --out grm_all_txt

# ----------------------- Detect delimiter in scan_lookup -------------------
head -n1 "$SCAN_LOOKUP" | tr -d '\r' > .scan_header.txt
if grep -q $'\t' .scan_header.txt; then
  FS_RE=$'\t'; DELIM_NAME="TAB"
elif grep -q ',' .scan_header.txt; then
  FS_RE=',';   DELIM_NAME="COMMA"
elif grep -q ';' .scan_header.txt; then
  FS_RE=';';   DELIM_NAME="SEMICOLON"
else
  FS_RE='[[:space:]]+'; DELIM_NAME="WHITESPACE"
fi
echo "[mgrm] scan_lookup delimiter detected: ${DELIM_NAME}"

# ----------------------- label -> scan_id map ------------------------------
# Normalize header names (lowercase, non-alnum -> '_'); accept common aliases.
awk -v FS="$FS_RE" '
  BEGIN{ OFS="\t" }
  NR==1{
    for(i=1;i<=NF;i++){
      h=$i; gsub(/\r/,"",h); sub(/^\xef\xbb\xbf/,"",h)
      l=tolower(h); gsub(/[^a-z0-9]+/,"_",l)
      if(!scan && (l=="scan_id" || l=="scan" || l=="iid")) scan=i
      if(!lab  && (l=="sample_label" || l=="label" || l=="sample" || l=="sample_id" || l=="rfid" || l=="id")) lab=i
    }
    if(!scan || !lab){
      printf("[mgrm] ERROR: could not detect scan/label columns in scan_lookup header (") >"/dev/stderr"
      for(i=1;i<=NF;i++){ l=tolower($i); gsub(/[^a-z0-9]+/,"_",l); printf(" %s", l) >"/dev/stderr" }
      print " )" >"/dev/stderr"; exit 2
    }
    printf("[mgrm] header indices: scan=%d label=%d\n", scan, lab) >"/dev/stderr"
    next
  }
  { print $lab, $scan }
' "$SCAN_LOOKUP" | awk -F'\t' 'BEGIN{OFS="\t"} !seen[$1]++' > .label_to_scan.tsv

nmap=$(wc -l < .label_to_scan.tsv || echo 0)
echo "[mgrm] label->scan map rows: $nmap"

# ----------------------- Cohort order (scan_ids) ---------------------------
awk '!seen[$1]++{print $1}' "$SCAN_IDS" > .target.scan
ntgt=$(wc -l < .target.scan || echo 0)
[[ "$ntgt" -gt 1 ]] || { echo "[mgrm] ERROR: cohort has <2 samples"; exit 3; }
echo "[mgrm] cohort scan_ids: $ntgt"

# ----------------------- Index bases ---------------------------------------
awk '{print NR, $2}' grm_all_txt.grm.id > .all.idx2scan
awk '{print $2}'      grm_all_txt.grm.id > .all.scan

awk '{x=$2; sub(/_.*/,"",x); print NR, x}' "$DAM_ID"  > .dam.idx2labelbase
awk '{x=$2; sub(/_.*/,"",x); print NR, x}' "$CAGE_ID" > .cage.idx2labelbase

awk 'NR==FNR{m[$1]=$2; next} { if($2 in m) print $1, m[$2] }' .label_to_scan.tsv .dam.idx2labelbase  > .dam.idx2scan
awk 'NR==FNR{m[$1]=$2; next} { if($2 in m) print $1, m[$2] }' .label_to_scan.tsv .cage.idx2labelbase > .cage.idx2scan

awk '{print $2}' .dam.idx2scan  | sort -u > .dam.scan
awk '{print $2}' .cage.idx2scan | sort -u > .cage.scan

# ----------------------- Intersection & final order ------------------------
# Order is cohort order; require presence in all three (allSNP, dam, cage).
awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $1}' .dam.scan .cage.scan | sort -u > .srm_both.scan
awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $1}' .srm_both.scan .all.scan > .in_all_three.scan
awk 'NR==FNR{keep[$1]=1; next} ($1 in keep){print $1}' .in_all_three.scan .target.scan > target.scan

ntarget=$(wc -l < target.scan || echo 0)
echo "[mgrm] final target n=$ntarget"
[[ "$ntarget" -gt 1 ]] || { echo "[mgrm] ERROR: <2 samples after intersection (cohort ∩ allSNP ∩ dam ∩ cage)"; exit 3; }

# ----------------------- Index maps old -> new -----------------------------
awk 'NR==FNR{ord[$1]=++k; next} { if($2 in ord) print $1, ord[$2] }' target.scan .all.idx2scan > .all.idx2new
awk 'NR==FNR{ord[$1]=++k; next} { if($2 in ord) print $1, ord[$2] }' target.scan .dam.idx2scan  > .dam.idx2new
awk 'NR==FNR{ord[$1]=++k; next} { if($2 in ord) print $1, ord[$2] }' target.scan .cage.idx2scan > .cage.idx2new

# ----------------------- Remapper (lower triangle text.gz) ----------------
remap_grm () {
  local src_gz="$1"; local map_file="$2"; local out_prefix="$3"
  gzip -dc "$src_gz" | awk '
    NR==FNR { new[$1]=$2; next }
    { i=$1; j=$2; if(!(i in new) || !(j in new)) next;
      ii=new[i]; jj=new[j]; if(ii<jj){t=ii; ii=jj; jj=t}
      printf("%d %d", ii, jj);
      for(k=3;k<=NF;k++) printf(" %s", $k);
      printf("\n");
    }' "$map_file" - | gzip -c > "${out_prefix}.grm.gz"
}

# ----------------------- Write harmonized IDs ------------------------------
# Keep grm_all_harmo.grm.id consistent with the original PLINK-based FID/IID (FID very often 0 from your pipeline)
awk 'NR==FNR{ord[$1]=++k; next} { if($2 in ord) print ord[$2] "\t" $0 }' target.scan grm_all_txt.grm.id \
  | sort -n | cut -f2- > grm_all_harmo.grm.id

# IMPORTANT: Force SRM IDs to the same two-column (FID,IID) convention as grm_all + pheno: FID=0, IID=<scan_id>
awk '{print 0, $1}' target.scan > srm_dam_harmo.grm.id
cp srm_dam_harmo.grm.id srm_cage_harmo.grm.id


# ----------------------- Write harmonized GRMs -----------------------------
remap_grm "grm_all_txt.grm.gz" "$PWD/.all.idx2new"  "grm_all_harmo"
remap_grm "$DAM_GZ"            "$PWD/.dam.idx2new"  "srm_dam_harmo"
remap_grm "$CAGE_GZ"           "$PWD/.cage.idx2new" "srm_cage_harmo"

# Quick sanity
expected_lines=$(( ntarget * (ntarget + 1) / 2 ))
for p in grm_all_harmo srm_dam_harmo srm_cage_harmo; do
  n=$(gzip -dc ${p}.grm.gz | wc -l || echo 0)
  echo "[mgrm] ${p}.grm.gz lines: $n (expected $expected_lines)"
done

# ----------------------- MGRM list ----------------------------------------
printf "%s\n%s\n%s\n" \
  "$(pwd)/grm_all_harmo" \
  "$(pwd)/srm_dam_harmo" \
  "$(pwd)/srm_cage_harmo" > "$MGRM_OUT"

echo "[mgrm] wrote $MGRM_OUT:"
nl -ba "$MGRM_OUT"
