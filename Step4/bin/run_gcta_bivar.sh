#!/usr/bin/env bash
set -euo pipefail

usage(){ cat <<'USAGE'
Usage:
  run_gcta_bivar.sh \
    --grm_prefix <prefix|mgrm.txt|mgrm_gz.txt> \
    --pheno <ORIGINAL phenotype with header> \
    [--pairs ALL | --pairs pairs.tsv] \
    [--expand-from pairs.tsv] \
    [--threads N] \
    [--out_prefix out__] \
    [--gcta_exec gcta64] \
    [--reml_maxit 100] \
    [--allow_neg_var] \
    [--allow_rg_outside] \
    [--nocove] \
    [--lrt_rg <value>]

Notes
- Point --pheno to the ORIGINAL file with a header; the script creates pheno.harmo internally.
- If pheno.harmo loses/garbles the header, we fall back to the original header for naming.
- --expand-from pairs.tsv: take the union of trait1/trait2, build all unique unordered pairs among them, and run those.
- Defaults: constrained REML, residual covariance modeled.
  * --allow_neg_var        -> GCTA --reml-no-constrain (variance components may be negative)
  * --allow_rg_outside     -> GCTA --reml-bivar-no-constrain (r_g not forced to [-1,1])
  * --nocove               -> GCTA --reml-bivar-nocove (drop residual covariance term)
  * --lrt_rg <v>           -> GCTA --reml-bivar-lrt-rg <v> (LRT of r_g fixed at v, e.g. 0)
USAGE
exit 1; }

# -------------------- args --------------------
GRM="" PHENO_IN="" PAIRS="ALL" EXPAND_FROM="" THREADS=1 OUT="rg_bivar_all__" GCTA_BIN="gcta64" REML_MAXIT=100
ALLOW_NEG_VAR=0 ALLOW_RG_OUTSIDE=0 NOCOVE=0 LRT_RG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --grm_prefix) GRM="$2"; shift 2;;
    --pheno)      PHENO_IN="$2"; shift 2;;
    --pairs)      PAIRS="$2"; shift 2;;
    --expand-from) EXPAND_FROM="$2"; shift 2;;
    --threads)    THREADS="$2"; shift 2;;
    --out_prefix) OUT="$2"; shift 2;;
    --gcta_exec)  GCTA_BIN="$2"; shift 2;;
    --reml_maxit) REML_MAXIT="$2"; shift 2;;
    --allow_neg_var)    ALLOW_NEG_VAR=1; shift 1;;
    --allow_rg_outside) ALLOW_RG_OUTSIDE=1; shift 1;;
    --nocove)            NOCOVE=1; shift 1;;
    --lrt_rg)           LRT_RG="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg $1"; usage;;
  esac
done

[[ -n "$GRM" && -n "$PHENO_IN" ]] || usage
[[ -r "$PHENO_IN" ]] || { echo "[bivar] ERROR: cannot read --pheno '$PHENO_IN'"; exit 2; }
[[ -z "${EXPAND_FROM}" || -r "$EXPAND_FROM" ]] || { echo "[bivar] ERROR: cannot read --expand-from '$EXPAND_FROM'"; exit 2; }

# -------------------- GRM handling --------------------
# Accept: single prefix (--grm*), list of prefixes (--mgrm or --mgrm-bin), or list of gz GRMs (--mgrm-gz)
if [[ -f "$GRM" && "${GRM##*.}" == "txt" ]]; then
  # Heuristic: if the first listed GRM ends with .grm.gz on disk, use --mgrm-gz; otherwise assume binary (--mgrm)
  first_line="$(head -n1 "$GRM" | tr -d '\r')"
  if [[ -f "${first_line}.grm.gz" ]]; then
    GRM_ARG=(--mgrm-gz "$GRM")
  else
    GRM_ARG=(--mgrm "$GRM")
  fi
  FIRST_PREFIX="$first_line"
  MGRM_LIST="$GRM"
else
  # Single GRM prefix (binary or gz text are both handled by GCTA via --grm / --grm-gz; prefer binary)
  if [[ -f "${GRM}.grm.bin" ]]; then
    GRM_ARG=(--grm "$GRM")
  elif [[ -f "${GRM}.grm.gz" ]]; then
    GRM_ARG=(--grm-gz "$GRM")
  else
    echo "[bivar] ERROR: cannot infer GRM files from prefix '$GRM'"; exit 2
  fi
  FIRST_PREFIX="$GRM"
  MGRM_LIST=""
fi
FIRST_ID="${FIRST_PREFIX}.grm.id"
[[ -s "$FIRST_ID" ]] || { echo "[bivar] ERROR: missing $FIRST_ID"; exit 2; }

# -------------------- Harmonize phenotype --------------------
harmonize_pheno_to_grm.sh --pheno "$PHENO_IN" --grm_id "$FIRST_ID" --out pheno.harmo
PHENO="pheno.harmo"

# -------------------- Helpers --------------------
get_traits_from_header(){ awk 'NR==1{for(i=3;i<=NF;i++){gsub(/\r/,"",$i); print $i}}' "$1"; }
sanitize(){ local s="$1"; s="${s//$'\r'/}"; s="${s// /_}"; s="${s//[^A-Za-z0-9_]/_}"; echo "$s"; }
normalize(){ local s="$1"; s="${s//$'\r'/}"; echo "$s"; }

# -------------------- Trait names (header rescue) --------------------
mapfile -t TRAITS_H < <(get_traits_from_header "$PHENO" || true)
if (( ${#TRAITS_H[@]} == 0 )) || awk 'NR==1{n=0; for(i=3;i<=NF;i++) if($i ~ /[A-Za-z_]/) n++; print (n?"OK":"NUM")}' "$PHENO" | grep -q NUM; then
  mapfile -t TRAITS_H < <(get_traits_from_header "$PHENO_IN")
fi
[[ ${#TRAITS_H[@]} -gt 1 ]] || { echo "[bivar] ERROR: could not infer trait names from headers"; exit 5; }

NT=$(awk 'NR>1{print NF-2; exit}' "$PHENO")
(( NT >= 2 )) || { echo "[bivar] ERROR: need >=2 traits in '$PHENO'"; exit 6; }
if (( NT != ${#TRAITS_H[@]} )); then
  echo "[bivar] ERROR: header name count (${#TRAITS_H[@]}) != trait columns in harmonized file ($NT)"; exit 7
fi

declare -A IDX
for ((i=1;i<=NT;i++)); do
  nm="$(normalize "${TRAITS_H[i-1]}")"
  [[ -z "$nm" ]] && { echo "[bivar] ERROR: empty trait name at position $i in header."; exit 12; }
  [[ -n "${IDX[$nm]:-}" ]] && { echo "[bivar] ERROR: duplicate trait name in header: '$nm'."; exit 13; }
  IDX["$nm"]="$i"
done

# -------------------- Build PAIRS --------------------
PAIRS="ALL"
pairs_file=""
if [[ "$PAIRS" == "ALL" && -z "${EXPAND_FROM}" ]]; then
  pairs_file="$(mktemp)"; printf "t1\tt2\tcol1\tcol2\n" > "$pairs_file"
  for ((i=1;i<=NT;i++)); do
    for ((j=i+1;j<=NT;j++)); do
      printf "%s\t%s\t%d\t%d\n" "${TRAITS_H[i-1]}" "${TRAITS_H[j-1]}" "$i" "$j" >> "$pairs_file"
    done
  done
elif [[ -n "${EXPAND_FROM}" ]]; then
  pairs_file="$(mktemp)"; printf "t1\tt2\tcol1\tcol2\n" > "$pairs_file"
  declare -A WANT
  { read -r _hdr || true
    while IFS=$'\t' read -r a b _rest; do
      [[ -z "${a:-}" && -z "${b:-}" ]] && continue
      a="$(normalize "${a:-}")"; b="$(normalize "${b:-}")"
      [[ -n "$a" ]] && WANT["$a"]=1
      [[ -n "$b" ]] && WANT["$b"]=1
    done
  } < "$EXPAND_FROM"
  mapfile -t WANT_NAMES < <(printf "%s\n" "${!WANT[@]}" | LC_ALL=C sort)
  (( ${#WANT_NAMES[@]} > 1 )) || { echo "[bivar] ERROR: need >=2 distinct trait names in $EXPAND_FROM"; exit 14; }
  missing=0
  for nm in "${WANT_NAMES[@]}"; do
    if [[ -z "${IDX[$nm]:-}" ]]; then echo "[bivar] ERROR: trait '$nm' from $EXPAND_FROM not found in phenotype header."; missing=1; fi
  done
  (( missing == 0 )) || exit 8
  for ((a=0;a<${#WANT_NAMES[@]};a++)); do
    for ((b=a+1;b<${#WANT_NAMES[@]};b++)); do
      nm1="${WANT_NAMES[a]}"; nm2="${WANT_NAMES[b]}"
      printf "%s\t%s\t%d\t%d\n" "$nm1" "$nm2" "${IDX[$nm1]}" "${IDX[$nm2]}" >> "$pairs_file"
    done
  done
else
  [[ -r "$PAIRS" ]] || { echo "[bivar] ERROR: pairs file not found: $PAIRS"; exit 9; }
  pairs_file="$(mktemp)"; printf "t1\tt2\tcol1\tcol2\n" > "$pairs_file"
  declare -A SEEN
  { read -r _hdr || true
    while IFS=$'\t' read -r t1 t2 c1 c2; do
      t1="$(normalize "${t1:-}")"; t2="$(normalize "${t2:-}")"
      [[ -z "$t1" || -z "$t2" || -z "${c1:-}" || -z "${c2:-}" ]] && continue
      [[ "$t1" == "$t2" ]] && continue
      if (( c1 < c2 )); then key="${c1},${c2}"; o1="$t1"; o2="$t2"; k1="$c1"; k2="$c2";
      else                   key="${c2},${c1}"; o1="$t2"; o2="$t1"; k1="$c2"; k2="$c1";
      fi
      [[ -n "${SEEN[$key]:-}" ]] && continue
      SEEN[$key]=1
      printf "%s\t%s\t%d\t%d\n" "$o1" "$o2" "$k1" "$k2" >> "$pairs_file"
    done
  } < "$PAIRS"
fi
[[ -s "$pairs_file" ]] || { echo "[bivar] ERROR: no pairs to run"; exit 10; }

# -------------------- gcta launcher --------------------
GCTA_PATH="$(command -v "$GCTA_BIN" || true)"
[[ -n "$GCTA_PATH" ]] || { echo "[mgrm] ERROR: cannot find $GCTA_BIN"; exit 127; }
if "$GCTA_PATH" --help >/dev/null 2>&1; then
  RUN_GCTA() { "$GCTA_PATH" "$@"; }
elif "$GCTA_PATH" --appimage-help >/dev/null 2>&1; then
  RUN_GCTA() { "$GCTA_PATH" --appimage-extract-and-run "$@"; }
else
  echo "[mgrm] ERROR: gcta not runnable"; exit 127
fi

# -------------------- run all pairs --------------------
tail -n +2 "$pairs_file" | while IFS=$'\t' read -r t1 t2 col1 col2; do
  tag="$(sanitize "$t1")__VS__$(sanitize "$t2")"
  echo "[bivar] Running ${t1} (col ${col1}) vs ${t2} (col ${col2}) → ${OUT}${tag}.hsq"
  args=( "${GRM_ARG[@]}"
         --pheno "$PHENO"
         --reml-bivar "$col1" "$col2"
         --reml-maxit "$REML_MAXIT"
         --thread-num "$THREADS"
         --out "${OUT}${tag}" )
  (( NOCOVE ))            && args+=( --reml-bivar-nocove )
  (( ALLOW_NEG_VAR ))     && args+=( --reml-no-constrain )
  (( ALLOW_RG_OUTSIDE ))  && args+=( --reml-bivar-no-constrain )
  [[ -n "$LRT_RG" ]]      && args+=( --reml-bivar-lrt-rg "$LRT_RG" )

  printf "[bivar] GCTA opts: %q " "${args[@]}"; echo
  RUN_GCTA "${args[@]}" || echo "[warn] Non-zero exit for ${t1} vs ${t2} (continuing)"
done

# -------------------- sanity: ID intersection --------------------
cut -f1-2 "$PHENO" | awk '{print $1" "$2}' > .phids
if [[ -n "$MGRM_LIST" ]]; then
  mapfile -t PREFS < <(cat "$MGRM_LIST")
  cut -d' ' -f1,2 "${PREFS[0]}.grm.id" > .ids_common
  for ((k=1;k<${#PREFS[@]};k++)); do
    cut -d' ' -f1,2 "${PREFS[k]}.grm.id" | sort -u > .ids_k
    comm -12 <(sort -u .ids_common) .ids_k > .ids_new && mv .ids_new .ids_common
  done
else
  cut -d' ' -f1,2 "${FIRST_PREFIX}.grm.id" > .ids_common
fi
if ! comm -12 <(sort -u .ids_common) <(sort -u .phids) | head -n1 | grep -q .; then
  echo "[check] ERROR: empty intersection across GRM IDs and phenotype IDs"; exit 11
fi
