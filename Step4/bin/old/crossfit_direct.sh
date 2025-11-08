#!/usr/bin/env bash
set -euo pipefail

# crossfit_direct.sh – builds a permissive (genus × gene × species) table and
# runs extract_crossfit_sentinels.py. No Excel dependency; no MAGMA proper.

assoc_dir=""
gene_loc=""
snp_loc=""
alpha="0.10"
outdir="."
genera_csv=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assoc_dir)  assoc_dir="$2";  shift 2;;
    --gene_loc)   gene_loc="$2";   shift 2;;
    --snp_loc)    snp_loc="$2";    shift 2;;
    --alpha)      alpha="$2";      shift 2;;
    --outdir)     outdir="$2";     shift 2;;
    --genera)     genera_csv="$2"; shift 2;;
    *) echo "[crossfit_direct] Unknown arg: $1" >&2; exit 2;;
  esac
done

log() { echo "[crossfit_direct] $*"; }
ts()  { date +"%F %T"; }

[[ -d "$assoc_dir" ]] || { echo "assoc_dir not found: $assoc_dir" >&2; exit 2; }
[[ -f "$gene_loc"  ]] || { echo "gene_loc not found: $gene_loc"   >&2; exit 2; }
[[ -f "$snp_loc"   ]] || { echo "snp_loc not found:  $snp_loc"    >&2; exit 2; }

mkdir -p "$outdir"
bundle="$outdir/magma_bundle"; mkdir -p "$bundle"
cp -f "$gene_loc" "$bundle/magma.genes.loc"
cp -f "$snp_loc"  "$bundle/magma.snps.loc"

# Inventory + permissive table
python3 - "$assoc_dir" "$genera_csv" <<'PY'
import os, sys, re, glob, pandas as pd
assoc_dir = sys.argv[1]
glist = set(sys.argv[2].split(",")) if len(sys.argv)>2 and sys.argv[2] else set()
paths = sorted(glob.glob(os.path.join(assoc_dir,"**","*.tsv.gz"), recursive=True))
print(f"[crossfit_direct] Assoc files (*.tsv.gz): {len(paths)}", flush=True)
traits=[]
for p in paths:
    t = re.sub(r"\.tsv\.gz$","", os.path.basename(p))
    if t.startswith("s__"):
        gen = t[3:].split("_",1)[0] if "_" in t[3:] else t[3:]
        if (not glist) or (gen in glist):
            traits.append((gen, t))
traits = sorted(set(traits))
print(f"[crossfit_direct] Species traits after genus filter: {len(traits)}", flush=True)
gl = pd.read_csv("magma_bundle/magma.genes.loc", sep=r"\s+", header=None, engine="python").iloc[:, :1]
gl.columns = ["GENE_ID"]; gl["GENE_ID"]=gl["GENE_ID"].astype(str)
if len(gl) and traits:
    df = pd.DataFrame(traits, columns=["genus","species"])
    out = df.merge(gl.rename(columns={"GENE_ID":"gene"}), how="cross")
    out[["genus","gene","species"]].to_csv("magma_bundle/magma_sig_species.tsv", sep="\t", index=False)
    print(f"[crossfit_direct] Rows in magma_sig_species.tsv: {len(out)}", flush=True)
else:
    open("magma_bundle/magma_sig_species.tsv","w").write("genus\tgene\tspecies\n")
PY

log "Selecting sentinels (variant.id-aware)"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PY_HELPER="${SCRIPT_DIR}/extract_crossfit_sentinels.py"
[[ -f "$PY_HELPER" ]] || { echo "[crossfit_direct] FATAL: Missing $PY_HELPER" >&2; exit 2; }

export PYTHONUNBUFFERED=1
python3 -u "$PY_HELPER" \
  --assoc_dir "$assoc_dir" \
  --magma_dir "$bundle" \
  --alpha "$alpha" \
  ${genera_csv:+--genera "$genera_csv"} \
  --outdir "$outdir" \
  2>&1 | tee "$outdir/crossfit_log.txt"

[[ -s "$outdir/sentinels_crossfit.tsv" ]] || { echo "[crossfit_direct] ERROR: sentinels_crossfit.tsv empty" >&2; exit 3; }
[[ -s "$outdir/beta_cf_by_species.tsv" ]] || { echo "[crossfit_direct] ERROR: beta_cf_by_species.tsv empty" >&2; exit 3; }

cp -f "$bundle/magma_sig_species.tsv" "$outdir/sig_species.tsv" || true
log "DONE at $(ts)"
