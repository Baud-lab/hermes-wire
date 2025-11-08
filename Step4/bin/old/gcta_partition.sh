#!/usr/bin/env bash
# Purpose: Partition heritability using LIMIX (multi-kernel variance decomposition)
#          Kernels: crypt, rest, dam SRM, cage SRM (+ IID noise)
# Inputs:  PLINK base, SNP lists, SRM GRMs (.grm.gz + .grm.id), srm_ids.txt, phenotypes from RData
# Output:  greml_summary.tsv (trait, V_crypt, SE_crypt, ..., Vp, h2_SNP)

set -euo pipefail
echo "[limix_partition] START"

# ---------------- Parse args ----------------
PLINK=; PHENORDA=; CRYPT=; REST=; DAM_GZ=; DAM_ID=; CAGE_GZ=; CAGE_ID=; IDS=; OUT=; GENERA=; TPREF=
while [[ $# -gt 0 ]]; do
  case $1 in
    --plink_base)           PLINK=$2; shift 2;;
    --pheno_rda)            PHENORDA=$2; shift 2;;
    --crypt_snps)           CRYPT=$2; shift 2;;
    --rest_snps)            REST=$2; shift 2;;
    --srm_dam_grm_gz)       DAM_GZ=$2; shift 2;;
    --srm_dam_grm_id)       DAM_ID=$2; shift 2;;
    --srm_cage_grm_gz)      CAGE_GZ=$2; shift 2;;
    --srm_cage_grm_id)      CAGE_ID=$2; shift 2;;
    --ids_file)             IDS=$2; shift 2;;
    --outdir)               OUT=$2; shift 2;;
    --genera)               GENERA=$2; shift 2;;
    --trait-prefix|--trait_prefix) TPREF=$2; shift 2;;
    *) echo "Unknown arg $1" >&2; exit 1;;
  esac
done

OUT=$(realpath "$OUT")
mkdir -p "$OUT/phenos"

echo "[limix_partition] PLINK base: $PLINK"
echo "[limix_partition] SNP lists: crypt=$(wc -l < "$CRYPT" 2>/dev/null || echo 0), rest=$(wc -l < "$REST" 2>/dev/null || echo 0)"
echo "[limix_partition] SRM GRMs: dam=($DAM_GZ $DAM_ID) cage=($CAGE_GZ $CAGE_ID)"

# ---------------- R: write per-trait .pheno files ----------------
echo "[limix_partition] Writing per-trait phenotype files..."
Rscript - "${PHENORDA}" "${IDS}" "${OUT}/phenos" "${GENERA}" "${TPREF}" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
pheno_rda <- args[1]; ids_file <- args[2]; outdir <- args[3]
genera <- strsplit(args[4], ",", fixed=TRUE)[[1]]
tpref  <- args[5]

load(pheno_rda)  # residuals or residuals_qned_counts_objs
ids <- readLines(ids_file)

if (exists("residuals_qned_counts_objs")) {
  residuals <- residuals_qned_counts_objs[[length(residuals_qned_counts_objs)]]
}
if (!exists("residuals")) stop("Object 'residuals' not found in RData")

colnames(residuals) <- sub("_.*", "", colnames(residuals))
residuals <- t(residuals)

traits_all <- colnames(residuals)
traits_keep <- character(0)
for (g in genera) {
  sel <- traits_all[grepl(paste0(tpref,g), traits_all)]
  if (g == "Bacteroides") sel <- sel[!grepl("s__Bacteroides_F", sel)]
  traits_keep <- c(traits_keep, sel)
}
traits_keep <- unique(traits_keep)
if (!length(traits_keep)) stop("No traits matched genera/prefix filters.")

residuals <- residuals[ids, traits_keep, drop=FALSE]

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
n <- 0L
for (tr in colnames(residuals)) {
  ph <- data.frame(FID=ids, IID=ids, PHE=residuals[,tr])
  if (all(is.na(ph$PHE)) || length(unique(na.omit(ph$PHE))) < 3) next
  fn <- file.path(outdir, paste0(tr, ".pheno"))
  write.table(ph, file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  n <- n + 1L
}
cat("[limix_partition/R] Traits written:", n, "\n")
RS

# ---------------- Python: build kernels & fit VarDec per trait ----------------
python3 - "$PLINK" "$CRYPT" "$REST" "$DAM_GZ" "$DAM_ID" "$CAGE_GZ" "$CAGE_ID" "$IDS" "$OUT" <<'PY'
import os, sys, gzip, numpy as np, pandas as pd
from limix.io import plink as limix_plink
from limix.stats import linear_kinship
from limix.vardec import VarDec
from limix.qc import normalise_covariance

plink_base, crypt_file, rest_file, dam_gz, dam_id, cage_gz, cage_id, ids_file, outdir = sys.argv[1:]
ids = pd.read_csv(ids_file, header=None, squeeze=True, dtype=str).tolist()

# Read PLINK
bim, fam, bed = limix_plink.read(plink_base, verbose=False)
fam['iid'] = fam['iid'].astype(str)
fam = fam.set_index('iid')

keep = [i for i in ids if i in fam.index]
if len(keep) < len(ids):
    missing = sorted(set(ids) - set(keep))
    print(f"[limix_partition] WARNING: dropping {len(missing)} IDs not in PLINK (e.g. {missing[:5]})", file=sys.stderr)

sidx = fam.loc[keep, 'i'].astype(int).values  # indices into bed axis=1

def read_snps(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return np.array([], dtype=int)
    s = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0], dtype=str).iloc[:,0].astype(str).values
    hit = bim[bim['snp'].astype(str).isin(s)]
    return hit['i'].astype(int).values

crypt_idx = read_snps(crypt_file)
rest_idx  = read_snps(rest_file)

def make_K(sub_idx):
    if sub_idx.size == 0:
        return None
    # bed shape: variants x samples (dask); select then compute to numpy
    X = bed[sub_idx, :][:, sidx].compute().T  # n x m
    K = linear_kinship(X, verbose=False)
    K = normalise_covariance(K)  # mean(diag)=1
    K[np.diag_indices_from(K)] += 1e-6
    return K

K_crypt = make_K(crypt_idx)
K_rest  = make_K(rest_idx)

def read_grm_gz(grm_gz, grm_id, order_keep):
    # GCTA text .grm.gz with lower-tri; last column is value
    ids_df = pd.read_csv(grm_id, sep=r"\s+", header=None, names=["fid","iid"])
    ids_df['iid'] = ids_df['iid'].astype(str)
    base = ids_df['iid'].tolist()
    n = len(base)
    K = np.zeros((n,n), float)
    with gzip.open(grm_gz, 'rt') as f:
        df = pd.read_csv(f, sep=r"\s+", header=None)
    val_col = df.columns[-1]
    for i,j,v in zip(df.iloc[:,0].astype(int), df.iloc[:,1].astype(int), df.iloc[:,val_col].astype(float)):
        i -= 1; j -= 1
        K[i,j] = v; K[j,i] = v
    # Reorder to 'keep'
    pos = [base.index(i) for i in order_keep]
    K = K[np.ix_(pos, pos)]
    K = normalise_covariance(K)
    K[np.diag_indices_from(K)] += 1e-6
    return K

K_dam  = read_grm_gz(dam_gz,  dam_id,  keep)
K_cage = read_grm_gz(cage_gz, cage_id, keep)

phenodir = os.path.join(outdir, "phenos")
rows = []
for fn in os.listdir(phenodir):
    if not fn.endswith(".pheno"): 
        continue
    trait = os.path.splitext(fn)[0]
    df = pd.read_csv(os.path.join(phenodir, fn), sep="\t")
    df['IID'] = df['IID'].astype(str)
    df = df.set_index('IID').loc[keep]
    y = df['PHE'].values.astype(float)
    M = np.ones((len(y), 1))  # intercept only
    vd = VarDec(y, "normal", M)
    if K_crypt is not None: vd.append(K_crypt, name="crypt")
    if K_rest  is not None: vd.append(K_rest,  name="rest")
    vd.append(K_dam, name="dam")
    vd.append(K_cage, name="cage")
    vd.append_iid(name="noise")
    vd.fit(verbose=False)

    # Extract component scales robustly
    try:
        comps = vd.covariance   # list-like of components
        names = [c.name for c in comps]
        scales= [float(getattr(c, "scale")) for c in comps]
        d = dict(zip(names, scales))
    except Exception as e:
        print(f"[limix_partition] WARN: could not read component scales ({e}); filling NA)", file=sys.stderr)
        d = {}

    V_crypt = float(d.get("crypt", 0.0 if K_crypt is not None else np.nan))
    V_rest  = float(d.get("rest",  0.0 if K_rest  is not None else np.nan))
    V_dam   = float(d.get("dam",   np.nan))
    V_cage  = float(d.get("cage",  np.nan))
    V_noise = float(d.get("noise", np.nan))
    Vp = np.nansum([V_crypt, V_rest, V_dam, V_cage, V_noise])
    h2 = (0.0 if np.isnan(Vp) or Vp==0 else (np.nan_to_num(V_crypt) + np.nan_to_num(V_rest)) / Vp)

    rows.append({
      "trait": trait,
      "V_crypt": V_crypt, "SE_crypt": np.nan,
      "V_rest":  V_rest,  "SE_rest":  np.nan,
      "V_dam":   V_dam,   "SE_dam":   np.nan,
      "V_cage":  V_cage,  "SE_cage":  np.nan,
      "Vp": Vp, "h2_SNP": h2
    })

out = pd.DataFrame(rows).sort_values("trait")
out.to_csv(os.path.join(outdir, "greml_summary.tsv"), sep="\t", index=False)
print(f"[limix_partition] Wrote greml_summary.tsv with {len(out)} traits.")
PY

echo "[limix_partition] DONE"
