# quick_genus_preview_ebm.py
# Empirical Brown’s Method (EBM) genus-level collector with verbose logging + heartbeat.
# Uses MAGMA gene files (*.genes.genes.out), picks P_MULTI for --gene-model multi (or top/mean).
# Outputs:
#   magma_all_genus_with_q.tsv
#   magma_sig_genus.tsv
#   family_fdr_by_genus/<genus>_families.tsv
#   species_calls_by_genus/<genus>_species_calls.tsv
import os, re, glob, sys, time, argparse
import numpy as np, pandas as pd

LOG = open("preview_ebm_debug.log","w", buffering=1)
def log(msg): LOG.write(msg+"\n")
def beat():   # heartbeat to stdout
    try: print(".", end="", flush=True)
    except Exception: pass

# -------- helpers --------
def trait_to_genus(t):
    t = re.sub(r"_assoc$", "", t)
    m = re.search(r"__([^_]+)", t)
    return m.group(1) if m else t.split("_")[0]

def simes_p(pvals):
    p = np.sort(np.asarray(pvals, float)); m = len(p)
    if m == 0: return 1.0
    return float(np.minimum(1.0, np.min(p * (m / np.arange(1, m+1)))))

def bh_q(p):
    p = np.asarray(p, float); m = len(p)
    o = np.argsort(p); r = np.empty(m,int); r[o]=np.arange(1,m+1)
    q = p[o] * m / r
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[o.argsort()] = np.clip(q,0,1)
    return out

# Empirical Brown’s Method (Poole et al., 2016): combine dependent p across traits
try:
    from scipy.stats import chi2
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

def ebm_from_vector(pvec, cov_sub):
    """pvec ordered to match cov_sub columns; NaNs allowed (ignored)."""
    p = np.asarray(pvec, float)
    ok = np.isfinite(p) & (p>0) & (p<=1)
    m  = int(ok.sum())
    if (not HAS_SCIPY) or m < 2:            # need >=2 signals + SciPy
        return simes_p(p[ok]), "simes_fallback"
    p = p[ok]
    S  = float(np.sum(-2.0*np.log(p)))
    mu = 2.0*m
    # take matching submatrix of cov for the ok traits
    Σ  = np.asarray(cov_sub)[ok][:, ok]
    varS = float(np.sum(Σ))
    if (not np.isfinite(varS)) or varS <= 0:
        return simes_p(p), "simes_badvar"
    nu = 2.0*(mu**2)/varS
    c  = varS/(2.0*mu)
    x  = S/c
    if (not np.isfinite(nu)) or nu <= 0 or (not np.isfinite(x)):
        return simes_p(p), "simes_num"
    return float(np.clip(1.0 - chi2.cdf(x, df=nu), 0.0, 1.0)), "ebm"

# -------- CLI --------
ap = argparse.ArgumentParser()
ap.add_argument("--q", type=float, default=0.10, help="BH FDR across families")
ap.add_argument("--gene-model", choices=["multi","top","mean"], default="multi",
               help="Which MAGMA gene p-value column to use (multi uses P_MULTI)")
args = ap.parse_args()

PMAP = {"multi": "P_MULTI", "top": "P_SNPWISE_TOP1", "mean": "P_SNPWISE_MEAN"}
PCOL = PMAP[args.gene_model]

files = sorted(glob.glob("*.genes.genes.out"))
log(f"[EBM] found {len(files)} MAGMA files")
if not files:
    print("\n[EBM] ERROR: no *.genes.genes.out here", flush=True); sys.exit(1)
if not HAS_SCIPY:
    print("\n[EBM] SciPy not found — will fall back to Simes where needed.", flush=True)

# -------- read + gather gene-level P by trait --------
rows=[]; t0=time.time()
for i, path in enumerate(files, 1):
    log(f"[read] {i}/{len(files)}  {path}")
    try:
        df = pd.read_csv(path, delim_whitespace=True, comment="#", engine="python")
    except Exception as e:
        log(f"  !! read failed: {e}"); beat(); continue
    cols=set(df.columns)
    gene_col = "GENE" if "GENE" in cols else ("Gene" if "Gene" in cols else None)
    if gene_col is None:
        log("  !! no GENE column; skip"); beat(); continue
    # preferred P column with fallbacks
    if   PCOL in cols:                pcol = PCOL
    elif "P_SNPWISE_TOP1" in cols:    pcol = "P_SNPWISE_TOP1"
    elif "P_SNPWISE_MEAN" in cols:    pcol = "P_SNPWISE_MEAN"
    elif "P" in cols:                 pcol = "P"
    else:
        log(f"  !! no usable P column; skip"); beat(); continue
    mini = df[[gene_col, pcol]].copy()
    mini.columns = ["GENE_ID","P"]
    mini["trait"] = os.path.basename(path).replace(".genes.genes.out","")
    mini["P"] = pd.to_numeric(mini["P"], errors="coerce")
    before=len(mini); mini = mini[np.isfinite(mini["P"]) & (mini["P"]>0) & (mini["P"]<=1)]
    log(f"  -> kept {len(mini)}/{before} using {pcol}")
    rows.append(mini); beat()
print("")  # newline after dots

if not rows:
    print("[EBM] ERROR: zero usable rows (unexpected headers?)", flush=True); sys.exit(2)

all_df = pd.concat(rows, ignore_index=True)
log(f"[EBM] total parsed rows: {len(all_df)}")

# optional pretty names
if os.path.exists("gene_name_map.tsv"):
    try:
        gmap = pd.read_csv("gene_name_map.tsv", sep="\t")
        if {"GENE_ID","GENE"}.issubset(gmap.columns):
            all_df = all_df.merge(gmap, how="left", on="GENE_ID")
            all_df["GENE"] = all_df["GENE"].fillna(all_df["GENE_ID"])
        else:
            all_df["GENE"] = all_df["GENE_ID"]
    except Exception as e:
        log(f"[EBM] WARN: couldn't read gene_name_map.tsv: {e}")
        all_df["GENE"] = all_df["GENE_ID"]
else:
    all_df["GENE"] = all_df["GENE_ID"]

all_df["genus"] = all_df["trait"].map(trait_to_genus)

# -------- covariance per genus (traits columns) --------
cov_by_genus = {}
log("[EBM] building per-genus covariance of Y = -2 ln P (rows=genes, cols=traits)")
for gn, sub in all_df.groupby("genus", sort=False):
    # one p per gene×trait → pivot genes x traits
    wide = sub.pivot_table(index="GENE_ID", columns="trait", values="P", aggfunc="min")
    if wide.shape[1] < 2 or wide.shape[0] < 10:
        log(f"  - {gn}: too few traits/genes for stable cov (genes={wide.shape[0]}, traits={wide.shape[1]}); EBM→Simes")
        cov_by_genus[gn] = None
        continue
    Y = -2.0*np.log(wide.clip(lower=np.finfo(float).tiny))   # avoid log(0)
    cov = Y.cov(min_periods=5)                               # DataFrame cov across traits
    # stabilize diagonals and replace non-finite
    cov = cov.replace([np.inf,-np.inf], np.nan)
    for i in range(cov.shape[0]):
        if not np.isfinite(cov.iat[i,i]) or cov.iat[i,i] <= 1e-8:
            cov.iat[i,i] = 1e-8
    cov = cov.fillna(0.0)
    cov_by_genus[gn] = cov
    log(f"  - {gn}: cov OK  genes={Y.shape[0]}  traits={Y.shape[1]}")

# -------- combine per (GENE, genus) with EBM (fallback Simes) ----------
fam_rows=[]
log("[EBM] combining per gene×genus ...")
for (gid, gname, gn), sub in all_df.groupby(["GENE_ID","GENE","genus"], sort=False):
    # ensure consistent trait order
    traits = list(sub["trait"].unique())
    # prefer cov’s trait order if available
    if (cov_by_genus.get(gn) is not None):
        trait_order = list(cov_by_genus[gn].columns)
    else:
        trait_order = sorted(traits)

    # build p-vector aligned to trait_order (NaN where missing)
    pvec = []
    sub_map = dict(zip(sub["trait"], pd.to_numeric(sub["P"], errors="coerce")))
    for tr in trait_order:
        pvec.append(sub_map.get(tr, np.nan))

    if cov_by_genus.get(gn) is not None and HAS_SCIPY:
        pstar, tag = ebm_from_vector(pvec, cov_by_genus[gn].loc[trait_order, trait_order])
    else:
        pstar, tag = simes_p(pvec), "simes"
    fam_rows.append({"GENE_ID":gid,"GENE":gname,"genus":gn,
                     "combinedP":pstar,"n_species":int(np.isfinite(pvec).sum()),
                     "combiner":tag})
    # light heartbeat per ~200 genes
    if len(fam_rows) % 200 == 0: beat()
print("")

fams = pd.DataFrame(fam_rows)
if fams.empty:
    print("[EBM] ERROR: no combined rows", flush=True); sys.exit(3)

# -------- BH across families + Benjamini–Bogomolov within families --------
fams["q"] = bh_q(fams["combinedP"].values)
fams["Significance"] = fams["q"] < args.q
fams.sort_values(["genus","q","combinedP"], ascending=[True,True,True], inplace=True)

os.makedirs("family_fdr_by_genus", exist_ok=True)
os.makedirs("species_calls_by_genus", exist_ok=True)
fams.to_csv("magma_all_genus_with_q.tsv", sep="\t", index=False)
fams[fams["Significance"]].to_csv("magma_sig_genus.tsv", sep="\t", index=False)

R = int(fams["Significance"].sum()); M = int(len(fams))
q_within = args.q * (R / M) if M>0 else 0.0
log(f"[EBM] BH across families q={args.q}: selected R={R} / M={M} → q_within={q_within:.4g}")

def within_family_calls(genus, fam_df, all_df, q_within):
    sel = fam_df.loc[fam_df["Significance"], "GENE_ID"].tolist()
    calls=[]
    for gid in sel:
        sp = all_df[(all_df["genus"]==genus) & (all_df["GENE_ID"]==gid)].copy()
        if sp.empty: continue
        p = pd.to_numeric(sp["P"], errors="coerce").values
        o = np.argsort(p); r = np.arange(1, len(p)+1)
        q = p[o] * len(p) / r
        q = np.minimum.accumulate(q[::-1])[::-1]
        sp = sp.iloc[o].copy()
        sp["q_species"] = q
        sp["Selected_in_family"] = sp["q_species"] < q_within
        calls.append(sp[["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"]])
    return pd.concat(calls, ignore_index=True) if calls else pd.DataFrame(columns=["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"])

for gn, fam_sub in fams.groupby("genus", sort=False):
    fam_sub.to_csv(f"family_fdr_by_genus/{gn}_families.tsv", sep="\t", index=False)
    species_calls = within_family_calls(gn, fam_sub, all_df, q_within) if q_within>0 else pd.DataFrame(columns=["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"])
    species_calls.sort_values(["GENE","q_species","P"], ascending=[True,True,True], inplace=True)
    species_calls.to_csv(f"species_calls_by_genus/{gn}_species_calls.tsv", sep="\t", index=False)
    log(f"[EBM] wrote species_calls for {gn} (n={len(species_calls)})")

elapsed = time.time() - t0
print(f"[EBM] DONE in {elapsed:.1f}s  | fams={len(fams)}  sig={R}  (see preview_ebm_debug.log)", flush=True)
