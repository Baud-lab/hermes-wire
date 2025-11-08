# save as quick_genus_preview.py
# usage:  python3 -u quick_genus_preview.py --q 0.10 --gene-model multi --combiner empirical_brown
# notes : uses P_MULTI (multi) | P_SNPWISE_TOP1 (top) | P_SNPWISE_MEAN (mean)
import os, re, glob, math, argparse
import numpy as np, pandas as pd

try:
    from scipy.stats import chi2
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

def log(msg): print(msg, flush=True)

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
    q = p*m/r
    qs = q[o]; qs = np.minimum.accumulate(qs[::-1])[::-1]
    out = np.empty_like(qs); out[o.argsort()] = np.clip(qs,0,1)
    return out

# Empirical Brown's Method (Poole et al. 2016). Falls back to Simes if ill-conditioned or SciPy missing.
def ebm_p_for_gene(pvals, cov_traits):
    pvals = np.asarray(pvals, float)
    traits = cov_traits.columns
    ok = np.isfinite(pvals) & (pvals>0) & (pvals<=1)
    if ok.sum() < 2:  # need at least 2 signals
        return simes_p(pvals[ok])
    use = np.where(ok)[0]
    S = float(np.sum(-2.0*np.log(pvals[use])))
    m = len(use)
    mu = 2.0*m
    Σ = np.array(cov_traits.values[np.ix_(use, use)])
    varS = float(np.sum(Σ))
    if not HAS_SCIPY or not np.isfinite(varS) or varS <= 0:
        return simes_p(pvals[ok])
    nu = 2.0*(mu**2)/varS
    c  = varS/(2.0*mu)
    x  = S/c
    if not all(np.isfinite([nu, x])) or nu <= 0:
        return simes_p(pvals[ok])
    return float(np.clip(1.0 - chi2.cdf(x, df=nu), 0.0, 1.0))

ap = argparse.ArgumentParser()
ap.add_argument("--q", type=float, default=0.10, help="BH FDR across families")
ap.add_argument("--gene-model", choices=["multi","top","mean"], default="multi",
               help="Which gene p-value column to use from MAGMA output")
ap.add_argument("--combiner", choices=["empirical_brown","simes"], default="empirical_brown",
               help="Combine species (traits) within genus")
args = ap.parse_args()

# pick p-value column based on gene model
PMAP = {
    "multi": "P_MULTI",
    "top":   "P_SNPWISE_TOP1",
    "mean":  "P_SNPWISE_MEAN",
}
pcol_wanted = PMAP[args.gene_model]

files = sorted(glob.glob("*.genes.genes.out"))
log(f"[preview] found {len(files)} MAGMA gene files")
if not files:
    raise SystemExit("[preview] ERROR: no *.genes.genes.out files here")

rows=[]
for i, path in enumerate(files, 1):
    log(f"[preview] [{i}/{len(files)}] reading {path}")
    df = pd.read_csv(path, delim_whitespace=True, comment="#", engine="python")
    cols = set(df.columns)
    # Robust header handling
    gene_col = "GENE" if "GENE" in cols else ("Gene" if "Gene" in cols else None)
    if gene_col is None:
        log("  !! no GENE column; skipping")
        continue
    # choose p column
    pcol = None
    for c in [pcol_wanted, "P", "p", "pval", "pvalue"]:
        if c in cols:
            pcol = c; break
    if pcol is None:
        log(f"  !! none of {[pcol_wanted,'P','p','pval','pvalue']} in columns; skipping")
        continue
    mini = df[[gene_col, pcol]].copy()
    mini.columns = ["GENE_ID","P"]
    mini["trait"] = os.path.basename(path).replace(".genes.genes.out","")
    mini["P"] = pd.to_numeric(mini["P"], errors="coerce")
    mini = mini[np.isfinite(mini["P"]) & (mini["P"]>0) & (mini["P"]<=1)]
    rows.append(mini)
    log(f"  -> kept {len(mini)} rows using column {pcol}")

if not rows:
    raise SystemExit("[preview] ERROR: couldn't parse any gene/p rows from outputs")

all_df = pd.concat(rows, ignore_index=True)
log(f"[preview] total parsed rows: {len(all_df)}")

# Optional pretty names
if os.path.exists("gene_name_map.tsv"):
    try:
        gmap = pd.read_csv("gene_name_map.tsv", sep="\t")
        if {"GENE_ID","GENE"}.issubset(gmap.columns):
            all_df = all_df.merge(gmap, how="left", on="GENE_ID")
            all_df["GENE"] = all_df["GENE"].fillna(all_df["GENE_ID"])
        else:
            all_df["GENE"] = all_df["GENE_ID"]
    except Exception as e:
        log(f"[preview] WARN: couldn't read gene_name_map.tsv: {e}")
        all_df["GENE"] = all_df["GENE_ID"]
else:
    all_df["GENE"] = all_df["GENE_ID"]

all_df["genus"] = all_df["trait"].map(trait_to_genus)

# --------- Build covariance per genus for EBM (over genes; columns = traits) ---------
cov_by_genus = {}
if args.combiner == "empirical_brown" and HAS_SCIPY:
    log("[preview] computing genus-level covariance for Empirical Brown ...")
    for gn, sub in all_df.groupby("genus", sort=False):
        wide = sub.pivot_table(index="GENE_ID", columns="trait", values="P", aggfunc="min")
        # transform and compute covariance across traits
        Y = -2.0*np.log(wide)
        cov = Y.cov(min_periods=5)  # needs enough genes
        # stabilize diagonals
        for i in range(cov.shape[0]): 
            if not np.isfinite(cov.iat[i,i]) or cov.iat[i,i] <= 1e-8: cov.iat[i,i] = 1e-8
        cov_by_genus[gn] = cov
        log(f"  - {gn}: traits={cov.shape[0]}  genes={Y.shape[0]}")
else:
    if args.combiner == "empirical_brown":
        log("[preview] SciPy not available; falling back to Simes for combining")

# --------- Combine per (GENE, genus) ----------
fam_rows=[]
for (gid, gname, gn), sub in all_df.groupby(["GENE_ID","GENE","genus"], sort=False):
    p = pd.to_numeric(sub["P"], errors="coerce")
    p = p[np.isfinite(p) & (p>0) & (p<=1)]
    if len(p)==0: continue
    if args.combiner == "empirical_brown" and HAS_SCIPY and gn in cov_by_genus and cov_by_genus[gn].shape[0] >= 2:
        # map p-values to trait order used in covariance
        traits = cov_by_genus[gn].columns.tolist()
        vec = []
        for tr in traits:
            val = sub.loc[sub["trait"]==tr, "P"]
            vec.append(val.iloc[0] if len(val) else np.nan)
        pstar = ebm_p_for_gene(np.array(vec, float), cov_by_genus[gn])
        tag = "ebm"
    else:
        pstar = simes_p(p.values)
        tag = "simes"
    fam_rows.append({"GENE_ID":gid,"GENE":gname,"genus":gn,"combinedP":pstar,"n_species":int(len(p)),"combiner":tag})

fams = pd.DataFrame(fam_rows)
if fams.empty:
    raise SystemExit("[preview] ERROR: no combined rows (all P missing?)")

# --------- BH across families, then BB within selected families ----------
fams["q"] = bh_q(fams["combinedP"].values)
fams["Significance"] = fams["q"] < args.q
fams.sort_values(["genus","q","combinedP"], ascending=[True,True,True], inplace=True)

os.makedirs("family_fdr_by_genus", exist_ok=True)
os.makedirs("species_calls_by_genus", exist_ok=True)

fams.to_csv("magma_all_genus_with_q.tsv", sep="\t", index=False)
fams[fams["Significance"]].to_csv("magma_sig_genus.tsv", sep="\t", index=False)

R = int(fams["Significance"].sum()); M = int(len(fams))
q_within = args.q * (R / M) if M>0 else 0.0
log(f"[preview] BH across families at q={args.q}: selected R={R} of M={M} families; BB within-family level={q_within:.4g}")

# Within-genus species calls per selected family (BB)
def within_family_calls(genus, fam_df, all_df, q_within):
    sel = fam_df.loc[fam_df["Significance"], "GENE_ID"].tolist()
    calls=[]
    for gid in sel:
        sp = all_df[(all_df["genus"]==genus) & (all_df["GENE_ID"]==gid)].copy()
        if sp.empty: continue
        # BH per family at q_within
        _, qv = np.unique(sp["P"].rank(method="min"), return_inverse=True)  # avoid confusion; use statsmodels if available
        # simple BH implementation:
        p = sp["P"].values
        o = np.argsort(p); r = np.arange(1, len(p)+1)
        q = p[o] * len(p) / r
        q = np.minimum.accumulate(q[::-1])[::-1]
        sp = sp.iloc[o].copy()
        sp["q_species"] = q
        sp["Selected_in_family"] = sp["q_species"] < q_within
        calls.append(sp[["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"]])
    return pd.concat(calls, ignore_index=True) if calls else pd.DataFrame(columns=["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"])

# write per-genus side tables
for gn, fam_sub in fams.groupby("genus", sort=False):
    fam_sub.to_csv(f"family_fdr_by_genus/{gn}_families.tsv", sep="\t", index=False)
    species_calls = within_family_calls(gn, fam_sub, all_df, q_within) if q_within>0 else pd.DataFrame(columns=["genus","GENE","GENE_ID","trait","P","q_species","Selected_in_family"])
    species_calls.sort_values(["GENE","q_species","P"], inplace=True, ascending=[True,True,True], na_position="last")
    species_calls.to_csv(f"species_calls_by_genus/{gn}_species_calls.tsv", sep="\t", index=False)

log("[preview] wrote magma_all_genus_with_q.tsv and magma_sig_genus.tsv")
log("[preview] done")
