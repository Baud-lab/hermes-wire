#!/usr/bin/env python3
import argparse, os, re, glob, math, time
import numpy as np
import pandas as pd
from statistics import NormalDist
from statsmodels.stats.multitest import fdrcorrection

ap = argparse.ArgumentParser()
ap.add_argument("--gene-map", required=True)
ap.add_argument("--alpha", type=float, default=0.10)
ap.add_argument("--combiner", choices=["stouffer","empirical_brown"], default="empirical_brown")
ap.add_argument("--min-support", type=int, default=2)
ap.add_argument("--out-all", default="magma_all_with_q.tsv")
ap.add_argument("--out-sig", default="magma_sig.tsv")
args = ap.parse_args()

def ts(): return time.strftime("%Y-%m-%d %H:%M:%S")
def info(m): print(f"{ts()} [INFO] {m}", flush=True)
def warn(m): print(f"{ts()} [WARNING] {m}", flush=True)

P_COL_PRIORITY = ["P","P_JOINT","P_MULTI","P_SNPWISE_MEAN","P_SNPWISE_TOP1"]
SAFE_MIN = 1e-300
ND = NormalDist()
_trait_re_assoc = re.compile(r"_assoc$")
_trait_re_genus = re.compile(r"__([^_]+)")
def trait_to_genus(t):
    t = _trait_re_assoc.sub("", t)
    m = _trait_re_genus.search(t)
    return m.group(1) if m else t.split("_")[0]

def read_magma_geneout(path):
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python", comment="#")
    except Exception as e:
        warn(f"Failed to read {os.path.basename(path)}: {e}"); return None, None
    if "GENE" not in df.columns:
        warn(f"{os.path.basename(path)} has no 'GENE' column"); return None, None
    used = next((c for c in P_COL_PRIORITY if c in df.columns), None)
    if not used:
        warn(f"{os.path.basename(path)} lacks p col ({', '.join(P_COL_PRIORITY)})"); return None, None
    P = pd.to_numeric(df[used], errors="coerce").clip(lower=SAFE_MIN).where(lambda x: x<=1.0)
    keep = P.notna()
    out = pd.DataFrame({"GENE": df.loc[keep, "GENE"].astype(str).values, "P": P.loc[keep].astype(float).values})
    info(f"{os.path.basename(path)} -> using '{used}', kept {int(keep.sum())} rows")
    return out, used

def stouffer_p(ps):
    z = np.sum(ND.inv_cdf(1.0 - ps)) / math.sqrt(len(ps))
    return float(1.0 - ND.cdf(z))

def chi2_sf_wh(x, v):
    if v <= 0: return 1.0
    y = (x / v) ** (1.0/3.0)
    mu = 1.0 - 2.0/(9.0*v)
    sd = math.sqrt(2.0/(9.0*v))
    z = (y - mu) / sd
    return float(1.0 - ND.cdf(z))

def empirical_brown_p(ps, R):
    m = len(ps)
    S = float(np.sum(-2.0 * np.log(np.clip(ps, SAFE_MIN, 1.0))))
    if m == 1:
        return min(1.0, math.exp(-S/2.0))
    Rm = R[:m, :m]
    tri = np.triu_indices(m, k=1)
    rho = Rm[tri]
    poly = 3.263*rho + 0.710*(rho**2) + 0.027*(rho**3)
    varS = 4.0*m + 2.0*np.sum(poly)
    muS  = 2.0*m
    c = varS / (2.0*muS)
    v = 2.0 * (muS**2) / varS
    x = S / c
    return float(chi2_sf_wh(x, v))

def build_genus_corr(D_sub):
    W = D_sub.pivot_table(index="GENE_ID", columns="trait", values="P", aggfunc="min")
    if W.shape[1] < 2:
        traits = list(W.columns)
        return pd.DataFrame(np.eye(len(traits)), index=traits, columns=traits)
    X = -np.log(W.clip(lower=SAFE_MIN, upper=1.0))
    R = X.corr(method="pearson", min_periods=2).fillna(0.0)
    np.fill_diagonal(R.values, 1.0)
    return R

def main():
    print("[phase4] Collector (per-genus combiner + BH)")
    info(f"alpha={args.alpha} combiner={args.combiner} min_support={args.min_support}")
    info(f"outputs: all={args.out_all} sig={args.out_sig}")

    gmap = pd.read_csv(args.gene_map, sep="\t")
    info(f"Loaded gene map: {len(gmap)} rows")

    files = sorted(glob.glob("*.genes.genes.out"))
    info(f"Found {len(files)} MAGMA outputs")
    rows = []
    for p in files:
        trait = os.path.basename(p).replace(".genes.genes.out","")
        df, _ = read_magma_geneout(p)
        if df is None or df.empty: continue
        df["trait"] = trait
        df["genus"] = trait_to_genus(trait)
        rows.append(df)
    if not rows:
        empty(); return

    D = pd.concat(rows, ignore_index=True)
    D = D.rename(columns={"GENE": "GENE_ID"}).merge(gmap, on="GENE_ID", how="left")

    out_rows = []
    for gn, sub in D.groupby("genus", sort=False):
        R = build_genus_corr(sub)
        for gid, gsub in sub.groupby("GENE_ID", sort=False):
            ps = gsub["P"].astype(float).values
            traits = gsub["trait"].tolist()
            if len(ps) < args.min_support: continue
            if args.combiner == "stouffer":
                p_comb = stouffer_p(ps)
            else:
                R_sub = R.loc[traits, traits].values
                p_comb = empirical_brown_p(ps, R_sub)
            i_min = int(np.argmin(ps))
            out_rows.append({
                "genus": gn,
                "GENE_ID": gid,
                "GENE": (gsub["GENE"].iloc[0] if "GENE" in gsub.columns else gid),
                "P": float(p_comb),
                "n_traits": int(len(ps)),
                "min_p": float(ps[i_min]),
                "trait": traits[i_min],
                "top_trait": traits[i_min]
            })

    if not out_rows:
        empty(); return

    G = pd.DataFrame(out_rows)
    blocks = []
    for gn, sub in G.groupby("genus", sort=False):
        rej, qv = fdrcorrection(sub["P"].values, alpha=args.alpha)
        t = sub.copy(); t["q"] = qv; t["Significant"] = rej
        info(f"Genus {gn}: tested={len(t)} sig={int(rej.sum())} at q<{args.alpha}")
        blocks.append(t)
    R = pd.concat(blocks, ignore_index=True)
    R = R[["genus","GENE_ID","GENE","P","q","Significant","n_traits","min_p","trait","top_trait"]]
    R.to_csv(args.out_all, sep="\t", index=False)
    R[R["Significant"]].to_csv(args.out_sig, sep="\t", index=False)
    info(f"Wrote {args.out_all} and {args.out_sig}")

def empty():
    cols = ["genus","GENE_ID","GENE","P","q","Significant","n_traits","min_p","trait","top_trait"]
    pd.DataFrame(columns=cols).to_csv(args.out_all, sep="\t", index=False)
    pd.DataFrame(columns=cols).to_csv(args.out_sig, sep="\t", index=False)

if __name__ == "__main__":
    main()
