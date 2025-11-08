#!/usr/bin/env python3
import glob, os, re, numpy as np, pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--alpha", type=float, default=0.10)
ap.add_argument("--out-all", default="genus_screen_all.tsv")
ap.add_argument("--out-pass", default="genus_screen_pass.tsv")
args = ap.parse_args()

P_COLS = ["P","P_MULTI","P_JOINT","P_SNPWISE_MEAN","P_SNPWISE_TOP1"]
SAFE = 1e-300

def read_geneout(p):
    df = pd.read_csv(p, sep=r"\s+")
    pcol = next((c for c in P_COLS if c in df.columns), None)
    if pcol is None: return None
    g = df[["GENE", pcol]].copy()
    g.columns = ["GENE_ID","P"]
    g["P"] = pd.to_numeric(g["P"], errors="coerce").clip(lower=SAFE, upper=1.0)
    g = g[g["P"].notna()]
    trait = os.path.basename(p).replace(".genes.genes.out","")
    g["trait"] = trait
    return g

files = sorted(glob.glob("*.genes.genes.out"))
if len(files)==0:
    pd.DataFrame(columns=["GENE_ID","P_simes","q","n_traits"]).to_csv(args.out_all, sep="\t", index=False)
    pd.DataFrame(columns=["GENE_ID","q"]).to_csv(args.out_pass, sep="\t", index=False)
    raise SystemExit

D = pd.concat([x for f in files if (x:=read_geneout(f)) is not None], ignore_index=True)
# Expect exactly two genus traits; still handle 1..2 gracefully
rows=[]
for gid, sub in D.groupby("GENE_ID", sort=False):
    ps = np.sort(sub["P"].values)
    m  = len(ps)
    if m==0: continue
    # Simes p = min_i { m * p_(i) / i }
    p_simes = np.min([ps[i]*m/(i+1) for i in range(m)])
    rows.append({"GENE_ID":gid, "P_simes":float(min(1.0, p_simes)), "n_traits":m})

G = pd.DataFrame(rows)
if len(G)==0:
    G.to_csv(args.out_all, sep="\t", index=False)
    G.to_csv(args.out_pass, sep="\t", index=False)
    raise SystemExit

rej, q = fdrcorrection(G["P_simes"].values, alpha=args.alpha)
G["q"] = q; G["Significant"] = rej
G.to_csv(args.out_all, sep="\t", index=False)
G[G["Significant"]].to_csv(args.out_pass, sep="\t", index=False)

