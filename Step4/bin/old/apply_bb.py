#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--bact", required=True)  # correlations_by_gene__Bacteroides.tsv
ap.add_argument("--prev", required=True)  # correlations_by_gene__Prevotella.tsv
ap.add_argument("--genes", required=True) # sig_genes_union.tsv (GENE_ID)
ap.add_argument("--sizes", required=True) # screen_sizes.txt with m_total,S_selected
ap.add_argument("--alpha", type=float, default=0.10)
ap.add_argument("--out", default="correlations_bb.tsv")
args = ap.parse_args()

def read_corr(p, genus):
    if not os.path.exists(p): 
        return pd.DataFrame(columns=["GENE_ID","p_spearman","p_pgls"])
    df = pd.read_csv(p, sep="\t")
    # assume columns: GENE_ID | ... | p_spearman | p_pgls
    keep = ["GENE_ID","p_spearman","p_pgls"]
    miss = [c for c in keep if c not in df.columns]
    if miss: 
        raise SystemExit(f"[BB] Missing columns in {p}: {miss}")
    df = df[keep].copy()
    df["genus"] = genus
    return df

B = read_corr(args.bact, "Bacteroides")
P = read_corr(args.prev, "Prevotella")
D = pd.concat([B,P], ignore_index=True)

# restrict to union genes
G = pd.read_csv(args.genes, sep="\t")["GENE_ID"].astype(str).tolist()
D = D[D["GENE_ID"].astype(str).isin(G)].copy()

# read m and S
sizes = dict([line.strip().split("\t") for line in open(args.sizes)])
m = int(float(sizes.get("m_total", "0"))); S = int(float(sizes.get("S_selected","0")))
alpha_prime = args.alpha * (m / S) if (S>0 and m>0) else args.alpha

# within-family BH (k<=2 tests per gene). We'll apply BH over the two p-values of *one* chosen endpoint.
# Use PGLS as primary; Spearman as supportive. You can adapt if you prefer a minP family-level test.
endpoint = "p_pgls"
out=[]
for gid, sub in D.groupby("GENE_ID", sort=False):
    ps = np.sort(sub[endpoint].values)
    k  = len(ps)
    if k==0: continue
    # BH within family: p_(i) <= alpha' * i / k
    sig = np.zeros(k, dtype=bool)
    thr = [alpha_prime * (i+1) / k for i in range(k)]
    for i,pv in enumerate(ps):
        if pv <= thr[i]:
            sig[:i+1] = True
    # map back to rows
    sub = sub.sort_values(endpoint).copy()
    sub["BB_alpha"] = alpha_prime
    sub["BB_sig_"+endpoint] = sig
    out.append(sub)

R = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["GENE_ID","genus",endpoint,"BB_alpha","BB_sig_"+endpoint])
R.to_csv(args.out, sep="\t", index=False)

