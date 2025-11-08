#!/usr/bin/env python3
import os, glob, pandas as pd
from statsmodels.stats.multitest import fdrcorrection

P_COLS = ["P","P_JOINT","P_MULTI","P_SNPWISE_MEAN","P_SNPWISE_TOP1"]

def read_one(p):
    df = pd.read_csv(p, sep=r"\s+", engine="python", comment="#")
    use = next((c for c in P_COLS if c in df.columns), None)
    if use is None: return None
    out = df[["GENE", use]].dropna().rename(columns={use:"P"})
    tr = os.path.basename(p).replace(".genes.genes.out","")
    out["trait"] = tr
    out["genus"] = "Prevotella" if "Prevotella" in tr else ("Bacteroides" if "Bacteroides" in tr else tr)
    return out

dfs=[]
for p in glob.glob("*.genes.genes.out"):
    d = read_one(p)
    if d is not None and len(d): dfs.append(d)

if not dfs:
    pd.DataFrame(columns=["genus","GENE_ID","GENE","P","q","Significant"]).to_csv("magma_genus_all_with_q.tsv", sep="\t", index=False)
    pd.DataFrame(columns=["genus","GENE_ID","GENE"]).to_csv("magma_genus_sig.tsv", sep="\t", index=False)
    pd.DataFrame(columns=["genus","GENE_ID","GENE"]).to_csv("magma_sig_species.tsv", sep="\t", index=False)
    raise SystemExit(0)

D = pd.concat(dfs, ignore_index=True).rename(columns={"GENE":"GENE_ID"})
blocks=[]
for gn, sub in D.groupby("genus", sort=False):
    rej, qv = fdrcorrection(sub["P"].values, alpha=float(os.environ.get("Q","0.10")))
    t = sub.copy(); t["q"]=qv; t["Significant"]=rej
    blocks.append(t)
R = pd.concat(blocks, ignore_index=True)

if os.path.exists("gene_name_map.tsv"):
    gmap = pd.read_csv("gene_name_map.tsv", sep="\t")
    R = R.merge(gmap, on="GENE_ID", how="left")
    R["GENE"] = R["GENE"].fillna(R["GENE_ID"])
else:
    R["GENE"] = R["GENE_ID"]

R[["genus","GENE_ID","GENE","P","q","Significant"]].to_csv("magma_genus_all_with_q.tsv", sep="\t", index=False)
sig = R[R["Significant"]][["genus","GENE_ID","GENE"]].drop_duplicates()
sig.to_csv("magma_genus_sig.tsv", sep="\t", index=False)
sig.to_csv("magma_sig_species.tsv", sep="\t", index=False)

