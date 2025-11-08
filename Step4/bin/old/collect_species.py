#!/usr/bin/env python3
"""
Combine MAGMA gene results across traits (species) WITHIN EACH GENUS, then do BH FDR per genus.

What this does
--------------
1) Read all '*.genes.genes.out' (one per trait/species) in CWD; ignore '#' comments.
2) For each file, pick the first available gene p-value column from:
      P, P_JOINT, P_MULTI, P_SNPWISE_MEAN, P_SNPWISE_TOP1
   Normalize it as column 'P'.
3) Add the 'trait' (file stem) and derive 'genus' from the trait name.
4) For each (genus, GENE_ID), COMBINE the set of P's across traits using the SIMES method:
      sort p_(1) <= ... <= p_(m);  p_simes = min_i (m/i)*p_(i), clipped to 1
   (Simes is valid and robust under positive dependence.)
5) Within each genus, apply BH FDR to the combined p’s and mark Significant = (q < alpha).
6) Write:
      magma_all.tsv  -> combined results for all (genus, GENE_ID)
      magma_sig.tsv  -> only rows with Significant == True

Notes
-----
* Filenames are generic (no "species" or "genes" words), as requested.
* We DO NOT drop tiny p-values. If a p <= 0 is encountered, we set it to the smallest positive float
  (machine epsilon) to avoid zeros that can break some transforms; this keeps “very small” p's intact.
"""

import os
import re
import glob
import time
import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection  # in your container

# ---------- small helpers ----------
def ts(): return time.strftime("%Y-%m-%d %H:%M:%S")
def info(msg): print(f"{ts()} [INFO] {msg}", flush=True)
def warn(msg): print(f"{ts()} [WARNING] {msg}", flush=True)
def err(msg):  print(f"{ts()} [ERROR] {msg}", flush=True)

# Columns we will accept from MAGMA outputs, in priority order
P_COL_PRIORITY = ["P", "P_JOINT", "P_MULTI", "P_SNPWISE_MEAN", "P_SNPWISE_TOP1"]

_trait_re_assoc = re.compile(r"_assoc$")
_trait_re_genus = re.compile(r"__([^_]+)")

def trait_to_genus(t: str) -> str:
    """ 's__Prevotella_spXXXX_assoc' -> 'Prevotella' """
    t = _trait_re_assoc.sub("", t)
    m = _trait_re_genus.search(t)
    return m.group(1) if m else t.split("_")[0]

def simes_p(pvals: np.ndarray) -> float:
    """
    Simes combined p-value. Assumes pvals are valid in (0,1].
    p_simes = min_i (m/i) * p_(i), where p_(i) are sorted ascending.
    """
    p = np.sort(pvals)
    m = float(len(p))
    idx = np.arange(1, len(p) + 1, dtype=float)
    ps = np.minimum(1.0, (m / idx) * p)
    return float(np.min(ps)) if len(ps) else 1.0

def read_geneout_one(path: str):
    """
    Read one MAGMA '.genes.genes.out', pick a usable p column, return tidy df with columns:
      ['GENE', 'P']
    """
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python", comment="#")
    except Exception as e:
        warn(f"Failed to read {os.path.basename(path)}: {e}")
        return None

    if "GENE" not in df.columns:
        warn(f"{os.path.basename(path)} missing 'GENE' column; skipping.")
        return None

    used = next((c for c in P_COL_PRIORITY if c in df.columns), None)
    if used is None:
        warn(f"{os.path.basename(path)} has no recognized p column ({', '.join(P_COL_PRIORITY)}); skipping.")
        return None

    # Coerce numeric p; keep tiny p's by clamping *zeros only* up to machine epsilon
    p = pd.to_numeric(df[used], errors="coerce")
    # Replace non-finite, negatives, and >1 with NaN; replace exact zeros with smallest positive float
    eps = np.finfo(float).tiny
    p = p.where(np.isfinite(p) & (p > 0) & (p <= 1.0), np.nan)
    p = p.fillna(0.0).replace(0.0, eps)  # exact zeros -> eps; NaNs were set to 0 then to eps
    keep = p.between(eps, 1.0, inclusive="both")

    out = pd.DataFrame({"GENE": df.loc[keep, "GENE"].astype(str).values,
                        "P":    p.loc[keep].astype(float).values})
    info(f"{os.path.basename(path)} -> using '{used}', kept {len(out)} rows.")
    return out

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-map", required=True, help="gene_name_map.tsv")
    ap.add_argument("--alpha", type=float, default=0.10, help="BH FDR threshold (within genus)")
    ap.add_argument("--out-all", default="magma_all.tsv", help="output: full combined table")
    ap.add_argument("--out-sig", default="magma_sig.tsv", help="output: significant subset")
    return ap.parse_args()

def main():
    args = parse_args()
    info("=== Collector (genus-combined) start ===")
    info(f"Alpha (FDR): {args.alpha:.4f}")
    info(f"Gene map   : {args.gene_map}")
    info(f"Outputs    : {args.out_all} | {args.out_sig}")

    # Load display names
    if not os.path.exists(args.gene_map):
        err("gene_name_map.tsv not found.")
        # write empty outputs for pipeline continuity
        cols = ["genus", "GENE_ID", "GENE", "P", "q", "Significant", "n_traits", "min_p"]
        pd.DataFrame(columns=cols).to_csv(args.out_all, sep="\t", index=False)
        pd.DataFrame(columns=cols).to_csv(args.out_sig, sep="\t", index=False)
        return

    gmap = pd.read_csv(args.gene_map, sep="\t")
    if not {"GENE_ID", "GENE"}.issubset(gmap.columns):
        err("gene_name_map.tsv must have columns: GENE_ID, GENE")
        return
    info(f"Loaded gene map with {len(gmap)} rows.")

    files = sorted(glob.glob("*.genes.genes.out"))
    info(f"Found {len(files)} MAGMA outputs.")
    if not files:
        warn("No .genes.genes.out files found; writing empty outputs.")
        cols = ["genus", "GENE_ID", "GENE", "P", "q", "Significant", "n_traits", "min_p"]
        pd.DataFrame(columns=cols).to_csv(args.out_all, sep="\t", index=False)
        pd.DataFrame(columns=cols).to_csv(args.out_sig, sep="\t", index=False)
        return

    # Stack all species-specific gene p-values
    stacks = []
    for fp in files:
        trait = os.path.basename(fp).replace(".genes.genes.out", "")
        df = read_geneout_one(fp)
        if df is None or df.empty:
            continue
        df["trait"] = trait
        df["genus"] = trait_to_genus(trait)
        stacks.append(df[["GENE", "P", "trait", "genus"]])

    if not stacks:
        warn("All files failed/empty; writing empty outputs.")
        cols = ["genus", "GENE_ID", "GENE", "P", "q", "Significant", "n_traits", "min_p"]
        pd.DataFrame(columns=cols).to_csv(args.out_all, sep="\t", index=False)
        pd.DataFrame(columns=cols).to_csv(args.out_sig, sep="\t", index=False)
        return

    D = pd.concat(stacks, ignore_index=True)
    D = D.rename(columns={"GENE": "GENE_ID"})
    D = D.merge(gmap, how="left", on="GENE_ID")  # add display name

    # 1) Combine across traits within (genus, GENE_ID) using Simes
    rows = []
    for (gn, gid), sub in D.groupby(["genus", "GENE_ID"], sort=False):
        pvals = sub["P"].astype(float).values
        pvals = pvals[(pvals > 0) & (pvals <= 1.0)]
        if pvals.size == 0:
            continue
        p_comb = simes_p(pvals)
        rows.append({
            "genus": gn,
            "GENE_ID": gid,
            "GENE": (sub["GENE"].dropna().iloc[0] if "GENE" in sub.columns and sub["GENE"].notna().any() else gid),
            "P": float(p_comb),
            "n_traits": int(pvals.size),
            "min_p": float(np.min(pvals))
        })

    if not rows:
        warn("After combining, no rows left; writing empty outputs.")
        cols = ["genus", "GENE_ID", "GENE", "P", "q", "Significant", "n_traits", "min_p"]
        pd.DataFrame(columns=cols).to_csv(args.out_all, sep="\t", index=False)
        pd.DataFrame(columns=cols).to_csv(args.out_sig, sep="\t", index=False)
        return

    C = pd.DataFrame(rows)

    # 2) BH FDR within each genus, on the combined p’s
    out_blocks = []
    for gn, sub in C.groupby("genus", sort=False):
        _, qv = fdrcorrection(sub["P"].values, alpha=args.alpha)
        tmp = sub.copy()
        tmp["q"] = qv
        tmp["Significant"] = tmp["q"] < args.alpha
        info(f"Genus {gn}: n={len(tmp)}, sig={int(tmp['Significant'].sum())} @ q<{args.alpha}")
        out_blocks.append(tmp)

    R = pd.concat(out_blocks, ignore_index=True)
    R = R[["genus", "GENE_ID", "GENE", "P", "q", "Significant", "n_traits", "min_p"]]

    R.to_csv(args.out_all, sep="\t", index=False)
    R[R["Significant"]].to_csv(args.out_sig, sep="\t", index=False)
    info(f"Wrote {args.out_all} (n={len(R)}) and {args.out_sig} (n={int(R['Significant'].sum())})")
    info("=== Collector (genus-combined) done ===")

if __name__ == "__main__":
    main()
