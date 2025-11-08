#!/usr/bin/env python3
import argparse, os, re, glob, gzip, sys
import numpy as np, pandas as pd

def trait_to_genus(t):
    t = re.sub(r"_assoc$", "", t)
    m = re.search(r"__([^_]+)", t)
    return m.group(1) if m else t.split("_")[0]

def read_assoc_min(path):
    """Read a GENESIS single-variant .tsv(.gz); return minimal standardized columns.
       Robust to empty files and mislabelled compression."""
    import gzip
    def _read(fh):
        return pd.read_csv(fh, sep="\t")
    # try gzip first if .gz; fallback to plain text on failure
    try:
        if path.endswith(".gz"):
            with gzip.open(path, "rt") as fh:
                df = _read(fh)
        else:
            with open(path, "rt") as fh:
                df = _read(fh)
    except Exception:
        try:
            with open(path, "rt") as fh:
                df = _read(fh)
        except Exception:
            return pd.DataFrame()
    if df is None or df.shape[0] == 0:
        return pd.DataFrame()
    lc = {c.lower(): c for c in df.columns}
    out = pd.DataFrame()
    for a in ("variant.id","vid","varid"):
        if a in lc: out["variant.id"] = df[lc[a]]; break
    for a in ("chr","chrom","chromosome"):
        if a in lc: out["chr"] = pd.to_numeric(df[lc[a]], errors="coerce"); break
    for a in ("pos","position","bp"):
        if a in lc: out["bp"]  = pd.to_numeric(df[lc[a]], errors="coerce"); break
    for a in ("score.pval","pval","p","p.value","pvalue"):
        if a in lc: out["pval"] = pd.to_numeric(df[lc[a]], errors="coerce"); break
    effect = None
    if "est" in lc: effect = df[lc["est"]]
    elif "beta" in lc: effect = df[lc["beta"]]
    elif "score.stat" in lc and "score.se" in lc:
        num = pd.to_numeric(df[lc["score.stat"]], errors="coerce")
        den = pd.to_numeric(df[lc["score.se"]], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            effect = num/den
    if effect is not None:
        out["effect"] = pd.to_numeric(effect, errors="coerce")
    return out


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc_dir", required=True)
    ap.add_argument("--magma_dir", required=True)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genera", default=None)
    return ap.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # MAGMA sidecars
    gene_loc_path = os.path.join(args.magma_dir, "magma.genes.loc")
    snp_loc_path  = os.path.join(args.magma_dir, "magma.snps.loc")
    sig_path = os.path.join(args.magma_dir, "magma_sig_species.tsv")
    if not os.path.exists(sig_path):
        raise SystemExit("magma_sig_species.tsv not found in --magma_dir")

    # gene windows (CHR, START, END for GENE_ID)
    gl = pd.read_csv(gene_loc_path, sep=r"\s+", header=None, engine="python").iloc[:, :4]
    gl.columns = ["GENE_ID","CHR","START","END"]
    gl["CHR"]   = pd.to_numeric(gl["CHR"], errors="coerce")
    gl["START"] = pd.to_numeric(gl["START"], errors="coerce")
    gl["END"]   = pd.to_numeric(gl["END"],   errors="coerce")

    # snp coords (SNP, CHR, BP) → used only to count / sanity check
    snp_loc = pd.read_csv(snp_loc_path, sep="\t", header=None, names=["SNP","CHR","BP"])
    snp_loc["CHR"] = pd.to_numeric(snp_loc["CHR"], errors="coerce")
    snp_loc["BP"]  = pd.to_numeric(snp_loc["BP"],  errors="coerce")

    # permissive (genus,gene,species) table
    sig = pd.read_csv(sig_path, sep="\t")
    if "GENE_ID" not in sig.columns:
        if "gene" in sig.columns: sig = sig.rename(columns={"gene":"GENE_ID"})
        elif "GENE" in sig.columns: sig = sig.rename(columns={"GENE":"GENE_ID"})
        else: raise SystemExit("sig table has no GENE_ID/gene/GENE column")
    sig["GENE_ID"] = sig["GENE_ID"].astype(str)
    sig["genus"]   = sig["genus"].astype(str)
    if args.genera:
        keep = set([x for x in args.genera.split(",") if x])
        sig = sig[sig["genus"].isin(keep)].copy()

    # enumerate folds
    folds = sorted([d for d in glob.glob(os.path.join(args.assoc_dir, "fold*")) if os.path.isdir(d)])
    if not folds: raise SystemExit("[crossfit] No foldXX directories under assoc_dir.")
    def traits_for(subdir):
        paths = sorted(glob.glob(os.path.join(subdir, "*.tsv.gz")))
        if not paths: return pd.DataFrame(columns=["trait","file","genus"])
        T = [(os.path.basename(p).replace(".tsv.gz",""), p) for p in paths]
        df = pd.DataFrame(T, columns=["trait","file"])
        df["genus"] = df["trait"].map(trait_to_genus)
        return df
    fold_info=[]
    for f in folds:
        tr_train = traits_for(os.path.join(f, "train"))
        tr_test  = traits_for(os.path.join(f, "test"))
        if tr_train.empty or tr_test.empty: continue
        fold_info.append((f, tr_train, tr_test))
    if not fold_info: raise SystemExit("[crossfit] No train/test files found under foldXX.")
    print(f"[crossfit] Folds detected: {len(fold_info)}", flush=True)

    sent_rows, beta_rows = [], []

    for (fold_path, train_map, test_map) in fold_info:
        fold_name = os.path.basename(fold_path)

        # Build a variant.id lookup once per fold from *one* non-empty train file
        varmap = None
        for p in train_map["file"]:
            df = read_assoc_min(p)
            if not df.empty and "variant.id" in df.columns and "chr" in df.columns and "bp" in df.columns:
                varmap = df[["variant.id","chr","bp"]].dropna().drop_duplicates()
                break
        if varmap is None or varmap.empty:
            print(f"[crossfit][{fold_name}] WARNING: could not build varmap (no non-empty train files).", flush=True)
            continue

        # pre-cache train/test minimal dfs by trait (to avoid re-reading)
        train_cache = {}
        for tr, p in zip(train_map["trait"], train_map["file"]):
            a = read_assoc_min(p)
            if not a.empty: train_cache[tr] = a
        test_cache = {}
        for tr, p in zip(test_map["trait"], test_map["file"]):
            a = read_assoc_min(p)
            if not a.empty: test_cache[tr] = a

        print(f"[crossfit][{fold_name}] train_nonempty={len(train_cache)} test_nonempty={len(test_cache)}", flush=True)
        if not train_cache or not test_cache:
            continue

        for gn in sorted(sig["genus"].unique()):
            sig_g = sig[sig["genus"]==gn]
            train_traits = [t for t in train_map["trait"] if t in train_cache and trait_to_genus(t)==gn]
            test_traits  = [t for t in test_map["trait"]  if t in test_cache  and trait_to_genus(t)==gn]
            if not train_traits or not test_traits:
                continue

            for gid in sig_g["GENE_ID"].unique():
                W = gl[gl["GENE_ID"]==gid]
                if W.empty: continue
                w = W.iloc[0]
                # gene window → candidate variant.ids (via varmap join on chr/bp)
                gene_vids = varmap[(varmap["chr"]==w["CHR"]) & (varmap["bp"]>=w["START"]) & (varmap["bp"]<=w["END"])]
                if gene_vids.empty: continue
                vids = set(gene_vids["variant.id"].astype(int))

                # TRAIN: pick sentinel by lowest p across species within vids
                best = (None, np.inf)  # (vid, p)
                for tr in train_traits:
                    a = train_cache[tr]
                    a2 = a[a["variant.id"].isin(vids)]
                    if a2.empty or "pval" not in a2.columns: continue
                    idx = int(np.nanargmin(a2["pval"].values))
                    p = a2["pval"].values[idx]
                    vid = int(a2["variant.id"].values[idx])
                    if np.isfinite(p) and p < best[1]:
                        best = (vid, p)
                if best[0] is None:
                    continue
                sentinel_vid, best_p = best
                sent_rows.append({"fold": fold_name, "genus": gn, "GENE_ID": gid,
                                  "sentinel_variant.id": sentinel_vid, "train_minP": best_p})

                # TEST: collect held-out effects at that variant.id across all species in genus
                for tr in test_traits:
                    a = test_cache[tr]
                    a2 = a[a["variant.id"]==sentinel_vid]
                    if a2.empty or "effect" not in a2.columns: continue
                    eff = pd.to_numeric(a2["effect"], errors="coerce")
                    if eff.size==0 or not np.isfinite(eff.iloc[0]): continue
                    beta_rows.append({
                        "fold": fold_name, "genus": gn, "GENE_ID": gid, "trait": tr,
                        "variant.id": sentinel_vid, "beta_test": float(eff.iloc[0])
                    })

    sent = pd.DataFrame(sent_rows)
    betas = pd.DataFrame(beta_rows)

    if not betas.empty:
        cf = betas.groupby(["genus","GENE_ID","trait"], as_index=False).agg(
            beta_cf=("beta_test","mean"), n_folds_used=("beta_test","count"))
    else:
        cf = pd.DataFrame(columns=["genus","GENE_ID","trait","beta_cf","n_folds_used"])

    sent.to_csv(os.path.join(args.outdir, "sentinels_crossfit.tsv"), sep="\t", index=False)
    cf.to_csv(os.path.join(args.outdir, "beta_cf_by_species.tsv"), sep="\t", index=False)
    print(f"[crossfit] Wrote sentinels_crossfit.tsv (n={len(sent)}) and beta_cf_by_species.tsv (n={len(cf)})", flush=True)

if __name__ == "__main__":
    main()
