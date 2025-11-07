#!/usr/bin/env python3
"""
umap_hdbscan_v2.py
──────────────────────────────────────────────────────────────
UMAP-HDBSCAN clustering, intrinsic validity, DBCV, bootstrap
Returns: embedding, outlier_mask, labels
"""

import argparse, glob, os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import umap, hdbscan
from tqdm import tqdm
from sklearn.metrics import (silhouette_score, davies_bouldin_score,
                             adjusted_rand_score)
from sklearn.utils import check_random_state
import matplotlib.pyplot as plt, seaborn as sns
from hdbscan.validity import validity_index

warnings.filterwarnings("ignore",
                        message="n_jobs value -1 overridden to 1 by setting random_state")

# ───────────────── Embedding + clustering ──────────────────
class EmbeddingOutput:
    def __init__(self, **kw):
        self.kw = kw
        self.random_state   = kw.get("random_state", 42)
        self.quantile_limit = kw.get("quantile_limit", 0.90)

    def fit(self, X):
        emb = umap.UMAP(
            n_neighbors   = self.kw["n_neighbors"],
            min_dist      = self.kw["min_dist"],
            n_components  = self.kw["n_components"],
            metric        = self.kw["metric_umap"],
            output_metric = "euclidean",
            random_state  = self.random_state
        ).fit_transform(X)

        cl = hdbscan.HDBSCAN(
            min_cluster_size         = self.kw["min_cluster_size"],
            min_samples              = self.kw["min_sample"],
            metric                   = self.kw["metric_hdb"],
            cluster_selection_epsilon= self.kw["cluster_epsilon"],
            core_dist_n_jobs         = 1
        ).fit(emb)

        thresh   = np.quantile(cl.outlier_scores_, self.quantile_limit)
        outliers = (cl.outlier_scores_ > thresh).astype(int)
        return emb, outliers, cl.labels_


# ─────────────────── helper functions ──────────────────────
def intrinsic_validity(X, labels, n_perm=1000, rng=42):
    keep   = labels != -1
    X_, lb = X[keep], labels[keep]
    sil    = silhouette_score(X_, lb)
    db     = davies_bouldin_score(X_, lb)

    r = check_random_state(rng)
    sil_p = [silhouette_score(X_, r.permutation(lb)) for _ in range(n_perm)]
    db_p  = [davies_bouldin_score(X_, r.permutation(lb)) for _ in range(n_perm)]
    p_sil = (sum(s >= sil for s in sil_p) + 1) / (n_perm + 1)
    p_db  = (sum(d <= db  for d in db_p)  + 1) / (n_perm + 1)
    return sil, db, p_sil, p_db


def bootstrap_ari(X, params, n_iter=30, frac=0.8, rng=42):
    r = check_random_state(rng)
    labs = []
    for _ in tqdm(range(n_iter), desc="  ↳ bootstraps", leave=False):
        idx = r.choice(X.shape[0], int(frac * X.shape[0]), replace=False)
        _, _, lab = EmbeddingOutput(**params).fit(X[idx])
        full = np.full(X.shape[0], -1, dtype=int)
        full[idx] = lab
        labs.append(full)

    ari = [adjusted_rand_score(labs[i], labs[j])
           for i in range(len(labs)) for j in range(i + 1, len(labs))]
    return np.mean(ari), np.std(ari)


# ───────────────────── CLI parsing ─────────────────────────
def cli():
    p = argparse.ArgumentParser(description="UMAP-HDBSCAN validator")
    p.add_argument("--tool", required=True)
    p.add_argument("--file_path", required=True)
    p.add_argument("--output_dir", default="./")
    # UMAP
    p.add_argument("--n_neighbors", type=int, default=30)
    p.add_argument("--min_dist",    type=float, default=0.05)
    p.add_argument("--n_components",type=int, default=20)
    p.add_argument("--metric_umap", default="cosine")
    # HDBSCAN
    p.add_argument("--metric_hdb",     default="euclidean")
    p.add_argument("--min_cluster_size", type=int, default=12)
    p.add_argument("--min_sample",       type=int, default=5)
    p.add_argument("--cluster_epsilon",  type=float, default=0.03)
    # misc
    p.add_argument("--quantile_limit",   type=float, default=0.90)
    p.add_argument("--random_state",     type=int, default=42)
    return p.parse_args()


# ─────────────────────────── main ──────────────────────────
def main():
    args = cli()
    if args.tool != "umap_hdbscan":
        sys.exit("Only --tool umap_hdbscan supported")

    mats = glob.glob(os.path.join(args.file_path, "*"))
    if not mats:
        sys.exit(f"No files in {args.file_path}")
    print(f"[1/5] Found {len(mats)} file(s)")

    params = vars(args)
    outdir = Path(args.output_dir); outdir.mkdir(parents=True, exist_ok=True)

    for f in mats:
        print(f"\n[2/5] Loading → {os.path.basename(f)}")
        X = (pd.read_pickle(f) if f.endswith(".pkl")
             else pd.read_csv(f, index_col=0)).to_numpy(float)

        print("[3/5] Clustering …")
        emb, outliers, labels = EmbeddingOutput(**params).fit(X)

        umap2d = umap.UMAP(n_neighbors=args.n_neighbors, min_dist=args.min_dist,
                           n_components=2, metric=args.metric_umap,
                           random_state=args.random_state).fit_transform(emb)
        palette = sns.color_palette("tab10", len(set(labels)))
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=umap2d[:, 0], y=umap2d[:, 1],
                        hue=labels, palette=palette, s=25, linewidth=0)
        plt.title("UMAP-HDBSCAN clusters")
        plt.tight_layout()
        plt.savefig(outdir / f"{Path(f).stem}_clusters.pdf")
        plt.close()

        print("[4/5] Intrinsic validity (1 000 perms) …")
        sil, db, p_sil, p_db = intrinsic_validity(X, labels)
        print(f"    Silhouette = {sil:.3f} (P={p_sil:.3g})")
        print(f"    DB-index   = {db:.3f} (P={p_db:.3g})")
        try:
            dbcv = validity_index(np.ascontiguousarray(emb, dtype=np.float64),
                                  labels.astype(np.intp), metric='euclidean')
            print(f"    DBCV       = {dbcv:.3f}")
        except Exception as e:
            print(f"    DBCV could not be computed: {e}")

        print("[5/5] Bootstrap stability (30 resamples) …")
        mean_ari, sd_ari = bootstrap_ari(X, params, n_iter=30)
        verdict = ("VALID" if (p_sil < 0.01 and p_db < 0.01 and mean_ari >= 0.70)
                   else "SUGGESTIVE" if mean_ari >= 0.60 else "UNSTABLE")
        print("────────────────────────────────────────────")
        print(f" Verdict: {verdict}")
        print(f" Mean ARI = {mean_ari:.2f} ± {sd_ari:.2f}")
        print("────────────────────────────────────────────")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    main()
