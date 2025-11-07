#!/usr/bin/env python3
"""
micnet.py – auto-tune UMAP-HDBSCAN parameters, then export clusters.
Stops as soon as DBCV and ARI reach their targets, re-uses that same
embedding to save Cluster_UMAP_HDBSCAN.csv and the 2-D PDF plot.
"""

# ── basic imports ────────────────────────────────────────────────────
import os, sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
import umap
import matplotlib.pyplot as plt, seaborn as sns

from utils import filter_otus
from umap_hdbscan import (EmbeddingOutput,
                             intrinsic_validity,
                             validity_index,
                             bootstrap_ari)

# ── tiny helpers ─────────────────────────────────────────────────────
def kind_file(x): return Path(x).suffix.lower() == ".txt"
def str2bool(v): return v.lower() in ("yes", "true", "t", "1")

# ── CLI parsing ──────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser("MicNet auto-tune driver")
    p.add_argument("--tool", required=True)
    p.add_argument("--file_path", required=True)
    p.add_argument("--output_dir", default="./")

    # base hyper-parameters
    p.add_argument("--n_neighbors",    type=int,   default=30)
    p.add_argument("--min_dist",       type=float, default=0.05)
    p.add_argument("--n_components",   type=int,   default=20)
    p.add_argument("--metric_umap",    default="cosine")
    p.add_argument("--metric_hdb",     default="euclidean")
    p.add_argument("--min_cluster_size", type=int, default=15)
    p.add_argument("--target_dbcv",       type=float, default=0.50)
    p.add_argument("--target_ari",       type=float, default=0.60)
    p.add_argument("--min_sample",       type=int, default=10)
    p.add_argument("--cluster_epsilon",  type=float, default=0.01)

    p.add_argument("--taxa",   type=str2bool, default=False)
    p.add_argument("--abundance_filter", type=str2bool, default=False)
    return p.parse_args()

# ── export results without extra bootstrap ───────────────────────────
def export_results(emb, outliers, labels, X, params,
                   file_in, out_dir, taxa_flag, taxa_labels,
                   dbcv, ari,target_dbcv, target_ari):
    sil, db, p_sil, p_db = intrinsic_validity(X.values, labels)

    verdict = ("VALID" if (p_sil < 0.01 and p_db < 0.01 and dbcv >= target_dbcv and ari >= target_ari)
               else "UNSTABLE")

    print("\n────────── FINAL CLUSTERING ──────────")
    print(f"Best params : nn={params['n_neighbors']}  "
          f"mcs={params['min_cluster_size']}  "
          f"eps={params['cluster_epsilon']:.2f}")
    print(f"[UMAP-HDBSCAN] Sil={sil:.3f} (P={p_sil:.3g})  "
          f"DB={db:.3f} (P={p_db:.3g})  DBCV={dbcv:.3f}")
    print(f"Bootstrap ARI (20× 90 %):  {ari:.2f}")
    print(f"Verdict : {verdict}")
    print("──────────────────────────────────────\n")

    res = pd.DataFrame({'Outliers': outliers, 'Cluster': labels})
    if taxa_flag: res['Taxa'] = taxa_labels
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "Cluster_UMAP_HDBSCAN.csv")

    umap2d = umap.UMAP(n_neighbors=params["n_neighbors"],
                       min_dist=params["min_dist"],
                       n_components=2,
                       metric=params["metric_umap"],
                       random_state=42).fit_transform(emb)
    palette = sns.color_palette("tab20", len(set(labels)))
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=umap2d[:, 0], y=umap2d[:, 1],
                    hue=labels, palette=palette, s=25, linewidth=0)
    plt.title("UMAP-HDBSCAN Clusters")
    plt.tight_layout()
    plt.savefig(out_dir / "UMAP_HDBSCAN_Clusters.pdf")
    plt.close()

# ── auto-tune loop ───────────────────────────────────────────────────
def autotune_export(base_params, X, file_in, taxa_flag,
                    taxa_labels, out_dir,
                    target_dbcv, target_ari):

    for nn in [30, 40, 50]:
        for mc in [15, 18, 22]:
            for ce in [0.01, 0.03, 0.05]:
                trial = base_params | {"n_neighbors": nn,
                                       "min_cluster_size": mc,
                                       "cluster_epsilon": ce}
                emb, outliers, lab = EmbeddingOutput(**trial).fit(X)

                dbcv = validity_index(np.ascontiguousarray(emb, dtype=np.float64),
                                      lab.astype(int), metric='euclidean')
                if dbcv < target_dbcv:
                    print(f"Grid nn={nn:2d} mcs={mc:2d} eps={ce:.2f} "
                          f"→ DBCV={dbcv:.3f}  (skip ARI)")
                    continue

                ari, _ = bootstrap_ari(X.values, trial, n_iter=20, frac=0.9)
                print(f"Grid nn={nn:2d} mcs={mc:2d} eps={ce:.2f} "
                      f"→ DBCV={dbcv:.3f}  ARI={ari:.2f}")

                if ari >= target_ari:
                    print("\n★★ TARGET MET – exporting ★★")
                    export_results(emb, outliers, lab, X, trial,
                                   file_in, out_dir,
                                   taxa_flag, taxa_labels,
                                   dbcv, ari, target_dbcv, target_ari)
                    return
    print("\nNo parameter met the thresholds.")

# ── main ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    a = parse_args()
    if a.tool != "umap_hdbscan":
        sys.exit("This file only supports --tool umap_hdbscan")

    # load matrix
    df = pd.read_table(a.file_path) if kind_file(a.file_path) else pd.read_csv(a.file_path)
    if a.taxa:
        X = df.iloc[:, 2:].astype(float)
        idx, X = filter_otus(X, a.abundance_filter)
        taxa_labels = (df.iloc[idx, 1].str.split(';').str.get(0) + '-' +
                       df.iloc[idx, 1].str.split(';').str.get(1) + '-' +
                       df.iloc[idx, 1].str.split(';').str.get(5))
    else:
        X = df.iloc[:, 1:].astype(float)
        idx, X = filter_otus(X, a.abundance_filter)
        taxa_labels = None

    base = dict(n_neighbors=a.n_neighbors,
                min_dist=a.min_dist,
                n_components=a.n_components,
                metric_umap=a.metric_umap,
                metric_hdb=a.metric_hdb,
                min_cluster_size=a.min_cluster_size,
                min_sample=a.min_sample,
                target_dbcv=a.target_dbcv,
                target_ari=a.target_ari,
                cluster_epsilon=a.cluster_epsilon,
                random_state=42)

    autotune_export(base, X, a.file_path, a.taxa, taxa_labels, a.output_dir,target_dbcv=a.target_dbcv, target_ari=a.target_ari)
