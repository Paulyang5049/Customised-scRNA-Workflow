"""
01 - Quality Control
====================
Filter low-quality cells using MAD-based outlier detection on QC metrics
(total counts, gene counts, mitochondrial %). Removes dead cells, empty
droplets, and doublets.


"""

import argparse
import os

import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation


def is_outlier(adata, metric: str, nmads: int):
    """Flag cells that deviate more than `nmads` MADs from the median."""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Annotate gene groups for QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    # Save QC plots
    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir
    import matplotlib.pyplot as plt

    sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    plt.savefig(os.path.join(args.figdir, "01_total_counts_hist.png"), bbox_inches="tight")
    plt.close()
    sc.pl.violin(adata, "pct_counts_mt", save="_01_pct_mt.png", show=False)
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
                  save="_01_counts_vs_genes.png", show=False)

    # MAD-based outlier detection
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > args.mt_threshold
    )
    print(f"  Outlier cells: {adata.obs['outlier'].sum()}")
    print(f"  MT outlier cells: {adata.obs['mt_outlier'].sum()}")

    # Filter
    n_before = adata.n_obs
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"  Cells after filtering: {adata.n_obs} (removed {n_before - adata.n_obs})")

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"  Genes after filtering: {adata.n_vars}")

    # Store raw counts for downstream steps
    adata.layers["counts"] = adata.X.copy()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="01 - Quality Control")
    parser.add_argument("--input", required=True, help="Path to raw .h5ad file")
    parser.add_argument("--output", default="results/01_qc.h5ad")
    parser.add_argument("--mt_threshold", type=float, default=8.0,
                        help="Max mitochondrial pct (default: 8)")
    parser.add_argument("--min_cells", type=int, default=20,
                        help="Min cells expressing a gene to keep it (default: 20)")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
