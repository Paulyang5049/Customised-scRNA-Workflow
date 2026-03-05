"""
02 - Normalization
==================
Apply multiple normalization strategies:
  1. Shifted logarithm (log1p after size-factor normalization)
  2. Scran pooling-based size factors (requires R/scran)
  3. Analytic Pearson residuals


"""

import argparse
import os

import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import csr_matrix, issparse


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # --- Method 1: Shifted logarithm ---
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    # --- Method 2: Scran pooling (pure Python fallback if R unavailable) ---
    try:
        import logging
        import rpy2.rinterface_lib.callbacks as rcb
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter

        rcb.logger.setLevel(logging.ERROR)

        # Preliminary clustering for deconvolution
        adata_pp = adata.copy()
        sc.pp.normalize_total(adata_pp)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp, n_comps=15)
        sc.pp.neighbors(adata_pp)
        sc.tl.leiden(adata_pp, key_added="groups", flavor="igraph", n_iterations=2, directed=False)

        data_mat = adata_pp.X.T
        if issparse(data_mat):
            if data_mat.nnz > 2**31 - 1:
                data_mat = data_mat.tocoo()
            else:
                data_mat = data_mat.tocsc()
            data_mat = data_mat.toarray()

        with localconverter(ro.default_converter + numpy2ri.converter):
            ro.globalenv["data_mat"] = data_mat
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["input_groups"] = adata_pp.obs["groups"]

        del adata_pp

        ro.r("""
            library(scran)
            library(BiocParallel)
            size_factors <- sizeFactors(
                computeSumFactors(
                    SingleCellExperiment(list(counts=data_mat)),
                    clusters = input_groups, min.mean = 0.1,
                    BPPARAM = MulticoreParam()
                )
            )
        """)
        size_factors = np.array(ro.globalenv["size_factors"])

        adata.obs["size_factors"] = size_factors
        scran = adata.X / adata.obs["size_factors"].values[:, None]
        adata.layers["scran_normalization"] = csr_matrix(np.log1p(scran))
        print("  Scran normalization: done (via R)")
    except Exception as e:
        print(f"  Scran normalization skipped (R not available): {e}")

    # --- Method 3: Analytic Pearson residuals ---
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

    # --- Comparison plots ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")
    plt.tight_layout()
    plt.savefig(os.path.join(args.figdir, "02_normalization_comparison.png"), bbox_inches="tight")
    plt.close()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="02 - Normalization")
    parser.add_argument("--input", default="results/01_qc.h5ad")
    parser.add_argument("--output", default="results/02_normalized.h5ad")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
