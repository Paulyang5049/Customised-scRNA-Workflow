"""
03 - Feature Selection
======================
Identify highly variable genes (HVGs) using Scanpy's cell_ranger flavor
and optionally deviance-based selection via scry (requires R).


"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # --- Deviance-based feature selection (via R/scry) ---
    try:
        import rpy2.robjects as ro
        import rpy2.robjects.packages as rpackages
        from rpy2.robjects import default_converter, numpy2ri, pandas2ri, r
        from rpy2.robjects.conversion import localconverter
        from scipy.sparse import issparse

        X_sparse = adata.X.T.tocoo()
        Matrix = rpackages.importr("Matrix")

        with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
            ro.globalenv["obs"] = adata.obs
            ro.globalenv["var"] = adata.var

        i, j = X_sparse.row, X_sparse.col
        x = X_sparse.data
        ro.globalenv["i"] = ro.IntVector((i + 1).tolist())
        ro.globalenv["j"] = ro.IntVector((j + 1).tolist())
        ro.globalenv["x"] = ro.FloatVector(x.tolist())
        r(f"X <- sparseMatrix(i = i, j = j, x = x, dims = c({X_sparse.shape[0]}, {X_sparse.shape[1]}))")

        ro.r("""
            library(scry)
            library(SingleCellExperiment)
            sce <- SingleCellExperiment(assays = list(X = X), colData = obs, rowData = var)
            sce <- devianceFeatureSelection(sce, assay = "X")
        """)

        with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
            binomial_deviance = ro.r("rowData(sce)$binomial_deviance")

        idx = binomial_deviance.argsort()[-args.n_top_genes:]
        mask = np.zeros(adata.var_names.shape, dtype=bool)
        mask[idx] = True
        adata.var["highly_deviant"] = mask
        adata.var["binomial_deviance"] = binomial_deviance
        print(f"  Deviance-based selection: {mask.sum()} HVGs")
    except Exception as e:
        print(f"  Deviance-based selection skipped (R not available): {e}")

    # --- Scanpy HVG selection ---
    norm_layer = "scran_normalization" if "scran_normalization" in adata.layers else None
    sc.pp.highly_variable_genes(adata, layer=norm_layer)
    print(f"  Scanpy HVGs: {adata.var['highly_variable'].sum()}")

    # Plot
    ax = sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5)
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    plt.savefig(os.path.join(args.figdir, "03_hvg_selection.png"), bbox_inches="tight")
    plt.close()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="03 - Feature Selection")
    parser.add_argument("--input", default="results/02_normalized.h5ad")
    parser.add_argument("--output", default="results/03_feature_selected.h5ad")
    parser.add_argument("--n_top_genes", type=int, default=4000)
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
