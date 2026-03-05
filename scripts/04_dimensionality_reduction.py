"""
04 - Dimensionality Reduction
=============================
Run PCA on highly variable genes, then compute t-SNE and UMAP embeddings
for 2D visualization.


"""

import argparse
import os

import scanpy as sc


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # Use scran-normalized layer if available, else log1p
    if "scran_normalization" in adata.layers:
        adata.X = adata.layers["scran_normalization"]
    elif "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"]

    # Mark HVGs: prefer deviance-based if available
    if "highly_deviant" in adata.var.columns:
        adata.var["highly_variable"] = adata.var["highly_deviant"]

    # PCA
    sc.pp.pca(adata, svd_solver="arpack", n_comps=args.n_pcs, mask_var="highly_variable")
    sc.pl.pca_scatter(adata, color="total_counts", save="_04_pca.png", show=False)

    # t-SNE
    sc.tl.tsne(adata, use_rep="X_pca")
    sc.pl.tsne(adata, color="total_counts", save="_04_tsne.png", show=False)

    # UMAP
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="total_counts", save="_04_umap.png", show=False)

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="04 - Dimensionality Reduction")
    parser.add_argument("--input", default="results/03_feature_selected.h5ad")
    parser.add_argument("--output", default="results/04_dimred.h5ad")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCs")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
