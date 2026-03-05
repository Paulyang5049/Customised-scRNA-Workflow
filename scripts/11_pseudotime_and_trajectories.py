"""
11 - Pseudotime & Trajectory Inference
=======================================
Compute diffusion pseudotime (DPT) via Scanpy. Compares with external
pseudotime estimates (e.g. Palantir) if available.


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

    # Preprocessing
    sc.pp.filter_genes(adata, min_counts=20)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)

    # Diffusion map
    sc.tl.diffmap(adata)

    # Set root cell: the cell with the lowest value on diffusion component 3
    root_ixs = adata.obsm["X_diffmap"][:, 3].argmin()
    adata.uns["iroot"] = root_ixs
    print(f"  Root cell index: {root_ixs}")

    # Diffusion pseudotime
    sc.tl.dpt(adata)

    # Visualizations
    color_cols = ["dpt_pseudotime"]
    if "palantir_pseudotime" in adata.obs.columns:
        color_cols.append("palantir_pseudotime")
    if args.cluster_key in adata.obs.columns:
        color_cols.append(args.cluster_key)

    basis = "tsne" if "X_tsne" in adata.obsm else "umap"
    sc.pl.scatter(adata, basis=basis, color=color_cols, color_map="gnuplot2",
                  save="_11_pseudotime.png", show=False)

    if args.cluster_key in adata.obs.columns:
        sc.pl.violin(adata, keys=["dpt_pseudotime"], groupby=args.cluster_key,
                     rotation=45, save="_11_pseudotime_violin.png", show=False)

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="11 - Pseudotime & Trajectories")
    parser.add_argument("--input", required=True, help="Path to .h5ad")
    parser.add_argument("--output", default="results/11_pseudotime.h5ad")
    parser.add_argument("--cluster_key", default="clusters")
    parser.add_argument("--n_pcs", type=int, default=10)
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
