"""
06 - Clustering
===============
Run Leiden graph-based clustering at multiple resolutions and visualize
on UMAP.

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
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

    # Compute neighbors if not already present
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_pcs=args.n_pcs)
        sc.tl.umap(adata)

    # Run Leiden at multiple resolutions
    resolutions = [float(r) for r in args.resolutions.split(",")]
    keys = []
    for res in resolutions:
        key = f"leiden_res{res}"
        sc.tl.leiden(adata, key_added=key, resolution=res, flavor="igraph", n_iterations=2)
        n_clusters = adata.obs[key].nunique()
        print(f"  Resolution {res}: {n_clusters} clusters")
        keys.append(key)

    # Plot
    sc.pl.umap(adata, color=keys, legend_loc="on data",
               save="_06_clustering.png", show=False)

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="06 - Clustering")
    parser.add_argument("--input", default="results/04_dimred.h5ad")
    parser.add_argument("--output", default="results/06_clustered.h5ad")
    parser.add_argument("--resolutions", default="0.25,0.5,1.0",
                        help="Comma-separated Leiden resolutions")
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
