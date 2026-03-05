"""
12 - RNA Velocity
=================
Estimate RNA velocity using scVelo (deterministic and dynamical modes).
Requires spliced/unspliced count layers in the AnnData object.


"""

import argparse
import os
import warnings

import scanpy as sc

warnings.filterwarnings("ignore", category=DeprecationWarning)


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    import scvelo as scv

    os.makedirs(args.figdir, exist_ok=True)
    scv.settings.figdir = args.figdir
    scv.settings.set_figure_params("scvelo")

    # Check that spliced/unspliced layers exist
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        print("  ERROR: spliced/unspliced layers not found. RNA velocity requires loom-derived data.")
        return

    # Filter and normalize
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

    # --- Deterministic mode ---
    scv.tl.velocity(adata, mode="deterministic")
    scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)

    cluster_key = args.cluster_key if args.cluster_key in adata.obs.columns else None
    scv.pl.velocity_embedding_stream(adata, basis="umap", color=cluster_key,
                                     save="12_velocity_deterministic.png", show=False)
    print("  Deterministic velocity: done")

    # --- Dynamical mode ---
    if args.run_dynamical:
        scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
        scv.pl.velocity_embedding_stream(adata, basis="umap", color=cluster_key,
                                         save="12_velocity_dynamical.png", show=False)

        # Top likelihood genes
        top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:5]
        scv.pl.scatter(adata, basis=top_genes, color=cluster_key, frameon=False,
                       save="12_top_likelihood_genes.png", show=False)
        print("  Dynamical velocity: done")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="12 - RNA Velocity")
    parser.add_argument("--input", required=True, help="Path to .h5ad with spliced/unspliced layers")
    parser.add_argument("--output", default="results/12_velocity.h5ad")
    parser.add_argument("--cluster_key", default="clusters")
    parser.add_argument("--n_jobs", type=int, default=4)
    parser.add_argument("--run_dynamical", action="store_true",
                        help="Also run the slower dynamical model")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
