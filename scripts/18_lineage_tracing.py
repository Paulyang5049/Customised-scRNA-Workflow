"""
18 - Lineage Tracing
=====================
Reconstruct phylogenetic trees from CRISPR barcoding data using
Cassiopeia. Fits maximum-parsimony or neighbor-joining trees and
plots the resulting lineage hierarchy annotated by cell metadata.

"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)

    try:
        import cassiopeia as cas
    except ImportError:
        print("  ERROR: cassiopeia not installed. Install via: pip install cassiopeia-lineage")
        return

    # Expect character matrix in obsm or X
    if "character_matrix" in adata.obsm:
        character_matrix = pd.DataFrame(
            adata.obsm["character_matrix"],
            index=adata.obs_names,
        )
    else:
        character_matrix = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.obs_names,
            columns=adata.var_names,
        )

    # Replace NaN / -1 with missing state
    character_matrix = character_matrix.fillna(-1).astype(int)

    # Build CassiopeiaTree
    tree = cas.data.CassiopeiaTree(
        character_matrix=character_matrix,
        cell_meta=adata.obs,
    )

    # Reconstruct phylogeny
    if args.solver == "nj":
        solver = cas.solver.NeighborJoiningSolver(
            dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance,
            add_root=True,
        )
    else:
        solver = cas.solver.VanillaGreedySolver()

    solver.solve(tree)
    print(f"  Tree reconstructed with {tree.n_cell} leaves")

    # Annotate tree with cell metadata
    color_col = args.color_col if args.color_col in adata.obs.columns else None
    if color_col:
        fig, ax = plt.subplots(figsize=(10, 10))
        cas.pl.plot_matplotlib(tree, meta_data=[color_col], ax=ax)
        plt.savefig(os.path.join(args.figdir, "18_lineage_tree.png"), bbox_inches="tight", dpi=150)
        plt.close()
    else:
        fig, ax = plt.subplots(figsize=(10, 10))
        cas.pl.plot_matplotlib(tree, ax=ax)
        plt.savefig(os.path.join(args.figdir, "18_lineage_tree.png"), bbox_inches="tight", dpi=150)
        plt.close()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="18 - Lineage Tracing (Cassiopeia)")
    parser.add_argument("--input", required=True, help="Path to .h5ad with character matrix")
    parser.add_argument("--output", default="results/18_lineage.h5ad")
    parser.add_argument("--solver", choices=["nj", "greedy"], default="nj",
                        help="Tree solver: 'nj' (Neighbor Joining) or 'greedy'")
    parser.add_argument("--color_col", default="cell_type",
                        help="obs column for coloring tree leaves")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
