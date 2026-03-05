"""
14 - Cell-Cell Communication
==============================
Infer receptor-ligand interactions using LIANA (LIgand-receptor ANAlysis).
LIANA wraps several CCC methods; we use CellPhoneDB-style scoring by default.


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)

    import liana as li

    # Overwrite X with normalized values if available
    if "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"].copy()

    # Ensure cell type column is clean
    if args.cell_type_key not in adata.obs.columns:
        print(f"  ERROR: {args.cell_type_key} not in adata.obs")
        return
    adata.obs[args.cell_type_key] = adata.obs[args.cell_type_key].astype("category")

    # Drop small groups
    min_cells = 10
    counts = adata.obs[args.cell_type_key].value_counts()
    keep = counts[counts >= min_cells].index
    adata = adata[adata.obs[args.cell_type_key].isin(keep)].copy()
    adata.obs[args.cell_type_key] = adata.obs[args.cell_type_key].cat.remove_unused_categories()

    # Run LIANA
    li.mt.rank_aggregate(
        adata,
        groupby=args.cell_type_key,
        n_perms=args.n_perms,
        use_raw=False,
        verbose=True,
    )

    # Export interactions table
    interactions = adata.uns["liana_res"]
    interactions.to_csv(os.path.join(args.figdir, "14_liana_interactions.csv"), index=False)
    print(f"  LIANA inferred {len(interactions)} interactions")

    # Dotplot of top interactions
    try:
        li.pl.dotplot(
            adata=adata,
            colour="lr_means",
            size="cellphone_pvals",
            inverse_size=True,
            source_labels=keep[:5].tolist(),
            target_labels=keep[:5].tolist(),
            top_n=20,
            orderby="lr_means",
            orderby_ascending=False,
        )
        plt.savefig(os.path.join(args.figdir, "14_liana_dotplot.png"), bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"  Dotplot skipped: {e}")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="14 - Cell-Cell Communication (LIANA)")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results/14_ccc.h5ad")
    parser.add_argument("--cell_type_key", default="cell_type")
    parser.add_argument("--n_perms", type=int, default=100, help="Number of permutations")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
