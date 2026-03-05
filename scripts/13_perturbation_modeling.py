"""
13 - Perturbation Modeling
===========================
Rank cell types by their responsiveness to experimental perturbation
using Augur (via pertpy). Augur trains classifiers per cell type to
distinguish control vs perturbed cells — higher AUC = more responsive.


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)

    import pertpy as pt

    # --- Augur ---
    ag = pt.tl.Augur("random_forest_classifier")

    loaded, results = ag.load(adata, label_col=args.condition_key,
                              cell_type_col=args.cell_type_key)

    adata_augur, results = ag.predict(
        loaded,
        subsample_size=args.subsample_size,
        n_threads=args.n_threads,
        select_variance_features=True,
        span=0.75,
    )

    # Summary
    summary = results["summary_metrics"]
    print(f"  Augur AUC per cell type:\n{summary}")
    summary.to_csv(os.path.join(args.figdir, "13_augur_auc.csv"))

    # Lollipop plot
    ag.plot_lollipop(results)
    plt.savefig(os.path.join(args.figdir, "13_augur_lollipop.png"), bbox_inches="tight")
    plt.close()

    # UMAP colored by Augur score
    if "X_umap" in adata.obsm:
        ag.plot_umap(results, adata=adata, embedding="X_umap")
        plt.savefig(os.path.join(args.figdir, "13_augur_umap.png"), bbox_inches="tight")
        plt.close()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="13 - Perturbation Modeling (Augur)")
    parser.add_argument("--input", required=True, help="Path to .h5ad")
    parser.add_argument("--output", default="results/13_perturbation.h5ad")
    parser.add_argument("--condition_key", default="label")
    parser.add_argument("--cell_type_key", default="cell_type")
    parser.add_argument("--subsample_size", type=int, default=20)
    parser.add_argument("--n_threads", type=int, default=4)
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
