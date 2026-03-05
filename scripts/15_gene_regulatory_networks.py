"""
15 - Gene Regulatory Networks (GRNs)
======================================
Infer gene regulatory networks using pySCENIC:
  1. Co-expression modules (GRNBoost2 / adjacencies)
  2. Regulon pruning via cis-regulatory motif enrichment (cisTarget)
  3. Regulon activity scoring per cell (AUCell)


"""

import argparse
import glob
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
        from pyscenic.utils import modules_from_adjacencies
        from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
        from pyscenic.prune import prune2df, df2regulons
        from pyscenic.aucell import aucell
        import loompy
    except ImportError:
        print("  ERROR: pySCENIC not installed. Install via: pip install pyscenic")
        return

    from arboreto.algo import grnboost2

    # Use raw counts
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    # Convert to DataFrame for arboreto
    sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.01))
    expr_mat = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names, columns=adata.var_names,
    )

    # Load TF list
    if not os.path.isfile(args.tf_list):
        print(f"  ERROR: TF list not found at {args.tf_list}")
        return
    tf_names = pd.read_csv(args.tf_list, header=None)[0].tolist()

    # Step 1: Co-expression modules
    print("  Step 1/3: GRNBoost2 co-expression inference...")
    adjacencies = grnboost2(expr_mat, tf_names=tf_names, verbose=True)
    modules = list(modules_from_adjacencies(adjacencies, expr_mat))
    print(f"    Found {len(modules)} modules")

    # Step 2: cisTarget motif enrichment
    print("  Step 2/3: cisTarget regulon pruning...")
    ranking_dbs = [RankingDatabase(fname=f, name=os.path.basename(f))
                   for f in glob.glob(os.path.join(args.db_dir, "*.feather"))]
    if not ranking_dbs:
        print(f"    No ranking databases found in {args.db_dir}")
        adata.write_h5ad(args.output)
        return

    df_motifs = prune2df(ranking_dbs, modules, args.motif_annotations)
    regulons = df2regulons(df_motifs)
    print(f"    Pruned to {len(regulons)} regulons")

    # Step 3: AUCell activity scoring
    print("  Step 3/3: AUCell regulon activity scoring...")
    auc_mtx = aucell(expr_mat, regulons)
    adata.obsm["X_aucell"] = auc_mtx.values

    # Save regulon names
    adata.uns["regulon_names"] = [r.name for r in regulons]

    # Heatmap of top regulon activities
    top_n = min(20, len(regulons))
    top_regulons = auc_mtx.var(axis=0).nlargest(top_n).index.tolist()

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.imshow(auc_mtx[top_regulons].T.values, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(top_regulons)))
    ax.set_yticklabels(top_regulons, fontsize=8)
    ax.set_xlabel("Cells")
    ax.set_title("Top regulon AUCell activities")
    plt.savefig(os.path.join(args.figdir, "15_grn_regulon_heatmap.png"), bbox_inches="tight")
    plt.close()

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="15 - Gene Regulatory Networks (SCENIC)")
    parser.add_argument("--input", required=True, help="Path to .h5ad")
    parser.add_argument("--output", default="results/15_grn.h5ad")
    parser.add_argument("--tf_list", default="allTFs_hg38.txt",
                        help="Path to TF gene list (one gene per line)")
    parser.add_argument("--db_dir", default="cisTarget_databases",
                        help="Directory with cisTarget .feather databases")
    parser.add_argument("--motif_annotations",
                        default="motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
                        help="Motif annotation file for cisTarget pruning")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
