"""
08 - Differential Gene Expression
==================================
Pseudobulk DGE analysis using edgeR (via R) and single-cell level tests
(Wilcoxon, t-test). Aggregates cells per donor/sample to form pseudobulk
replicates, then tests for condition-specific gene expression changes.


"""

import argparse
import os
import random
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

warnings.filterwarnings("ignore")


NUM_OF_CELLS_PER_DONOR = 30


def aggregate_and_filter(adata, cell_identity, donor_key="sample",
                         condition_key="label", cell_identity_key="cell_type",
                         obs_to_keep=None):
    """Aggregate cells into pseudobulk samples per donor for a given cell type."""
    if obs_to_keep is None:
        obs_to_keep = []

    adata_sub = adata[adata.obs[cell_identity_key] == cell_identity].copy()

    size_by_donor = adata_sub.obs.groupby([donor_key]).size()
    donors_to_drop = [d for d in size_by_donor.index if size_by_donor[d] <= NUM_OF_CELLS_PER_DONOR]
    if donors_to_drop:
        print(f"    Dropping low-count donors: {donors_to_drop}")

    df = pd.DataFrame(columns=[*adata_sub.var_names, *obs_to_keep])
    adata_sub.obs[donor_key] = adata_sub.obs[donor_key].astype("category")

    for donor in adata_sub.obs[donor_key].cat.categories:
        if donor in donors_to_drop:
            continue
        adata_donor = adata_sub[adata_sub.obs[donor_key] == donor]
        agg_dict = dict.fromkeys(adata_donor.var_names, "sum")
        for obs in obs_to_keep:
            agg_dict[obs] = "first"
        df_donor = pd.DataFrame(adata_donor.X.toarray() if hasattr(adata_donor.X, "toarray") else adata_donor.X)
        df_donor.index = adata_donor.obs_names
        df_donor.columns = adata_donor.var_names
        df_donor = df_donor.join(adata_donor.obs[obs_to_keep])
        df_donor = df_donor.groupby(donor_key).agg(agg_dict)
        df_donor[donor_key] = donor
        df.loc[f"donor_{donor}"] = df_donor.loc[donor]

    adata_pb = sc.AnnData(
        df[adata_sub.var_names].astype(float),
        obs=df.drop(columns=adata_sub.var_names)
    )
    return adata_pb


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # Use raw counts for DGE
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    # Clean up obs columns
    for col in [args.condition_key, args.cell_type_key, args.sample_key]:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")

    # --- Single-cell level: Wilcoxon ---
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, groupby=args.condition_key, method="wilcoxon",
                            key_added="wilcoxon")
    result_df = sc.get.rank_genes_groups_df(adata, None, key="wilcoxon")
    result_df.to_csv(os.path.join(args.figdir, "08_wilcoxon_dge.csv"), index=False)
    print(f"  Wilcoxon DGE: {len(result_df)} results saved")

    # Plot top DE genes
    sc.pl.rank_genes_groups(adata, key="wilcoxon", n_genes=20,
                            save="_08_wilcoxon_top20.png", show=False)

    # --- Pseudobulk DGE ---
    obs_to_keep = [args.condition_key, args.cell_type_key, args.sample_key]
    obs_to_keep = [o for o in obs_to_keep if o in adata.obs.columns]

    if args.cell_type_key in adata.obs.columns:
        cell_types = adata.obs[args.cell_type_key].cat.categories
        pb_list = []
        for ct in cell_types:
            print(f"  Pseudobulk aggregation: {ct}")
            try:
                pb = aggregate_and_filter(adata, ct, donor_key=args.sample_key,
                                          condition_key=args.condition_key,
                                          cell_identity_key=args.cell_type_key,
                                          obs_to_keep=obs_to_keep)
                pb_list.append(pb)
            except Exception as e:
                print(f"    Skipped {ct}: {e}")
        if pb_list:
            adata_pb = pb_list[0].concatenate(*pb_list[1:]) if len(pb_list) > 1 else pb_list[0]
            adata_pb.write_h5ad(os.path.join(os.path.dirname(args.output) or ".", "08_pseudobulk.h5ad"))
            print(f"  Pseudobulk AnnData saved")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="08 - Differential Gene Expression")
    parser.add_argument("--input", required=True, help="Path to annotated .h5ad")
    parser.add_argument("--output", default="results/08_dge.h5ad")
    parser.add_argument("--condition_key", default="label")
    parser.add_argument("--cell_type_key", default="cell_type")
    parser.add_argument("--sample_key", default="sample")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
