"""
05 - Data Integration & Batch Correction
=========================================
Compare multiple batch correction methods:
  - scVI / scANVI (deep learning)
  - BBKNN (graph-based)
  - Seurat CCA (requires R)
Then benchmark with scIB metrics.


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir
    batch_key = args.batch_key
    label_key = args.label_key

    # Prepare: normalize and select HVGs per batch
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    # Unintegrated baseline
    sc.tl.pca(adata_hvg)
    sc.pp.neighbors(adata_hvg)
    sc.tl.umap(adata_hvg)
    sc.pl.umap(adata_hvg, color=[label_key, batch_key], wspace=1,
               save="_05_unintegrated.png", show=False)

    results = {"Unintegrated": adata_hvg}

    # --- scVI ---
    try:
        import scvi
        adata_scvi = adata_hvg.copy()
        scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_key)
        model = scvi.model.SCVI(adata_scvi)
        max_epochs = int(np.min([round((20000 / adata.n_obs) * 400), 400]))
        model.train()
        adata_scvi.obsm["X_scVI"] = model.get_latent_representation()
        sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
        sc.tl.umap(adata_scvi)
        sc.pl.umap(adata_scvi, color=[label_key, batch_key], wspace=1,
                   save="_05_scVI.png", show=False)
        results["scVI"] = adata_scvi
        print("  scVI: done")
    except Exception as e:
        print(f"  scVI skipped: {e}")

    # --- BBKNN ---
    try:
        import bbknn
        adata_bbknn = adata_hvg.copy()
        adata_bbknn.X = adata_bbknn.layers["logcounts"].copy()
        sc.pp.pca(adata_bbknn)
        nwb = 25 if adata_hvg.n_obs > 100000 else 3
        bbknn.bbknn(adata_bbknn, batch_key=batch_key, neighbors_within_batch=nwb)
        sc.tl.umap(adata_bbknn)
        sc.pl.umap(adata_bbknn, color=[label_key, batch_key], wspace=1,
                   save="_05_bbknn.png", show=False)
        results["BBKNN"] = adata_bbknn
        print("  BBKNN: done")
    except Exception as e:
        print(f"  BBKNN skipped: {e}")

    # --- Benchmarking with scIB ---
    try:
        import scib
        metrics_list = []
        for name, ad in results.items():
            embed = "X_scVI" if name == "scVI" else None
            m = scib.metrics.metrics_fast(adata, ad, batch_key, label_key, embed=embed)
            metrics_list.append(m.rename(columns={m.columns[0]: name}))
        metrics = pd.concat(metrics_list, axis=1)
        metrics.to_csv(os.path.join(args.figdir, "05_integration_metrics.csv"))
        print(f"  Metrics saved to {args.figdir}/05_integration_metrics.csv")
    except Exception as e:
        print(f"  scIB benchmarking skipped: {e}")

    # Save the best integrated result (or unintegrated fallback)
    best = results.get("scVI", results["Unintegrated"])
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    best.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="05 - Data Integration & Batch Correction")
    parser.add_argument("--input", required=True, help="Path to .h5ad with batch info")
    parser.add_argument("--output", default="results/05_integrated.h5ad")
    parser.add_argument("--batch_key", default="batch")
    parser.add_argument("--label_key", default="cell_type")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
