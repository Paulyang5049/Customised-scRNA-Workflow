"""
07 - Cell Type Annotation
=========================
Annotate clusters using:
  1. Manual marker gene inspection (dotplots, UMAPs)
  2. CellTypist automated annotation
  3. scArches label transfer from a reference atlas


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore", category=DeprecationWarning)


# Common marker genes for bone-marrow / PBMC cell types
MARKER_GENES = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    "cDC1": ["CLEC9A", "CADM1"],
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Normoblast": ["SLC4A1", "SLC25A37", "HBB", "HBA2", "HBA1", "TFRC"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    "CD4+ T naive": ["CD4", "IL7R", "TRBC2", "CCR7"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
    "HSC": ["NRIP1", "MECOM", "PROM1", "CD34"],
}


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir
    sc.set_figure_params(figsize=(5, 5))

    # Use normalized layer if available
    if "scran_normalization" in adata.layers:
        adata.X = adata.layers["scran_normalization"]
    elif "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"]

    # Ensure embedding exists
    if "X_umap" not in adata.obsm:
        if "highly_deviant" in adata.var.columns:
            adata.var["highly_variable"] = adata.var["highly_deviant"]
        sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    # Cluster at the annotation resolution
    cluster_key = f"leiden_{args.resolution}"
    sc.tl.leiden(adata, resolution=args.resolution, key_added=cluster_key)
    sc.pl.umap(adata, color=cluster_key, legend_loc="on data",
               save=f"_07_clusters_res{args.resolution}.png", show=False)

    # Filter marker genes to those present in the data
    marker_genes_in_data = {
        ct: [m for m in markers if m in adata.var_names]
        for ct, markers in MARKER_GENES.items()
    }
    marker_genes_in_data = {k: v for k, v in marker_genes_in_data.items() if v}

    # Dotplot of markers
    if marker_genes_in_data:
        sc.pl.dotplot(adata, groupby=cluster_key, var_names=marker_genes_in_data,
                      standard_scale="var", save="_07_markers_dotplot.png", show=False)

    # --- CellTypist automated annotation ---
    try:
        import celltypist
        from celltypist import models

        model = models.download_model(args.celltypist_model, force_update=False)
        model = models.Model.load(model=args.celltypist_model)
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        adata.obs["celltypist_pred"] = predictions.predicted_labels["majority_voting"]
        sc.pl.umap(adata, color="celltypist_pred",
                   save="_07_celltypist.png", show=False)
        print(f"  CellTypist annotation: done ({adata.obs['celltypist_pred'].nunique()} types)")
    except Exception as e:
        print(f"  CellTypist skipped: {e}")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="07 - Cell Type Annotation")
    parser.add_argument("--input", default="results/06_clustered.h5ad")
    parser.add_argument("--output", default="results/07_annotated.h5ad")
    parser.add_argument("--resolution", type=float, default=2.0)
    parser.add_argument("--celltypist_model", default="Immune_All_Low.pkl")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
