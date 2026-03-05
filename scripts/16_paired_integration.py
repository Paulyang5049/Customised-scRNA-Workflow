"""
16 - Paired Multi-Modal Integration (CITE-seq)
================================================
Integrate paired RNA + protein (ADT) data from CITE-seq using:
  1. muon (MuData) for joint representation
  2. Weighted Nearest Neighbors (WNN) via Seurat (R)


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")


def main(args):
    print(f"Loading RNA from {args.rna_input}")
    adata_rna = sc.read_h5ad(args.rna_input)
    print(f"  RNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # Load protein (ADT) data
    if args.adt_input and os.path.isfile(args.adt_input):
        adata_adt = sc.read_h5ad(args.adt_input)
        print(f"  ADT: {adata_adt.n_obs} cells x {adata_adt.n_vars} proteins")
    else:
        print("  No ADT file provided or found. Running RNA-only workflow.")
        adata_rna.write_h5ad(args.output)
        return

    try:
        import muon as mu

        # Build MuData
        mdata = mu.MuData({"rna": adata_rna, "adt": adata_adt})
        print(f"  MuData: {mdata}")

        # Process RNA modality
        mu.pp.intersect_obs(mdata)
        sc.pp.normalize_total(mdata.mod["rna"], target_sum=1e4)
        sc.pp.log1p(mdata.mod["rna"])
        sc.pp.highly_variable_genes(mdata.mod["rna"])
        sc.tl.pca(mdata.mod["rna"])

        # Process ADT modality (CLR normalization)
        mu.prot.pp.clr(mdata.mod["adt"])
        sc.tl.pca(mdata.mod["adt"])

        # Multi-omics factor analysis or WNN
        mu.pp.neighbors(mdata.mod["rna"])
        mu.pp.neighbors(mdata.mod["adt"])
        mu.tl.umap(mdata)

        # Plot joint embedding
        mu.pl.umap(mdata, color=[args.cluster_key] if args.cluster_key in
                   mdata.obs.columns else None)
        plt.savefig(os.path.join(args.figdir, "16_muon_umap.png"), bbox_inches="tight")
        plt.close()

        # WNN via muon
        try:
            mu.pp.neighbors(mdata, key_added="wnn")
            mu.tl.umap(mdata, neighbors_key="wnn")
            mu.pl.umap(mdata, color=[args.cluster_key] if args.cluster_key in
                       mdata.obs.columns else None, neighbors_key="wnn")
            plt.savefig(os.path.join(args.figdir, "16_wnn_umap.png"), bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"  WNN skipped: {e}")

        # Save as MuData
        mdata.write(args.output)
        print(f"  Saved MuData to {args.output}")

    except ImportError:
        print("  muon not installed. Running basic concatenation.")
        adata_rna.write_h5ad(args.output)

    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="16 - Paired Integration (CITE-seq)")
    parser.add_argument("--rna_input", required=True, help="Path to RNA .h5ad")
    parser.add_argument("--adt_input", default=None, help="Path to ADT .h5ad")
    parser.add_argument("--output", default="results/16_paired_integration.h5mu")
    parser.add_argument("--cluster_key", default="cell_type")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
