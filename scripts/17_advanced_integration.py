"""
17 - Advanced Multi-Modal Integration
=======================================
Integrate unpaired or partially overlapping modalities using:
  1. TOTALVI  (scvi-tools) — joint RNA + protein
  2. GLUE     (scglue)     — graph-linked unified embedding
  3. MultiGrate           — multi-modal autoencoder

"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")


def run_totalvi(adata, batch_key, protein_key, figdir):
    """Joint RNA+protein integration with TOTALVI."""
    import scvi

    # Assume protein expression stored in adata.obsm[protein_key]
    if protein_key not in adata.obsm:
        print(f"    TOTALVI skipped: {protein_key} not in adata.obsm")
        return adata

    scvi.model.TOTALVI.setup_anndata(
        adata, batch_key=batch_key,
        protein_expression_obsm_key=protein_key,
    )
    model = scvi.model.TOTALVI(adata, latent_distribution="normal")
    model.train(max_epochs=100, early_stopping=True)

    adata.obsm["X_totalvi"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_totalvi")
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[batch_key], save="_17_totalvi.png", show=False)
    print("    TOTALVI: done")
    return adata


def run_glue(adata_rna, adata_atac, figdir):
    """Cross-modality integration with GLUE."""
    try:
        import scglue
    except ImportError:
        print("    GLUE skipped: scglue not installed")
        return adata_rna

    scglue.data.get_gene_annotation(adata_rna, gtf=None, gtf_by="gene_name")
    guidance = scglue.genomics.rna_anchored_guidance_graph(adata_rna, adata_atac)

    scglue.models.configure_dataset(adata_rna, "NB", use_highly_variable=True,
                                    use_rep="X_pca")
    scglue.models.configure_dataset(adata_atac, "NB", use_highly_variable=True,
                                    use_rep="X_lsi")

    glue = scglue.models.fit_SCGLUE(
        {"rna": adata_rna, "atac": adata_atac},
        guidance, fit_kws={"directory": "glue_model"}
    )
    adata_rna.obsm["X_glue"] = glue.encode_data("rna", adata_rna)
    sc.pp.neighbors(adata_rna, use_rep="X_glue")
    sc.tl.umap(adata_rna)
    sc.pl.umap(adata_rna, save="_17_glue.png", show=False)
    print("    GLUE: done")
    return adata_rna


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    # --- TOTALVI ---
    if args.run_totalvi:
        print("  Running TOTALVI...")
        adata = run_totalvi(adata, args.batch_key, args.protein_key, args.figdir)

    # --- GLUE (requires separate ATAC file) ---
    if args.run_glue and args.atac_input and os.path.isfile(args.atac_input):
        print("  Running GLUE...")
        adata_atac = sc.read_h5ad(args.atac_input)
        adata = run_glue(adata, adata_atac, args.figdir)

    # --- MultiGrate ---
    if args.run_multigrate:
        try:
            import multigrate as mtg

            print("  Running MultiGrate...")
            mtg.model.MultiVAE.setup_anndata(adata, batch_key=args.batch_key)
            model = mtg.model.MultiVAE(adata)
            model.train(max_epochs=100)
            adata.obsm["X_multigrate"] = model.get_latent_representation()
            sc.pp.neighbors(adata, use_rep="X_multigrate")
            sc.tl.umap(adata)
            sc.pl.umap(adata, color=[args.batch_key],
                       save="_17_multigrate.png", show=False)
            print("    MultiGrate: done")
        except ImportError:
            print("    MultiGrate skipped: multigrate not installed")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="17 - Advanced Multi-Modal Integration")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results/17_advanced_integration.h5ad")
    parser.add_argument("--batch_key", default="batch")
    parser.add_argument("--protein_key", default="protein_expression",
                        help="obsm key for protein data (TOTALVI)")
    parser.add_argument("--atac_input", default=None, help="Path to ATAC .h5ad (GLUE)")
    parser.add_argument("--run_totalvi", action="store_true")
    parser.add_argument("--run_glue", action="store_true")
    parser.add_argument("--run_multigrate", action="store_true")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
