"""
10 - Compositional Analysis
============================
Test whether cell-type proportions change between conditions using
scCODA (Bayesian compositional model via pertpy).


"""

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

warnings.filterwarnings("ignore")


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)

    import pertpy as pt

    # Build scCODA model
    sccoda_model = pt.tl.Sccoda()
    sccoda_data = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=args.cell_type_key,
        sample_identifier=args.sample_key,
        covariate_obs=[args.condition_key],
    )
    print(f"  scCODA data prepared: {sccoda_data}")

    # Boxplots of cell type compositions
    sccoda_model.plot_boxplots(
        sccoda_data, modality_key="coda", feature_name=args.condition_key,
        add_dots=True
    )
    plt.savefig(os.path.join(args.figdir, "10_composition_boxplots.png"), bbox_inches="tight")
    plt.close()

    # Stacked barplot
    sccoda_model.plot_stacked_barplot(
        sccoda_data, modality_key="coda", feature_name=args.condition_key
    )
    plt.savefig(os.path.join(args.figdir, "10_composition_stacked.png"), bbox_inches="tight")
    plt.close()

    # Run the model
    sccoda_data = sccoda_model.prepare(
        sccoda_data, modality_key="coda", formula=f"C({args.condition_key})",
        reference_cell_type=args.reference_cell_type
    )
    sccoda_model.run_nuts(sccoda_data, modality_key="coda")
    sccoda_model.set_fdr(sccoda_data, 0.4)

    # Credible effects
    try:
        results = sccoda_model.credible_effects(sccoda_data, modality_key="coda")
        print(f"  Credible effects:\n{results}")
        results.to_csv(os.path.join(args.figdir, "10_credible_effects.csv"))
    except Exception as e:
        print(f"  Could not extract credible effects: {e}")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="10 - Compositional Analysis")
    parser.add_argument("--input", required=True, help="Path to annotated .h5ad")
    parser.add_argument("--output", default="results/10_compositional.h5ad")
    parser.add_argument("--cell_type_key", default="cell_label")
    parser.add_argument("--sample_key", default="batch")
    parser.add_argument("--condition_key", default="condition")
    parser.add_argument("--reference_cell_type", default="automatic",
                        help="Reference cell type for scCODA")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
