"""
09 - Gene Set Enrichment Analysis (GSEA)
=========================================
Run pathway enrichment analysis using decoupler:
  - Preranked GSEA on DE statistics
  - AUCell activity scoring per cell
Uses Reactome pathways by default.


"""

import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")


def gmt_to_decoupler(pth):
    """Parse a GMT file into a decoupler-compatible DataFrame."""
    from itertools import chain, repeat
    pathways = {}
    with Path(pth).open("r") as f:
        for line in f:
            name, _, *genes = line.strip().split("\t")
            pathways[name] = genes
    return pd.DataFrame.from_records(
        chain.from_iterable(zip(repeat(k), v) for k, v in pathways.items()),
        columns=["geneset", "genesymbol"],
    )


def main(args):
    print(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    os.makedirs(args.figdir, exist_ok=True)
    sc.settings.figdir = args.figdir

    import decoupler

    # Normalize if not already
    if "counts" in adata.layers:
        adata.layers["counts"] = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Find HVGs for GSEA
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3",
                                    subset=False, layer="counts")

    # DE for GSEA input
    group_key = args.group_key
    if group_key not in adata.obs.columns:
        # Build group column from condition + cell_type
        adata.obs["group"] = (
            adata.obs[args.condition_key].astype(str) + "_" + adata.obs[args.cell_type_key].astype(str)
        )
        group_key = "group"

    sc.tl.rank_genes_groups(adata, group_key, method="t-test", key_added="t-test")

    # Load pathway database
    if not Path(args.gmt_file).is_file():
        print(f"  GMT file not found at {args.gmt_file}, skipping GSEA.")
        adata.write_h5ad(args.output)
        return

    reactome = gmt_to_decoupler(args.gmt_file)
    geneset_size = reactome.groupby("geneset").size()
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 500)]
    reactome_filtered = reactome[reactome["geneset"].isin(gsea_genesets)]

    # Run preranked GSEA for each group
    groups = adata.obs[group_key].unique()
    all_results = []
    for grp in groups[:5]:  # limit to first 5 for speed
        try:
            df = sc.get.rank_genes_groups_df(adata, grp, key="t-test")
            t_stats = df.set_index("names").loc[adata.var["highly_variable"]].sort_values(
                "scores", key=np.abs, ascending=False)[["scores"]].rename_axis([grp], axis=1)

            scores, norm, pvals = decoupler.run_gsea(
                t_stats.T, reactome_filtered, source="geneset", target="genesymbol"
            )
            res = pd.concat({"score": scores.T, "norm": norm.T, "pval": pvals.T}, axis=1).droplevel(1, axis=1)
            res["group"] = grp
            all_results.append(res.sort_values("pval").head(20))
        except Exception as e:
            print(f"    GSEA for {grp} failed: {e}")

    if all_results:
        combined = pd.concat(all_results)
        combined.to_csv(os.path.join(args.figdir, "09_gsea_results.csv"))
        print(f"  GSEA results saved to {args.figdir}/09_gsea_results.csv")

    # AUCell per-cell activity scores
    try:
        decoupler.run_aucell(adata, reactome_filtered, source="geneset",
                             target="genesymbol", use_raw=False)
        print(f"  AUCell scores computed for {adata.obsm['aucell_estimate'].shape[1]} gene sets")
    except Exception as e:
        print(f"  AUCell skipped: {e}")

    # Save
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="09 - Gene Set Enrichment Analysis")
    parser.add_argument("--input", required=True, help="Path to .h5ad with DE results or raw data")
    parser.add_argument("--output", default="results/09_gsea.h5ad")
    parser.add_argument("--gmt_file", default="c2.cp.reactome.v7.5.1.symbols.gmt",
                        help="Path to Reactome GMT file")
    parser.add_argument("--group_key", default="group", help="obs column used for group comparisons")
    parser.add_argument("--condition_key", default="condition")
    parser.add_argument("--cell_type_key", default="cell_type")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
