"""Correct droplet-derived ambient RNA with SoupX using an unfiltered matrix."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread, mmwrite
from workflow_utils import get_counts, record_provenance


def _clusters(adata, counts, cluster_key: str, seed: int) -> pd.Series:
    if cluster_key in adata.obs:
        return adata.obs[cluster_key].astype(str)
    if adata.n_obs < 30:
        raise ValueError(
            f"obs[{cluster_key!r}] is absent and there are too few cells for preliminary clustering."
        )

    working = sc.AnnData(counts.copy(), obs=adata.obs[[column for column in adata.obs if column]])
    working.var_names = adata.var_names.copy()
    sc.pp.normalize_total(working, target_sum=1e4)
    sc.pp.log1p(working)
    sc.pp.highly_variable_genes(
        working,
        n_top_genes=min(2000, working.n_vars),
        flavor="seurat",
        subset=False,
    )
    n_comps = min(30, working.n_obs - 1, int(working.var["highly_variable"].sum()) - 1)
    if n_comps < 2:
        raise ValueError("Insufficient cells or variable genes for preliminary SoupX clustering.")
    sc.tl.pca(working, n_comps=n_comps, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(working, n_pcs=n_comps, random_state=seed)
    sc.tl.leiden(
        working,
        resolution=0.5,
        key_added=cluster_key,
        flavor="igraph",
        n_iterations=2,
        random_state=seed,
    )
    return working.obs[cluster_key].astype(str)


def main(args: argparse.Namespace) -> None:
    if shutil.which(args.rscript) is None:
        raise RuntimeError("Rscript is unavailable; SoupX cannot run. Activate environment.yml.")

    adata = sc.read_h5ad(args.input)
    raw = sc.read_h5ad(args.raw_droplets)
    if not adata.obs_names.isin(raw.obs_names).all():
        raise ValueError("Every filtered-cell barcode must exist in the unfiltered droplet object.")
    if not adata.var_names.isin(raw.var_names).all():
        missing = adata.var_names[~adata.var_names.isin(raw.var_names)].tolist()[:10]
        raise ValueError(f"Filtered genes are absent from the raw droplet object: {missing}")

    filtered_counts = get_counts(adata, args.counts_layer)
    raw_counts = get_counts(raw, args.raw_counts_layer)
    filtered_counts = sparse.csr_matrix(filtered_counts)
    raw_counts = sparse.csr_matrix(raw_counts[:, raw.var_names.get_indexer(adata.var_names)])
    clusters = _clusters(adata, filtered_counts, args.cluster_key, args.seed)
    script = Path(__file__).with_name("run_soupx.R")

    with tempfile.TemporaryDirectory(prefix="scrna_soupx_") as temporary:
        workdir = Path(temporary)
        paths = {
            "raw": workdir / "raw.mtx",
            "filtered": workdir / "filtered.mtx",
            "genes": workdir / "genes.txt",
            "raw_cells": workdir / "raw_cells.txt",
            "cells": workdir / "cells.txt",
            "clusters": workdir / "clusters.tsv",
            "corrected": workdir / "corrected.mtx",
            "contamination": workdir / "contamination.tsv",
        }
        mmwrite(paths["raw"], raw_counts.T.tocoo())
        mmwrite(paths["filtered"], filtered_counts.T.tocoo())
        paths["genes"].write_text("\n".join(adata.var_names.astype(str)) + "\n")
        paths["raw_cells"].write_text("\n".join(raw.obs_names.astype(str)) + "\n")
        paths["cells"].write_text("\n".join(adata.obs_names.astype(str)) + "\n")
        pd.DataFrame(
            {"barcode": adata.obs_names.astype(str), "cluster": clusters.to_numpy()}
        ).to_csv(paths["clusters"], sep="\t", index=False)

        subprocess.run(
            [
                args.rscript,
                str(script),
                str(paths["raw"]),
                str(paths["filtered"]),
                str(paths["genes"]),
                str(paths["raw_cells"]),
                str(paths["cells"]),
                str(paths["clusters"]),
                str(paths["corrected"]),
                str(paths["contamination"]),
                str(args.seed),
            ],
            check=True,
        )
        corrected = sparse.csr_matrix(mmread(paths["corrected"]).T)
        contamination = pd.read_csv(paths["contamination"], sep="\t").set_index("barcode")

    if corrected.shape != adata.shape:
        raise RuntimeError(f"SoupX returned shape {corrected.shape}; expected {adata.shape}.")
    contamination = contamination.reindex(adata.obs_names.astype(str))
    if contamination["contamination_fraction"].isna().any():
        raise RuntimeError("SoupX did not return contamination estimates for every cell.")

    adata.layers["soupx_counts"] = corrected
    adata.obs["ambient_contamination_fraction"] = contamination[
        "contamination_fraction"
    ].to_numpy(dtype=float)
    if args.max_contamination is None:
        adata.obs["passes_ambient_qc"] = True
    else:
        adata.obs["passes_ambient_qc"] = (
            adata.obs["ambient_contamination_fraction"] <= args.max_contamination
        )

    if args.set_active_counts:
        if "counts_pre_ambient" not in adata.layers:
            adata.layers["counts_pre_ambient"] = filtered_counts.copy()
        adata.layers["counts"] = corrected.copy()

    record_provenance(
        adata,
        "ambient_rna_correction",
        {
            "method": "SoupX",
            "raw_droplets": str(Path(args.raw_droplets).resolve()),
            "cluster_key": args.cluster_key,
            "max_contamination": args.max_contamination,
            "set_active_counts": args.set_active_counts,
            "seed": args.seed,
            "median_contamination": float(
                np.median(adata.obs["ambient_contamination_fraction"])
            ),
        },
        ["anndata", "scanpy", "numpy", "pandas", "scipy"],
    )

    if args.filter_high_contamination:
        if args.max_contamination is None:
            raise ValueError("--filter-high-contamination requires --max-contamination.")
        adata = adata[adata.obs["passes_ambient_qc"]].copy()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(
        adata.obs["ambient_contamination_fraction"].describe(
            percentiles=[0.5, 0.9, 0.95, 0.99]
        )
    )
    print(f"Saved: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct ambient RNA with SoupX")
    parser.add_argument("--input", required=True, help="Filtered-cell .h5ad")
    parser.add_argument("--raw-droplets", required=True, help="Unfiltered droplet .h5ad")
    parser.add_argument("--output", default="results/01c_soupx.h5ad")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--raw-counts-layer", default="counts")
    parser.add_argument("--cluster-key", default="soupx_groups")
    parser.add_argument("--max-contamination", type=float, default=None)
    parser.add_argument("--filter-high-contamination", action="store_true")
    parser.add_argument("--set-active-counts", action="store_true")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--rscript", default="Rscript")
    main(parser.parse_args())

