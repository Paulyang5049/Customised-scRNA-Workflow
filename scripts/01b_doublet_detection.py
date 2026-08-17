"""Sample-aware doublet detection using the Bioconductor scDblFinder package."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite
from workflow_utils import get_counts, record_provenance, require_columns


def main(args: argparse.Namespace) -> None:
    if shutil.which(args.rscript) is None:
        raise RuntimeError(
            "Rscript is unavailable. Create and activate environment.yml; "
            "scDblFinder will not be silently replaced with another caller."
        )

    adata = sc.read_h5ad(args.input)
    require_columns(adata, [args.sample_key])
    if adata.obs[args.sample_key].isna().any():
        raise ValueError(f"obs[{args.sample_key!r}] contains missing sample assignments.")

    counts = get_counts(adata, args.counts_layer)
    counts = counts if sparse.issparse(counts) else sparse.csr_matrix(counts)
    script = Path(__file__).with_name("run_scdblfinder.R")

    with tempfile.TemporaryDirectory(prefix="scrna_scdblfinder_") as temporary:
        workdir = Path(temporary)
        matrix_path = workdir / "counts_genes_by_cells.mtx"
        samples_path = workdir / "samples.tsv"
        result_path = workdir / "doublets.tsv"

        mmwrite(matrix_path, counts.T.tocoo())
        pd.DataFrame(
            {
                "barcode": adata.obs_names.astype(str),
                "sample": adata.obs[args.sample_key].astype(str).to_numpy(),
            }
        ).to_csv(samples_path, sep="\t", index=False)

        expected_rate = "auto" if args.expected_rate is None else str(args.expected_rate)
        command = [
            args.rscript,
            str(script),
            str(matrix_path),
            str(samples_path),
            str(result_path),
            str(args.seed),
            expected_rate,
        ]
        subprocess.run(command, check=True)
        calls = pd.read_csv(result_path, sep="\t").set_index("barcode")

    if calls.index.tolist() != adata.obs_names.astype(str).tolist():
        raise RuntimeError("scDblFinder returned barcodes in an unexpected order.")

    adata.obs["scDblFinder_score"] = calls["score"].to_numpy(dtype=float)
    adata.obs["scDblFinder_class"] = pd.Categorical(calls["class"].astype(str))
    adata.obs["scDblFinder_partition"] = pd.Categorical(calls["sample"].astype(str))
    adata.obs["passes_doublet_qc"] = adata.obs["scDblFinder_class"].eq("singlet")

    summary = (
        adata.obs.groupby(["scDblFinder_partition", "scDblFinder_class"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)

    figure_path = Path(args.figure)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    for sample, frame in adata.obs.groupby("scDblFinder_partition", observed=True):
        ax.hist(frame["scDblFinder_score"], bins=50, alpha=0.45, label=str(sample))
    ax.set(xlabel="scDblFinder score", ylabel="Cells", title="Doublet scores by partition")
    if adata.obs["scDblFinder_partition"].nunique() <= 12:
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)

    record_provenance(
        adata,
        "doublet_detection",
        {
            "method": "scDblFinder",
            "sample_key": args.sample_key,
            "counts_layer": args.counts_layer,
            "expected_rate": expected_rate,
            "seed": args.seed,
            "filter_doublets": args.filter_doublets,
            "n_doublets": int((~adata.obs["passes_doublet_qc"]).sum()),
        },
        ["anndata", "scanpy", "numpy", "pandas", "scipy"],
    )

    if args.filter_doublets:
        adata = adata[adata.obs["passes_doublet_qc"]].copy()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(summary.to_string(index=False))
    print(f"Saved: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Call doublets per sample with scDblFinder")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results/01b_doublets.h5ad")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--expected-rate", type=float, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--filter-doublets", action="store_true")
    parser.add_argument("--summary", default="results/01b_doublet_summary.csv")
    parser.add_argument("--figure", default="figures/01b_doublet_scores.png")
    parser.add_argument("--rscript", default="Rscript")
    main(parser.parse_args())

