"""Validate an AnnData input before running the scRNA-seq workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc
from workflow_utils import get_counts, write_json


def main(args: argparse.Namespace) -> None:
    input_path = Path(args.input)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    adata = sc.read_h5ad(input_path, backed="r")
    problems: list[str] = []
    warnings: list[str] = []

    if adata.n_obs == 0 or adata.n_vars == 0:
        problems.append("The object must contain at least one barcode and one gene.")
    if not adata.obs_names.is_unique:
        problems.append("Cell/barcode names are not unique.")
    if not adata.var_names.is_unique:
        problems.append("Gene names are not unique; resolve identifiers before analysis.")

    required = [item for item in args.required_obs.split(",") if item]
    for column in required:
        if column not in adata.obs:
            problems.append(f"Missing required obs column: {column}")
        elif adata.obs[column].isna().any():
            problems.append(f"obs[{column!r}] contains missing values.")

    if args.sample_key in adata.obs:
        sizes = adata.obs[args.sample_key].value_counts(dropna=False)
        if len(sizes) < 2:
            warnings.append("Only one sample/capture partition is present.")
        if (sizes < args.min_cells_per_sample).any():
            small = sizes[sizes < args.min_cells_per_sample].to_dict()
            warnings.append(f"Small sample partitions detected: {small}")

    try:
        get_counts(adata, args.counts_layer)
    except (KeyError, ValueError) as exc:
        problems.append(str(exc))

    report = {
        "input": str(input_path.resolve()),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": adata.obs.columns.tolist(),
        "layers": list(adata.layers.keys()),
        "problems": problems,
        "warnings": warnings,
        "valid": not problems,
    }
    write_json(args.report, report)

    for warning in warnings:
        print(f"WARNING: {warning}")
    if problems:
        for problem in problems:
            print(f"ERROR: {problem}")
        raise SystemExit(2)
    print(f"Input is valid: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Validation report: {args.report}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a raw-count AnnData input")
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--report", default="results/00_validation.json")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument(
        "--required-obs",
        default="sample",
        help="Comma-separated required adata.obs columns",
    )
    parser.add_argument("--min-cells-per-sample", type=int, default=100)
    main(parser.parse_args())

