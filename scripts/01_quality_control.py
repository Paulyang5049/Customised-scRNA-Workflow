"""Data-driven, sample-aware cell QC with provenance-preserving flags."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
from workflow_utils import get_counts, record_provenance


def mad_flag(values: pd.Series, nmads: float, side: str) -> pd.Series:
    """Flag robust outliers on one side of a distribution."""
    median = float(np.median(values))
    mad = float(median_abs_deviation(values))
    if not np.isfinite(mad) or mad == 0:
        return pd.Series(False, index=values.index)
    if side == "low":
        return values < median - nmads * mad
    if side == "high":
        return values > median + nmads * mad
    raise ValueError("side must be 'low' or 'high'")


def grouped_flag(adata, metric: str, sample_key: str, nmads: float, side: str) -> pd.Series:
    if sample_key in adata.obs:
        return (
            adata.obs.groupby(sample_key, observed=True, group_keys=False)[metric]
            .apply(lambda values: mad_flag(values, nmads, side))
            .reindex(adata.obs_names)
            .fillna(False)
            .astype(bool)
        )
    return mad_flag(adata.obs[metric], nmads, side)


def threshold_table(adata, sample_key: str, nmads: float, mt_nmads: float, max_pct_mt):
    rows = []
    for sample, frame in adata.obs.groupby(sample_key, observed=True):
        counts_median = np.median(frame["log1p_total_counts"])
        counts_mad = median_abs_deviation(frame["log1p_total_counts"])
        genes_median = np.median(frame["log1p_n_genes_by_counts"])
        genes_mad = median_abs_deviation(frame["log1p_n_genes_by_counts"])
        mt_median = np.median(frame["pct_counts_mt"])
        mt_mad = median_abs_deviation(frame["pct_counts_mt"])
        mt_high = mt_median + mt_nmads * mt_mad
        if max_pct_mt is not None:
            mt_high = min(mt_high, max_pct_mt)
        rows.append(
            {
                sample_key: str(sample),
                "min_total_counts": max(0.0, float(np.expm1(counts_median - nmads * counts_mad))),
                "max_total_counts": float(np.expm1(counts_median + nmads * counts_mad)),
                "min_genes": max(0.0, float(np.expm1(genes_median - nmads * genes_mad))),
                "max_genes": float(np.expm1(genes_median + nmads * genes_mad)),
                "max_pct_mt": float(mt_high),
            }
        )
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    adata = sc.read_h5ad(args.input)
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell and gene names must be unique; run 00_validate_input.py first.")
    if args.sample_key not in adata.obs:
        adata.obs[args.sample_key] = "all_cells"
        print(f"WARNING: obs[{args.sample_key!r}] absent; QC thresholds are global.")

    counts = get_counts(adata, args.counts_layer)
    if args.counts_layer not in adata.layers:
        adata.layers[args.counts_layer] = counts.copy()

    prefixes = tuple(prefix for prefix in args.mt_prefixes.split(",") if prefix)
    adata.var["mt"] = adata.var_names.str.startswith(prefixes)
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        layer=args.counts_layer,
        percent_top=None,
        inplace=True,
        log1p=True,
    )

    adata.obs["qc_low_counts"] = grouped_flag(
        adata, "log1p_total_counts", args.sample_key, args.nmads, "low"
    )
    adata.obs["qc_high_counts"] = grouped_flag(
        adata, "log1p_total_counts", args.sample_key, args.nmads, "high"
    )
    adata.obs["qc_low_genes"] = grouped_flag(
        adata, "log1p_n_genes_by_counts", args.sample_key, args.nmads, "low"
    )
    adata.obs["qc_high_genes"] = grouped_flag(
        adata, "log1p_n_genes_by_counts", args.sample_key, args.nmads, "high"
    )
    adata.obs["qc_high_mito"] = grouped_flag(
        adata, "pct_counts_mt", args.sample_key, args.mt_nmads, "high"
    )
    if args.max_pct_mt is not None:
        adata.obs["qc_high_mito"] |= adata.obs["pct_counts_mt"] > args.max_pct_mt

    reason_columns = [
        "qc_low_counts",
        "qc_high_counts",
        "qc_low_genes",
        "qc_high_genes",
        "qc_high_mito",
    ]
    adata.obs["passes_qc"] = ~adata.obs[reason_columns].any(axis=1)
    adata.obs["qc_reason"] = adata.obs[reason_columns].apply(
        lambda row: ";".join(column for column in reason_columns if bool(row[column])) or "pass",
        axis=1,
    )
    thresholds = threshold_table(
        adata, args.sample_key, args.nmads, args.mt_nmads, args.max_pct_mt
    )
    thresholds_path = Path(args.thresholds)
    thresholds_path.parent.mkdir(parents=True, exist_ok=True)
    thresholds.to_csv(thresholds_path, index=False)

    figure_dir = Path(args.figdir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    plot_frame = adata.obs[
        [args.sample_key, "total_counts", "n_genes_by_counts", "pct_counts_mt", "passes_qc"]
    ].copy()
    for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        grid = sns.displot(
            data=plot_frame,
            x=metric,
            hue="passes_qc",
            col=args.sample_key,
            col_wrap=4,
            bins=60,
            facet_kws={"sharex": False, "sharey": False},
        )
        threshold_index = thresholds.set_index(args.sample_key)
        line_columns = {
            "total_counts": ("min_total_counts", "max_total_counts"),
            "n_genes_by_counts": ("min_genes", "max_genes"),
            "pct_counts_mt": (None, "max_pct_mt"),
        }
        for ax, sample in zip(grid.axes.flat, grid.col_names):
            low_column, high_column = line_columns[metric]
            row = threshold_index.loc[str(sample)]
            if low_column:
                ax.axvline(row[low_column], color="#a61b1b", linestyle="--", linewidth=1)
            ax.axvline(row[high_column], color="#a61b1b", linestyle="--", linewidth=1)
        grid.fig.savefig(figure_dir / f"01_{metric}_by_sample.png", dpi=150, bbox_inches="tight")
        plt.close(grid.fig)

    summary = (
        adata.obs.groupby(args.sample_key, observed=True)
        .agg(
            input_cells=("passes_qc", "size"),
            passing_cells=("passes_qc", "sum"),
            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),
            median_pct_mt=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    summary["failed_cells"] = summary["input_cells"] - summary["passing_cells"]
    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)

    record_provenance(
        adata,
        "quality_control",
        {
            "sample_key": args.sample_key,
            "counts_layer": args.counts_layer,
            "nmads": args.nmads,
            "mt_nmads": args.mt_nmads,
            "max_pct_mt": args.max_pct_mt,
            "mt_prefixes": prefixes,
            "filter_failed": args.filter_failed,
            "input_cells": int(adata.n_obs),
            "passing_cells": int(adata.obs["passes_qc"].sum()),
            "thresholds_file": str(thresholds_path),
        },
        ["anndata", "scanpy", "numpy", "pandas", "scipy", "seaborn"],
    )

    if args.filter_failed:
        adata = adata[adata.obs["passes_qc"]].copy()
    if args.min_cells_per_gene > 1:
        active_counts = adata.layers[args.counts_layer]
        detected = np.asarray((active_counts > 0).sum(axis=0)).ravel()
        adata = adata[:, detected >= args.min_cells_per_gene].copy()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(summary.to_string(index=False))
    print(f"Saved: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Flag low-quality cells per sample")
    parser.add_argument("--input", required=True, help="Raw-count .h5ad")
    parser.add_argument("--output", default="results/01_qc.h5ad")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--nmads", type=float, default=3.0)
    parser.add_argument("--mt-nmads", type=float, default=3.0)
    parser.add_argument(
        "--max-pct-mt",
        type=float,
        default=None,
        help="Optional dataset-justified absolute mitochondrial percentage ceiling",
    )
    parser.add_argument("--mt-prefixes", default="MT-,mt-")
    parser.add_argument("--min-cells-per-gene", type=int, default=1)
    parser.add_argument("--filter-failed", action="store_true")
    parser.add_argument("--summary", default="results/01_qc_summary.csv")
    parser.add_argument("--thresholds", default="results/01_qc_thresholds.csv")
    parser.add_argument("--figdir", default="figures")
    main(parser.parse_args())
