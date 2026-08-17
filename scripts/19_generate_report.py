"""Generate a portable HTML audit report from a workflow AnnData object."""

from __future__ import annotations

import argparse
import base64
import html
import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from workflow_utils import package_versions, write_json

QC_METRICS = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "scDblFinder_score",
    "ambient_contamination_fraction",
]


def _figure_data_uri(fig) -> str:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _qc_figure(adata) -> str | None:
    metrics = [metric for metric in QC_METRICS if metric in adata.obs]
    if not metrics:
        return None
    fig, axes = plt.subplots(len(metrics), 1, figsize=(8, 3 * len(metrics)), squeeze=False)
    for ax, metric in zip(axes[:, 0], metrics):
        values = pd.to_numeric(adata.obs[metric], errors="coerce").dropna()
        ax.hist(values, bins=60, color="#3478bf", alpha=0.85)
        ax.set(title=metric, ylabel="Cells")
    axes[-1, 0].set_xlabel("Value")
    fig.tight_layout()
    return _figure_data_uri(fig)


def _umap_figure(adata, color_key: str) -> str | None:
    basis = "X_umap_integrated" if "X_umap_integrated" in adata.obsm else "X_umap"
    if basis not in adata.obsm or color_key not in adata.obs:
        return None
    coordinates = np.asarray(adata.obsm[basis])
    values = adata.obs[color_key].astype(str)
    counts = values.value_counts()
    if len(counts) > 20:
        keep = set(counts.head(19).index)
        values = values.where(values.isin(keep), "Other")
    categories = pd.Categorical(values)

    fig, ax = plt.subplots(figsize=(7, 6))
    scatter = ax.scatter(
        coordinates[:, 0], coordinates[:, 1], c=categories.codes, cmap="tab20", s=3, alpha=0.75
    )
    del scatter
    handles = []
    for index, category in enumerate(categories.categories):
        handles.append(
            plt.Line2D(
                [], [], marker="o", linestyle="", markersize=5, color=plt.cm.tab20(index % 20),
                label=str(category)
            )
        )
    if len(handles) <= 20:
        ax.legend(handles=handles, frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set(title=f"UMAP colored by {color_key}", xticks=[], yticks=[])
    fig.tight_layout()
    return _figure_data_uri(fig)


def _existing_figures(directory: Path, limit: int = 24) -> list[tuple[str, str]]:
    found = []
    if not directory.is_dir():
        return found
    for path in sorted(directory.glob("*.png"))[:limit]:
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        found.append((path.name, f"data:image/png;base64,{encoded}"))
    return found


def _audit_flags(adata) -> list[str]:
    flags = []
    if "counts" not in adata.layers:
        flags.append("No immutable raw-count layer named 'counts' was found.")
    if "scDblFinder_class" not in adata.obs:
        flags.append("Doublet detection has not been recorded.")
    if "ambient_contamination_fraction" not in adata.obs:
        flags.append("Ambient-RNA correction or assessment has not been recorded.")
    if "X_integrated" not in adata.obsm:
        flags.append("No batch-corrected latent representation is present (this may be intentional).")
    if not any(key in adata.obs for key in ["cell_type", "cell_type_coarse", "celltypist_pred"]):
        flags.append("No cell-type annotation column was found.")
    if "scrna_workflow" not in adata.uns:
        flags.append("Workflow provenance is absent.")
    return flags


def main(args: argparse.Namespace) -> None:
    adata = sc.read_h5ad(args.input)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "input": str(Path(args.input).resolve()),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "layers": list(adata.layers.keys()),
        "embeddings": list(adata.obsm.keys()),
        "obs_columns": adata.obs.columns.tolist(),
        "audit_flags": _audit_flags(adata),
        "versions": package_versions(["anndata", "scanpy", "numpy", "pandas", "scipy"]),
        "provenance": adata.uns.get("scrna_workflow", {}),
    }
    write_json(args.summary_json, summary)

    table_dir = output.parent / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    available_metrics = [metric for metric in QC_METRICS if metric in adata.obs]
    if available_metrics:
        adata.obs[available_metrics].describe().T.to_csv(table_dir / "qc_metric_summary.csv")
    for key in [args.sample_key, args.batch_key, args.cell_type_key, "scDblFinder_class"]:
        if key in adata.obs:
            adata.obs[key].value_counts(dropna=False).rename("n_cells").to_csv(
                table_dir / f"counts_{key}.csv"
            )

    images: list[tuple[str, str]] = []
    qc_image = _qc_figure(adata)
    if qc_image:
        images.append(("QC distributions", qc_image))
    for key in [args.batch_key, args.cell_type_key, "scDblFinder_class", "passes_qc"]:
        image = _umap_figure(adata, key)
        if image:
            images.append((f"UMAP: {key}", image))
    images.extend(_existing_figures(Path(args.figures_dir)))

    flags_html = "".join(f"<li>{html.escape(flag)}</li>" for flag in summary["audit_flags"])
    if not flags_html:
        flags_html = "<li>No automatic audit warnings.</li>"
    images_html = "".join(
        f"<section><h3>{html.escape(title)}</h3><img src=\"{uri}\" alt=\"{html.escape(title)}\"></section>"
        for title, uri in images
    )
    provenance = html.escape(json.dumps(summary["provenance"], indent=2, default=str))
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>scRNA-seq workflow report</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;color:#1d2733}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem}}
.card,section{{border:1px solid #d9e0e7;border-radius:10px;padding:1rem;margin:1rem 0;background:#fff}}
.warning{{background:#fff8e6;border-left:5px solid #e6a700}} img{{max-width:100%;height:auto}}
pre{{overflow:auto;background:#f5f7f9;padding:1rem}} small{{color:#52606d}}
</style></head><body>
<h1>scRNA-seq workflow report</h1><small>{html.escape(summary['input'])}</small>
<div class="cards"><div class="card"><strong>{summary['n_cells']:,}</strong><br>cells</div>
<div class="card"><strong>{summary['n_genes']:,}</strong><br>genes</div>
<div class="card"><strong>{len(summary['layers'])}</strong><br>layers</div>
<div class="card"><strong>{len(summary['embeddings'])}</strong><br>embeddings</div></div>
<section class="warning"><h2>Automatic audit</h2><ul>{flags_html}</ul></section>
<section><h2>Available data</h2><p><strong>Layers:</strong> {html.escape(', '.join(summary['layers']) or 'none')}</p>
<p><strong>Embeddings:</strong> {html.escape(', '.join(summary['embeddings']) or 'none')}</p></section>
<h2>Figures</h2>{images_html or '<p>No reportable figures were found.</p>'}
<section><h2>Provenance</h2><pre>{provenance}</pre></section>
<p><small>Generated from the saved AnnData object. Automatic warnings require analyst review.</small></p>
</body></html>"""
    output.write_text(document)
    print(f"HTML report: {output}")
    print(f"JSON summary: {args.summary_json}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an HTML audit report")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="reports/index.html")
    parser.add_argument("--summary-json", default="reports/summary.json")
    parser.add_argument("--figures-dir", default="figures")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--cell-type-key", default="cell_type")
    main(parser.parse_args())

