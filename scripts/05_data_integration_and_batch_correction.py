"""Optional batch correction with Harmony or scVI and auditable outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from workflow_utils import get_counts, record_provenance, require_columns


def _prepare_features(adata, batch_key: str, counts_layer: str, n_top_genes: int, seed: int):
    counts = get_counts(adata, counts_layer)
    if "log1p_norm" not in adata.layers:
        normalized = sc.pp.normalize_total(adata, layer=counts_layer, target_sum=1e4, inplace=False)
        adata.layers["log1p_norm"] = sc.pp.log1p(normalized["X"], copy=True)

    sc.pp.highly_variable_genes(
        adata,
        layer=counts_layer,
        n_top_genes=min(n_top_genes, adata.n_vars),
        flavor="seurat_v3",
        batch_key=batch_key,
        subset=False,
    )
    work = adata[:, adata.var["highly_variable"]].copy()
    work.X = work.layers["log1p_norm"].copy()
    n_pcs = min(50, work.n_obs - 1, work.n_vars - 1)
    if n_pcs < 2:
        raise ValueError("At least three cells and three selected genes are required.")
    sc.tl.pca(work, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    return work, counts, n_pcs


def _run_harmony(work, batch_key: str, seed: int):
    try:
        import harmonypy  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("Harmony requires harmonypy; install environment.yml.") from exc
    sc.external.pp.harmony_integrate(
        work,
        key=batch_key,
        basis="X_pca",
        adjusted_basis="X_harmony",
        random_state=seed,
    )
    return work.obsm["X_harmony"], "X_harmony"


def _run_scvi(work, batch_key: str, counts_layer: str, epochs: int, seed: int):
    try:
        import scvi
    except ImportError as exc:
        raise RuntimeError("scVI requires environment-advanced.yml.") from exc
    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(work, layer=counts_layer, batch_key=batch_key)
    model = scvi.model.SCVI(work, gene_likelihood="nb")
    model.train(max_epochs=epochs, early_stopping=True)
    return model.get_latent_representation(), "X_scVI", model


def main(args: argparse.Namespace) -> None:
    adata = sc.read_h5ad(args.input)
    require_columns(adata, [args.batch_key])
    if adata.obs[args.batch_key].isna().any():
        raise ValueError(f"obs[{args.batch_key!r}] contains missing batch assignments.")
    n_batches = adata.obs[args.batch_key].nunique()
    if n_batches < 2 and args.method != "none":
        raise ValueError("Batch correction requires at least two batches; use --method none.")

    work, _, n_pcs = _prepare_features(
        adata, args.batch_key, args.counts_layer, args.n_top_genes, args.seed
    )
    sc.pp.neighbors(work, use_rep="X_pca", key_added="uncorrected_neighbors")
    sc.tl.umap(work, neighbors_key="uncorrected_neighbors", random_state=args.seed)
    uncorrected_umap = work.obsm["X_umap"].copy()

    selected = "harmony" if args.method == "auto" else args.method
    model = None
    if selected == "harmony":
        latent, representation = _run_harmony(work, args.batch_key, args.seed)
    elif selected == "scvi":
        latent, representation, model = _run_scvi(
            work, args.batch_key, args.counts_layer, args.epochs, args.seed
        )
    elif selected == "none":
        latent, representation = work.obsm["X_pca"], "X_pca"
    else:
        raise ValueError(f"Unknown method: {selected}")

    work.obsm["X_integrated"] = np.asarray(latent)
    sc.pp.neighbors(work, use_rep="X_integrated", key_added="integrated_neighbors")
    sc.tl.umap(work, neighbors_key="integrated_neighbors", random_state=args.seed)

    adata.obsm["X_pca_uncorrected"] = np.asarray(work.obsm["X_pca"])
    adata.obsm["X_umap_uncorrected"] = uncorrected_umap
    adata.obsm["X_integrated"] = np.asarray(work.obsm["X_integrated"])
    adata.obsm["X_umap_integrated"] = np.asarray(work.obsm["X_umap"])
    adata.obsm["X_umap"] = np.asarray(work.obsm["X_umap"])
    connectivities_key = work.uns["integrated_neighbors"]["connectivities_key"]
    distances_key = work.uns["integrated_neighbors"]["distances_key"]
    adata.uns["neighbors"] = dict(work.uns["integrated_neighbors"])
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    adata.uns["neighbors"]["distances_key"] = "distances"
    adata.obsp["connectivities"] = work.obsp[connectivities_key]
    adata.obsp["distances"] = work.obsp[distances_key]

    figure = Path(args.figure)
    figure.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    categories = adata.obs[args.batch_key].astype("category")
    colors = categories.cat.codes.to_numpy()
    for ax, coordinates, title in [
        (axes[0], adata.obsm["X_umap_uncorrected"], "Uncorrected"),
        (axes[1], adata.obsm["X_umap_integrated"], f"Corrected: {selected}"),
    ]:
        ax.scatter(coordinates[:, 0], coordinates[:, 1], c=colors, s=3, cmap="tab20", alpha=0.7)
        ax.set_title(title)
        ax.set(xticks=[], yticks=[])
    fig.suptitle(f"Batch comparison: {args.batch_key}")
    fig.tight_layout()
    fig.savefig(figure, dpi=150)
    plt.close(fig)

    if model is not None and args.model_dir:
        model_path = Path(args.model_dir)
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model.save(model_path, overwrite=True)

    record_provenance(
        adata,
        "batch_correction",
        {
            "method_requested": args.method,
            "method_selected": selected,
            "batch_key": args.batch_key,
            "n_batches": int(n_batches),
            "counts_layer": args.counts_layer,
            "n_top_genes": args.n_top_genes,
            "n_pcs": n_pcs,
            "representation": representation,
            "epochs": args.epochs if selected == "scvi" else None,
            "seed": args.seed,
            "rationale": args.rationale,
        },
        ["anndata", "scanpy", "numpy", "harmonypy", "scvi-tools"],
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"Batch representation: {selected} ({latent.shape[1]} dimensions)")
    print(f"Saved: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct an explicitly defined batch effect")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results/05_integrated.h5ad")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--method", choices=["auto", "harmony", "scvi", "none"], default="auto")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--rationale",
        default="Batch correction requested by the analyst after reviewing batch-associated structure.",
    )
    parser.add_argument("--figure", default="figures/05_batch_comparison.png")
    parser.add_argument("--model-dir", default="models/scvi")
    main(parser.parse_args())
