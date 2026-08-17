# User guide

This guide covers the supported core path. Scripts 09–18 are optional research modules and should not be treated as a single mandatory sequence.

## 1. System requirements

- Linux or macOS. Windows users should use WSL2.
- Conda, Miniforge, Mambaforge, or Micromamba.
- At least 16 GB RAM for a small-to-medium dataset. Large atlases may need substantially more.
- Enough disk space for several `.h5ad` checkpoints.
- R is installed by the Conda environment; a separate system R is not required.

The core environment supports validation, QC, SoupX, scDblFinder, normalization, feature selection, PCA/UMAP, Harmony, clustering, and reporting. The advanced environment additionally includes scVI, velocity, automated annotation, perturbation, and communication packages.

## 2. Install with Conda

Clone the repository, then run:

```bash
conda env create -f environment.yml
conda activate scrna-workflow
python scripts/00_validate_input.py --help
```

Mamba is faster and can use the same file:

```bash
mamba env create -f environment.yml
conda activate scrna-workflow
```

For scripts that require scVI, scVelo, CellTypist, LIANA, decoupler, or pertpy:

```bash
conda env create -f environment-advanced.yml
conda activate scrna-workflow-advanced
```

To update an existing environment after the file changes:

```bash
conda env update --name scrna-workflow --file environment.yml --prune
```

## 3. Required input contract

The primary input is an AnnData `.h5ad` object with:

- cells/barcodes in rows and genes in columns;
- unique `obs_names` and `var_names`;
- raw, non-negative, integer-like UMI counts in `layers["counts"]` or `.X`;
- `obs["sample"]` identifying the biological sample or capture channel;
- `obs["batch"]` only when a technical batch is known;
- `obs["condition"]` for downstream condition comparisons;
- `obs["donor"]` when multiple samples belong to a donor.

Do not use the condition column as the batch key. A covariate should only be corrected when it represents unwanted technical variation rather than the biology being tested.

Ambient-RNA correction additionally requires an unfiltered droplet matrix. Convert the filtered and raw matrices into two AnnData objects with matching gene identifiers and cell barcodes.

## 4. Validate before analysis

```bash
python scripts/00_validate_input.py \
  --input data/raw_counts.h5ad \
  --sample-key sample \
  --required-obs sample,batch,condition \
  --report results/00_validation.json
```

Validation exits with status 2 if counts or required metadata are invalid.

## 5. QC and threshold review

The first pass flags cells without deleting them. Thresholds are calculated within sample from counts, detected genes, and mitochondrial percentage.

```bash
python scripts/01_quality_control.py \
  --input data/raw_counts.h5ad \
  --output results/01_qc_review.h5ad \
  --sample-key sample

python scripts/19_generate_report.py \
  --input results/01_qc_review.h5ad \
  --output reports/qc_review.html
```

Review the plots and sample summaries. If the thresholds are appropriate, create the filtered checkpoint explicitly:

The exact per-sample lower and upper limits are saved in `results/01_qc_thresholds.csv` and drawn as dashed lines on the QC plots.

```bash
python scripts/01_quality_control.py \
  --input data/raw_counts.h5ad \
  --output results/01_qc_filtered.h5ad \
  --sample-key sample \
  --filter-failed
```

There is no universal mitochondrial cutoff. Use `--max-pct-mt` only when the observed distributions and tissue biology justify it.

## 6. Ambient RNA

SoupX is optional but recommended for droplet data when the unfiltered matrix is available:

```bash
python scripts/01c_ambient_rna_correction.py \
  --input results/01_qc_filtered.h5ad \
  --raw-droplets data/raw_unfiltered_droplets.h5ad \
  --output results/01c_soupx.h5ad \
  --set-active-counts
```

The original counts are retained in `layers["counts_pre_ambient"]`; corrected counts are stored in `layers["soupx_counts"]` and, with `--set-active-counts`, become `layers["counts"]` for downstream analysis.

## 7. Doublet detection

Run scDblFinder by capture channel or sample:

```bash
python scripts/01b_doublet_detection.py \
  --input results/01c_soupx.h5ad \
  --output results/01b_doublet_review.h5ad \
  --sample-key sample
```

This adds `scDblFinder_score`, `scDblFinder_class`, `scDblFinder_partition`, and `passes_doublet_qc`. Review doublet locations and sample rates before filtering. To create an explicit singlet-only checkpoint, rerun with `--filter-doublets`.

If `scDblFinder` or R fails, the script exits with an error. It will not silently substitute a different method.

## 8. Normalization through clustering

Choose the input checkpoint appropriate to your review:

```bash
python scripts/02_normalization.py \
  --input results/01b_doublet_review.h5ad \
  --output results/02_normalized.h5ad

python scripts/03_feature_selection.py
python scripts/04_dimensionality_reduction.py
```

Batch correction is optional. Compare the uncorrected and corrected plots and ensure condition-specific biology is not erased.

```bash
python scripts/05_data_integration_and_batch_correction.py \
  --input results/04_dimred.h5ad \
  --output results/05_integrated.h5ad \
  --batch-key batch \
  --method harmony \
  --rationale "Library-preparation day separates otherwise matched donors"

python scripts/06_clustering.py \
  --input results/05_integrated.h5ad \
  --output results/06_clustered.h5ad
```

Use `--method none` when correction is unwarranted. Use `--method scvi` only in the advanced environment and when dataset size and batch complexity justify it.

## 9. Automatic report

Generate a self-contained report from any checkpoint:

```bash
python scripts/19_generate_report.py \
  --input results/06_clustered.h5ad \
  --output reports/index.html \
  --summary-json reports/summary.json \
  --figures-dir figures
```

The report includes QC distributions, available UMAP views, cell counts, workflow provenance, package versions, and warnings for missing stages. CSV summaries are written to `reports/tables/`.

## 10. Reproducibility checklist

- Keep raw counts immutable.
- Save every major checkpoint rather than overwriting inputs.
- Record why thresholds, batch correction, references, and exclusions were chosen.
- Preserve donor/sample identifiers for replicate-aware differential expression.
- Use the same Conda file and random seed when reproducing a run.
- Archive the HTML report, JSON summary, figures, and environment file with results.

## Troubleshooting

`Rscript is unavailable`: activate `scrna-workflow` and check `which Rscript`.

`Counts are not integer-like`: the input contains normalized expression. Recover the original UMI matrix before proceeding.

`Missing sample column`: add a capture-channel identifier to `.obs`; do not combine independent runs for doublet detection.

Out-of-memory errors: keep matrices sparse, analyze capture channels separately, and avoid scripts that explicitly convert matrices to dense data frames.
