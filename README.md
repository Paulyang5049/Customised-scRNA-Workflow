# Customised scRNA Workflow

When dealing with single-cell RNA-seq (scRNA-seq) data, starting at the raw count matrix can be overwhelming. Why are there so many steps just to understand what cells we have? This repository is my customized workflow, structured sequentially like a student's personal notes, to explain exactly *why* we do what we do at each stage of the analysis.

## My Analysis Notes & Workflow

I've structured the analytical code into 18 numbered scripts in the `scripts/` directory. Each script is production-ready with `argparse` CLI, function wrapping, and file-based I/O (read `.h5ad` → process → write `.h5ad`). Here is my thought process as I go through the pipeline:

### 1. Cleaning up the mess (Prep and QC)
Biology is messy, and sequencing captures that mess. 
* **`01_quality_control`**: First, I need to filter out dead or dying cells (which have too much mitochondrial RNA) and doublets (when two cells get caught in one droplet). If I don't remove ambient RNA floating around, it will confuse my gene expression later.
* **`02_normalization`**: Not all droplets get sequenced equally deeply. Some cells have 10,000 counts, others only 2,000. I need to normalize these counts so I'm comparing apples to apples across cells.
* **`03_feature_selection`**: Human cells have ~20,000 genes, but most of them are just housekeeping genes doing the same thing everywhere. I run this to find the "highly variable genes" (HVGs) that tell me what makes one cell different from another.

### 2. Giving the data shape (Structure and Typing)
* **`04_dimensionality_reduction`**: Humans can't visualize 2,000-dimensional space (the HVGs). I use PCA to compress the math, and then UMAP/t-SNE to squash it down to 2D so I can actually look at the data on my screen.
* **`05_data_integration_and_batch_correction`**: If I processed samples on different days, or from different donors, the technical noise might overshadow the biology. This step mathematically aligns the datasets so a T-cell from Donor A plots next to a T-cell from Donor B.
* **`06_clustering`**: Now that the cells are plotted properly, I use graph-based algorithms to draw boundaries around them, grouping similar cells together.
* **`07_cell_type_annotation`**: A cluster is just a number (e.g., "Cluster 3"). I look at the marker genes expressed in that cluster to figure out its biological identity (e.g., "Oh, Cluster 3 is a Macrophage!").

### 3. Asking the real biological questions (Downstream Analysis)
Once I know who the cells are, I can figure out what they are doing.
* **`08_differential_gene_expression`**: What genes are turned on in Disease vs. Control within a specific cell type?
* **`09_gene_set_enrichment_analysis`**: Instead of looking at single genes, are entire biological pathways (like "inflammatory response") activated?
* **`10_compositional_analysis`**: Did the *proportion* of a cell type change? (e.g., "Are there more macrophages in the treated sample?")
* **`11_pseudotime_and_trajectories` & `12_rna_velocity`**: Cells aren't static. These steps help me map out how stem cells differentiate into mature cells, or how cells transition from state A to state B over time.

### 4. Advanced Modalities & Mechanisms
Beyond standard scRNA-seq, biology works in complex networks.
* **`13_perturbation_modeling`**: How do specific external changes (like drug treatments or CRISPR knockouts) alter the cell state?
* **`14_cell_cell_communication`**: Cells talk to each other. By looking at receptor-ligand pairs, I can infer which groups are signaling one another.
* **`15_gene_regulatory_networks`**: Which transcription factors are pulling the strings behind the scenes?
* **`16_paired_integration` & `17_advanced_integration`**: Working with multi-omics (like CITE-seq or ATAC-seq + RNA-seq in the same cell) to get a multi-layered view of the biology.
* **`18_lineage_tracing`**: Following the exact genealogical family tree of cells alongside their transcriptomes.

## Setup and Environments
All required packages (`scanpy`, `anndata`, `scvi-tools`, etc.) are covered in the environment files. You can just replicate what I did using:
```bash
# Clone and set up the conda environment
conda env create -f environment.yml
conda activate sc-tutorial-env
```

---
*Attribution: The core code and tutorials in this pipeline are adapted from the [Theis Lab Single-Cell Best Practices](https://github.com/theislab/single-cell-best-practices) open-source repository and the accompanying Nature Reviews Genetics publication [Heumos et al., 2023].*