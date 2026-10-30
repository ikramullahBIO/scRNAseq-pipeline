# scRNA-seq analysis pipeline (R / Bioconductor) — single-file version

A general-purpose, reproducible single-cell RNA-seq pipeline, delivered as
**one master script**: `scRNAseq_pipeline.R`. It covers the full standard
workflow — QC → normalisation → feature selection → PCA/UMAP → (optional)
batch correction → clustering → marker genes → (optional) cell-type
annotation → (optional) pseudobulk differential expression → summary report.

It follows the same Bioconductor ecosystem used in the [GenE2 scRNA-seq
workshop](https://jserizay.com/GenE2-scRNAseq-intro/devel/pages/workshop.html)
(`SingleCellExperiment`, `scran`/`scuttle`, `scater`, `SingleR`, `DESeq2`),
generalised so it runs on **any** dataset and **any** metadata — every input
path, column name, and parameter is read from one YAML config file, so
adapting it to a new dataset never requires editing the script.

## Files

```
scRNAseq-pipeline/
├── README.md
├── .gitignore
├── scRNAseq_pipeline.R           # the entire pipeline, one file, 11 titled steps
└── config/
    ├── pipeline_config.yml       # edit this for your own dataset
    └── toy_dataset_config.yml    # matches the synthetic test dataset (see below)
```

Everything the script creates (`data/`, `output/`, `logs/`) is generated at
run time and is not part of the repo — see `.gitignore`.

## Requirements

- R >= 4.3
- Internet access the first time you run it (to install packages, and if
  `annotation: run: true`, to download the SingleR reference)
- No manual package installation needed — the script installs everything it
  needs (CRAN + Bioconductor) the first time it runs.

## Quick start

**1. Try it on a synthetic dataset first (no data of your own needed):**

```bash
Rscript scRNAseq_pipeline.R config/toy_dataset_config.yml --generate-toy
```

The `--generate-toy` flag simulates ~600 cells / 2,000 genes with
`scuttle::mockSCE()`, writes them to `data/raw/`, and then runs the complete
pipeline — including batch correction and differential expression — end to
end in under a couple of minutes. Check `output_toy/pipeline_summary.md`
when it finishes.

**2. Run it on your own data:**

1. Copy `config/pipeline_config.yml` (e.g. to `config/my_dataset.yml`).
2. Edit `paths: input_path` / `input_type` to point at your data — a 10x mtx
   directory, a 10x `.h5` file, a plain counts CSV, or an `.rds` file (exact
   format details are in the comments of the config file).
3. If you have per-cell/per-sample metadata (condition, batch, patient,
   timepoint, ...), set `paths: metadata_path`.
4. Adjust `qc: mito_pattern` for your species (`"^MT-"` human, `"^mt-"`
   mouse — check your reference's gene-naming convention).
5. Turn on `batch_correction: run`, `annotation: run`, and/or
   `differential_expression: run` as needed, and set their variable/
   reference names to match your metadata.
6. Run it:

```bash
Rscript scRNAseq_pipeline.R config/my_dataset.yml
```

## What the script does, step by step

| Step | What it does | Key tools |
|---|---|---|
| 00 | Installs/loads every package, pins the Bioconductor release, sets the random seed | `BiocManager` |
| 01 | Reads raw counts in any supported format into a `SingleCellExperiment`; attaches external metadata | `DropletUtils`, `scuttle` |
| 02 | Computes per-cell QC metrics, flags outlier cells adaptively (MAD-based), filters cells and low-detection genes | `scuttle::perCellQCMetrics`, `scuttle::isOutlier`, `scater::plotColData` |
| 03 | Pooled deconvolution (or library-size) size factors, log-normalisation | `scran::quickCluster`, `scran::computeSumFactors`, `scuttle::logNormCounts` |
| 04 | Models mean-variance trend, selects top highly variable genes (HVGs) | `scran::modelGeneVar`, `scran::getTopHVGs` |
| 05 | PCA on HVGs, then UMAP/t-SNE on the PCA space | `scater::runPCA`, `scater::runUMAP`, `scater::runTSNE` |
| 06 | Optional fastMNN batch correction across batches; recomputes UMAP on the corrected space | `batchelor::fastMNN` |
| 07 | Graph-based clustering (Walktrap or Louvain) on the chosen embedding | `bluster::clusterRows`, `bluster::NNGraphParam` |
| 08 | Scores marker genes per cluster, writes tables, plots top markers | `scran::scoreMarkers`, `scater::plotDots` |
| 09 | Optional reference-based cell-type annotation | `celldex`, `SingleR::SingleR` |
| 10 | Optional pseudobulk differential expression between two condition levels, per cluster | `scuttle::aggregateAcrossCells`, `DESeq2` |
| 11 | Combined summary figure + a markdown run report | `scater`, `patchwork` |

Batch correction, annotation, and differential expression are **off by
default**, and each checks at runtime that the metadata column it needs
actually exists before running — so pointing the script at a dataset that
lacks a `batch` or `condition` column never crashes it; that step simply
logs that it was skipped.

## Reproducibility

- **One config, one run.** Every parameter lives in the YAML config — none
  of it is hard-coded in the script.
- **One seed controls every stochastic step** (`project: seed`): pooled
  size-factor estimation, UMAP, and graph-based clustering.
- **Pinned Bioconductor release** (`BIOC_VERSION` near the top of the
  script) so everyone running it installs the same package versions.
- **Full run log** (`logs/pipeline.log`): every step is timestamped, with
  the parameters used and cell/gene/cluster counts at each stage.
- **`sessionInfo()` is captured** at the end of every run
  (`logs/sessionInfo.txt`) — exact R and package versions, for a methods
  section or for debugging a discrepancy later.
- **The config used is copied into `output/`** alongside the results, so a
  results folder is always self-describing.
- For even stronger reproducibility, initialise
  [`renv`](https://rstudio.github.io/renv/) in this project
  (`renv::init()`, then `renv::snapshot()` after a successful run) to lock
  every package to an exact version — see the comment block at the end of
  `scRNAseq_pipeline.R`.

## Attribution

The analysis logic (QC → normalisation → HVGs → PCA/UMAP → batch
correction → clustering → markers → annotation → pseudobulk DE) follows the
structure taught in the GenE2 scRNA-seq introduction workshop
(jserizay.com/GenE2-scRNAseq-intro), generalised here into a single,
reusable, config-driven script.
