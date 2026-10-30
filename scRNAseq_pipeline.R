## =============================================================================
## Script:   scRNAseq_pipeline.R
## Title:    Reproducible scRNA-seq analysis pipeline (single-file master script)
## =============================================================================
## A general-purpose, reproducible single-cell RNA-seq pipeline covering the
## full standard Bioconductor workflow: QC -> normalisation -> feature
## selection -> PCA/UMAP -> (optional) batch correction -> clustering ->
## marker genes -> (optional) cell type annotation -> (optional) pseudobulk
## differential expression -> summary report.
##
## It follows the structure taught in the GenE2 scRNA-seq introduction
## workshop (jserizay.com/GenE2-scRNAseq-intro) - SingleCellExperiment,
## scran/scuttle, scater, SingleR, DESeq2 - generalised here so it runs on
## ANY dataset and ANY metadata: every path, column name, and parameter is
## read from one YAML config file (config/pipeline_config.yml), so adapting
## the pipeline to a new dataset never means editing this script.
##
## USAGE
## -----
##   Rscript scRNAseq_pipeline.R [path/to/config.yml] [--generate-toy]
##
##   - config path defaults to "config/pipeline_config.yml" if omitted.
##   - pass --generate-toy as the second argument to first simulate a small
##     synthetic test dataset (via scuttle::mockSCE()) and write it to
##     data/raw/ in every format this script understands, so you can test
##     the whole pipeline end-to-end with no data of your own. Use it with
##     config/toy_dataset_config.yml, e.g.:
##       Rscript scRNAseq_pipeline.R config/toy_dataset_config.yml --generate-toy
##
## REPRODUCIBILITY
## ----------------
##   - Every parameter lives in the YAML config - nothing is hard-coded here.
##   - A single seed (project$seed) controls every stochastic step (pooled
##     size-factor estimation, UMAP, graph-based clustering).
##   - The Bioconductor release is pinned (BIOC_VERSION below), so every
##     collaborator installs the same package versions.
##   - Every step is logged with a timestamp to logs/pipeline.log.
##   - sessionInfo() and a copy of the config used are written to the output
##     folder at the end of the run, so results are always traceable back to
##     the exact parameters and package versions that produced them.
##   - For a fully version-locked environment, initialise renv in this
##     project (renv::init(); after a successful run, renv::snapshot()) -
##     see the end of this file for details.
##
## OPTIONAL STEPS
## ---------------
##   Batch correction, cell type annotation, and differential expression are
##   OFF by default and each checks, at runtime, that the metadata column it
##   needs actually exists before running - so pointing this script at a
##   dataset that lacks a "batch" or "condition" column never crashes it;
##   the step simply logs that it was skipped. Turn them on in the config
##   once your metadata has the columns they need.
## =============================================================================


## =============================================================================
## SECTION 0A - Helper functions (logging, config, safe package install)
## =============================================================================

#' Timestamped logging to console and a log file
log_msg <- function(..., log_file = LOG_FILE) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  message(msg)
  if (!is.null(log_file)) {
    ensure_dir(dirname(log_file))
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  invisible(msg)
}

#' Create a directory (recursively) if it does not already exist
ensure_dir <- function(path) {
  if (!is.null(path) && nzchar(path) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Install (if needed) and attach CRAN packages
install_and_load_cran <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(TRUE)
}

#' Install (if needed) and attach Bioconductor packages
install_and_load_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(TRUE)
}

#' Resolve a colData() column by name, returning FALSE instead of erroring if
#' it is absent - every OPTIONAL metadata column (batch, condition, ...) is
#' checked this way so the pipeline degrades gracefully on datasets that
#' don't have that column, rather than crashing.
has_coldata_column <- function(sce, colname) {
  !is.null(colname) && nzchar(colname) && colname %in% colnames(SummarizedExperiment::colData(sce))
}

#' Simulate a small synthetic test dataset in every input format this script
#' understands (see --generate-toy in the usage notes above).
generate_toy_dataset <- function() {
  install_and_load_bioc(c("SingleCellExperiment", "scuttle", "DropletUtils"))
  install_and_load_cran(c("Matrix"))
  set.seed(1234)

  log_msg("Simulating a toy SingleCellExperiment with scuttle::mockSCE() ...")
  sce <- scuttle::mockSCE(ncells = 600, ngenes = 2000)

  rownames(sce)[1:13] <- paste0("MT-GENE", 1:13)
  colData(sce)$sample_id <- paste0("cell_", seq_len(ncol(sce)))
  colData(sce)$batch     <- sample(c("batch1", "batch2"), ncol(sce), replace = TRUE)
  colData(sce)$condition <- sample(c("control", "treatment"), ncol(sce), replace = TRUE)

  out_dir <- file.path("data", "raw")
  ensure_dir(out_dir)

  DropletUtils::write10xCounts(
    path = file.path(out_dir, "10x_mtx"), x = counts(sce),
    gene.id = rownames(sce), gene.symbol = rownames(sce),
    barcodes = colnames(sce), overwrite = TRUE
  )
  saveRDS(sce, file.path(out_dir, "toy_sce.rds"))
  saveRDS(as.matrix(counts(sce)), file.path(out_dir, "toy_matrix.rds"))
  write.csv(as.data.frame(as.matrix(counts(sce))), file.path(out_dir, "toy_counts.csv"))

  ensure_dir("data")
  write.csv(as.data.frame(colData(sce)), file.path("data", "sample_metadata.csv"), row.names = FALSE)

  log_msg("Toy dataset written to ", out_dir, " (10x_mtx, rds_sce, rds_matrix, csv) ",
          "plus data/sample_metadata.csv. Use config/toy_dataset_config.yml to run ",
          "the full pipeline on it (input_type: 'rds_sce').")
}


## =============================================================================
## SECTION 0B - Resolve arguments, read config, start logging
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1 && !startsWith(args[[1]], "--")) args[[1]] else "config/pipeline_config.yml"
GENERATE_TOY <- "--generate-toy" %in% args

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml", repos = "https://cloud.r-project.org")
}
if (!file.exists(config_path)) {
  stop("Config file not found at: ", config_path,
       "\nCopy config/pipeline_config.yml, edit it for your dataset, and pass its path.")
}
cfg <- yaml::read_yaml(config_path)

ensure_dir(cfg$paths$output_dir)
ensure_dir(cfg$paths$log_dir)
LOG_FILE <- file.path(cfg$paths$log_dir, "pipeline.log")
FIG_DIR  <- file.path(cfg$paths$output_dir, "figures")
ensure_dir(FIG_DIR)

log_msg("=================================================================")
log_msg("Starting pipeline run: '", cfg$project$name, "' (config: ", config_path, ")")
log_msg("=================================================================")

RUN_START <- Sys.time()

if (GENERATE_TOY) {
  log_msg("--generate-toy flag set: simulating test data before the run.")
  generate_toy_dataset()
}


## =============================================================================
## STEP 00 - Environment setup: install/load every package, pin Bioconductor,
##           set the global random seed
## =============================================================================
log_msg("STEP 00: environment setup")

BIOC_VERSION <- "3.19"  # pinned release -> same package versions for every collaborator
BiocManager::install(version = BIOC_VERSION, update = FALSE, ask = FALSE)

install_and_load_cran(c("yaml", "ggplot2", "patchwork", "igraph", "Matrix"))
install_and_load_bioc(c(
  "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
  "DropletUtils",   # 10x mtx/h5 readers
  "scuttle",         # QC metrics, normalisation, pseudobulk aggregation
  "scran",           # HVG modelling, pooled size factors, marker detection
  "scater",          # PCA/UMAP/t-SNE runners, QC & dim-reduction plots
  "bluster",         # graph-based clustering (NNGraphParam)
  "batchelor",       # fastMNN batch correction
  "SingleR",         # reference-based cell type annotation
  "celldex",         # curated reference expression datasets for SingleR
  "DESeq2"           # pseudobulk differential expression
))

set.seed(cfg$project$seed)
log_msg("Random seed set to ", cfg$project$seed, ". Packages ready.")


## =============================================================================
## STEP 01 - Load raw counts (any supported format) into a SingleCellExperiment
## =============================================================================
## Supports input_type: "10x_mtx", "10x_h5", "csv", "rds_sce", "rds_matrix".
## Optional external metadata (paths$metadata_path) is attached by matching
## paths$metadata_join_column, or by row order if that column isn't found.
log_msg("STEP 01: loading data (input_type = '", cfg$paths$input_type, "') from ",
        cfg$paths$input_path)

input_path <- cfg$paths$input_path
input_type <- cfg$paths$input_type

sce <- switch(
  input_type,
  "10x_mtx" = DropletUtils::read10xCounts(input_path, col.names = TRUE),
  "10x_h5"  = DropletUtils::read10xCounts(input_path, col.names = TRUE),
  "csv" = {
    mat <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))
    SingleCellExperiment(assays = list(counts = Matrix::Matrix(mat, sparse = TRUE)))
  },
  "rds_sce" = {
    obj <- readRDS(input_path)
    if (!is(obj, "SingleCellExperiment")) {
      stop("input_type is 'rds_sce' but the object in ", input_path,
           " is not a SingleCellExperiment (class: ", class(obj)[1], ").")
    }
    obj
  },
  "rds_matrix" = {
    mat <- readRDS(input_path)
    SingleCellExperiment(assays = list(counts = Matrix::Matrix(as.matrix(mat), sparse = TRUE)))
  },
  stop("Unknown paths$input_type '", input_type, "'. Must be one of: ",
       "10x_mtx, 10x_h5, csv, rds_sce, rds_matrix.")
)

rownames(sce) <- scuttle::uniquifyFeatureNames(
  ID = if (!is.null(rowData(sce)$ID)) rowData(sce)$ID else rownames(sce),
  names = if (!is.null(rowData(sce)$Symbol)) rowData(sce)$Symbol else rownames(sce)
)
colnames(sce) <- make.unique(colnames(sce))

meta_path <- cfg$paths$metadata_path
if (!is.null(meta_path) && nzchar(meta_path) && file.exists(meta_path)) {
  meta <- read.csv(meta_path, check.names = FALSE)
  join_col <- cfg$paths$metadata_join_column

  if (!is.null(join_col) && nzchar(join_col) &&
      join_col %in% colnames(meta) && join_col %in% colnames(colData(sce))) {
    log_msg("Joining metadata on column '", join_col, "'")
    matched <- meta[match(colData(sce)[[join_col]], meta[[join_col]]), , drop = FALSE]
    for (cn in setdiff(colnames(matched), colnames(colData(sce)))) {
      colData(sce)[[cn]] <- matched[[cn]]
    }
  } else if (nrow(meta) == ncol(sce)) {
    log_msg("Metadata join column not found on both sides; attaching by row order.")
    for (cn in colnames(meta)) colData(sce)[[cn]] <- meta[[cn]]
  } else {
    warning("Could not attach metadata_path: no usable join column and row ",
            "counts differ. Continuing without external metadata.")
  }
} else {
  log_msg("No metadata_path provided or file not found - continuing without it.")
}

log_msg("Loaded: ", nrow(sce), " genes x ", ncol(sce), " cells.")


## =============================================================================
## STEP 02 - Quality control: per-cell QC metrics, outlier-based filtering,
##           diagnostic plots
## =============================================================================
## Tools: scuttle::perCellQCMetrics(), scuttle::isOutlier(), scater::plotColData()
## Mitochondrial genes are flagged by regex (qc$mito_pattern); cells are
## flagged as outliers adaptively (median +/- nmads MADs) rather than with
## hard-coded cutoffs, so this works across datasets with very different
## sequencing depth.
log_msg("STEP 02: quality control")

mito_idx <- grep(cfg$qc$mito_pattern, rownames(sce))
log_msg(length(mito_idx), " mitochondrial gene(s) matched by pattern '",
        cfg$qc$mito_pattern, "'")

qc_metrics <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = mito_idx))
colData(sce) <- cbind(colData(sce), qc_metrics)

low_lib    <- scuttle::isOutlier(qc_metrics$sum, log = TRUE, type = "lower", nmads = cfg$qc$nmads)
low_genes  <- scuttle::isOutlier(qc_metrics$detected, log = TRUE, type = "lower", nmads = cfg$qc$nmads)
high_mito  <- scuttle::isOutlier(qc_metrics$subsets_Mito_percent, type = "higher", nmads = cfg$qc$nmads)
hard_floor <- qc_metrics$detected < cfg$qc$min_genes_per_cell
discard    <- low_lib | low_genes | high_mito | hard_floor
colData(sce)$discard <- discard

log_msg("QC filtering: low lib size = ", sum(low_lib), ", low genes = ", sum(low_genes),
        ", high mito = ", sum(high_mito), ", below floor = ", sum(hard_floor),
        ", TOTAL discarded = ", sum(discard), " / ", ncol(sce))

p1 <- scater::plotColData(sce, y = "sum", colour_by = "discard") +
  ggplot2::scale_y_log10() + ggplot2::ggtitle("Library size")
p2 <- scater::plotColData(sce, y = "detected", colour_by = "discard") +
  ggplot2::scale_y_log10() + ggplot2::ggtitle("Genes detected")
p3 <- scater::plotColData(sce, y = "subsets_Mito_percent", colour_by = "discard") +
  ggplot2::ggtitle("% mitochondrial counts")
ggplot2::ggsave(file.path(FIG_DIR, "02_qc_library_size.png"), p1, width = 5, height = 4, dpi = 150)
ggplot2::ggsave(file.path(FIG_DIR, "02_qc_genes_detected.png"), p2, width = 5, height = 4, dpi = 150)
ggplot2::ggsave(file.path(FIG_DIR, "02_qc_percent_mito.png"), p3, width = 5, height = 4, dpi = 150)
write.csv(as.data.frame(colData(sce)), file.path(cfg$paths$output_dir, "qc_metrics.csv"), row.names = FALSE)

sce <- sce[, !discard]
genes_detected_in <- Matrix::rowSums(counts(sce) > 0)
sce <- sce[genes_detected_in >= cfg$qc$min_cells_per_gene, ]
log_msg("After QC: ", nrow(sce), " genes x ", ncol(sce), " cells.")


## =============================================================================
## STEP 03 - Normalisation: size factors + log-transform
## =============================================================================
## Tools: scran::quickCluster(), scran::computeSumFactors(), scuttle::logNormCounts()
## method "deconvolution": pools cells into quick low-resolution clusters
## first so size factors are estimated within homogeneous groups before
## deconvolution (Lun, Bach & Marioni 2016) - recommended for heterogeneous
## populations. "libsize": simple library-size normalisation (faster, fine
## for more homogeneous datasets).
log_msg("STEP 03: normalisation (method = '", cfg$normalization$method, "')")

set.seed(cfg$project$seed)
if (identical(cfg$normalization$method, "deconvolution")) {
  quick_clusters <- scran::quickCluster(sce)
  log_msg(length(unique(quick_clusters)), " quick-cluster(s) found for size-factor pooling.")
  sce <- scran::computeSumFactors(sce, clusters = quick_clusters, min.mean = 0.1)
} else if (identical(cfg$normalization$method, "libsize")) {
  sce <- scuttle::computeLibraryFactors(sce)
} else {
  stop("Unknown normalization$method '", cfg$normalization$method, "'. Must be 'deconvolution' or 'libsize'.")
}
log_msg("Size factors: min = ", round(min(sizeFactors(sce)), 3),
        ", median = ", round(median(sizeFactors(sce)), 3),
        ", max = ", round(max(sizeFactors(sce)), 3))

sce <- scuttle::logNormCounts(sce)
log_msg("'logcounts' assay added.")


## =============================================================================
## STEP 04 - Feature selection: model gene-level variance, pick top HVGs
## =============================================================================
## Tools: scran::modelGeneVar(), scran::getTopHVGs()
## Fits a trend to the mean-variance relationship of log-expression and
## decomposes each gene's variance into technical (explained by the trend)
## and biological (residual) components; genes with the highest biological
## variance are kept for PCA.
log_msg("STEP 04: feature selection (top ", cfg$feature_selection$n_hvgs, " HVGs)")

gene_var <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(gene_var, n = cfg$feature_selection$n_hvgs)
metadata(sce)$hvgs <- hvgs
log_msg("Selected ", length(hvgs), " highly variable genes.")
write.csv(as.data.frame(gene_var), file.path(cfg$paths$output_dir, "gene_variance.csv"))

fit_df <- data.frame(mean = gene_var$mean, total = gene_var$total)
curve_df <- data.frame(mean = sort(gene_var$mean), trend = metadata(gene_var)$trend(sort(gene_var$mean)))
p <- ggplot2::ggplot(fit_df, ggplot2::aes(x = mean, y = total)) +
  ggplot2::geom_point(alpha = 0.4, size = 0.8) +
  ggplot2::geom_line(data = curve_df, ggplot2::aes(x = mean, y = trend), colour = "dodgerblue", linewidth = 1) +
  ggplot2::labs(x = "Mean log-expression", y = "Variance of log-expression",
                title = "Mean-variance trend and fitted technical noise") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(FIG_DIR, "04_mean_variance_trend.png"), p, width = 6, height = 4.5, dpi = 150)


## =============================================================================
## STEP 05 - Dimensionality reduction: PCA on HVGs, then UMAP / t-SNE
## =============================================================================
## Tools: scater::runPCA(), scater::runUMAP(), scater::runTSNE()
## UMAP/t-SNE are run on the PCA embedding, not the raw expression matrix -
## faster and less noisy, and the standard two-stage approach in Bioconductor.
log_msg("STEP 05: dimensionality reduction (", cfg$dimensionality_reduction$n_pcs, " PCs)")

set.seed(cfg$project$seed)
sce <- scater::runPCA(sce, subset_row = metadata(sce)$hvgs, ncomponents = cfg$dimensionality_reduction$n_pcs)

pct_var <- attr(reducedDim(sce, "PCA"), "percentVar")
if (!is.null(pct_var)) {
  p <- ggplot2::ggplot(data.frame(PC = seq_along(pct_var), variance = pct_var),
                        ggplot2::aes(x = PC, y = variance)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Principal component", y = "% variance explained", title = "PCA variance explained") +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(FIG_DIR, "05_pca_variance_explained.png"), p, width = 6, height = 4, dpi = 150)
}

if (isTRUE(cfg$dimensionality_reduction$run_umap)) {
  sce <- scater::runUMAP(sce, dimred = "PCA")
  log_msg("UMAP computed.")
}
if (isTRUE(cfg$dimensionality_reduction$run_tsne)) {
  sce <- scater::runTSNE(sce, dimred = "PCA")
  log_msg("t-SNE computed.")
}


## =============================================================================
## STEP 06 - (Optional) Batch correction across samples/runs
## =============================================================================
## Tools: batchelor::fastMNN(), scater::runUMAP()
## Only runs if batch_correction$run is TRUE AND the batch_var column exists.
## Otherwise the PCA embedding is simply carried forward unchanged.
log_msg("STEP 06: batch correction (run = ", isTRUE(cfg$batch_correction$run), ")")

batch_var <- cfg$batch_correction$batch_var
if (isTRUE(cfg$batch_correction$run) && has_coldata_column(sce, batch_var)) {
  log_msg("Running fastMNN on '", batch_var, "' (",
          length(unique(colData(sce)[[batch_var]])), " batches).")
  set.seed(cfg$project$seed)
  mnn_out <- batchelor::fastMNN(sce, batch = colData(sce)[[batch_var]], subset.row = metadata(sce)$hvgs)
  reducedDim(sce, "MNN") <- reducedDim(mnn_out, "corrected")
  USE_DIMRED <- "MNN"
  sce <- scater::runUMAP(sce, dimred = "MNN", name = "UMAP")

  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = batch_var) +
    ggplot2::ggtitle("UMAP after fastMNN batch correction, coloured by batch")
  ggplot2::ggsave(file.path(FIG_DIR, "06_umap_by_batch.png"), p, width = 6, height = 5, dpi = 150)
  log_msg("Batch correction complete.")
} else {
  if (isTRUE(cfg$batch_correction$run)) {
    log_msg("batch_correction$run is TRUE but column '", batch_var, "' not found - skipping.")
  } else {
    log_msg("Batch correction disabled - PCA embedding passed through unchanged.")
  }
  USE_DIMRED <- "PCA"
}


## =============================================================================
## STEP 07 - Graph-based clustering
## =============================================================================
## Tools: bluster::clusterRows() with bluster::NNGraphParam() (igraph backend)
## Builds a shared/k-NN graph over cells in the chosen embedding (PCA or MNN),
## partitioned with Walktrap or Louvain community detection - deterministic
## given the same graph and seed, so results are reproducible.
cluster_fun <- switch(
  cfg$clustering$method,
  "graph_walktrap" = "walktrap",
  "graph_louvain"  = "louvain",
  stop("Unknown clustering$method '", cfg$clustering$method, "'. Must be 'graph_walktrap' or 'graph_louvain'.")
)
log_msg("STEP 07: clustering on '", USE_DIMRED, "' (", cluster_fun, ", k = ", cfg$clustering$k_neighbors, ")")

set.seed(cfg$project$seed)
clusters <- bluster::clusterRows(
  reducedDim(sce, USE_DIMRED),
  bluster::NNGraphParam(k = cfg$clustering$k_neighbors, cluster.fun = cluster_fun)
)
colData(sce)$cluster <- factor(clusters)
log_msg(nlevels(colData(sce)$cluster), " cluster(s) found.")

cluster_sizes <- as.data.frame(table(cluster = colData(sce)$cluster))
colnames(cluster_sizes) <- c("cluster", "n_cells")
write.csv(cluster_sizes, file.path(cfg$paths$output_dir, "cluster_sizes.csv"), row.names = FALSE)

if ("UMAP" %in% reducedDimNames(sce)) {
  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "cluster") + ggplot2::ggtitle("UMAP by cluster")
  ggplot2::ggsave(file.path(FIG_DIR, "07_umap_by_cluster.png"), p, width = 6, height = 5, dpi = 150)
}


## =============================================================================
## STEP 08 - Marker gene detection
## =============================================================================
## Tools: scran::scoreMarkers(), scater::plotDots()
## Scores every gene, per cluster, on complementary effect-size statistics
## (summary logFC, AUC, Cohen's d) from pairwise comparisons against every
## other cluster; mean.AUC gives a threshold-free ranking of how well each
## gene distinguishes its cluster from the rest.
log_msg("STEP 08: marker gene detection")

N_TOP_MARKERS <- 5
marker_list <- scran::scoreMarkers(sce, groups = colData(sce)$cluster)

top_markers_combined <- do.call(rbind, lapply(names(marker_list), function(cl) {
  df <- as.data.frame(marker_list[[cl]])
  df$gene <- rownames(df)
  df <- df[order(-df$mean.AUC), ]
  df$cluster <- cl
  head(df[, c("cluster", "gene", "mean.AUC", "mean.logFC.cohen")], N_TOP_MARKERS)
}))
rownames(top_markers_combined) <- NULL

for (cl in names(marker_list)) {
  write.csv(as.data.frame(marker_list[[cl]]),
            file.path(cfg$paths$output_dir, paste0("markers_cluster_", cl, ".csv")))
}
write.csv(top_markers_combined, file.path(cfg$paths$output_dir, "top_markers_combined.csv"), row.names = FALSE)
log_msg("Wrote per-cluster marker tables and top-", N_TOP_MARKERS, "-per-cluster summary.")

top_genes <- unique(top_markers_combined$gene)
if (length(top_genes) > 0) {
  p <- scater::plotDots(sce, features = top_genes, group = "cluster") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(file.path(FIG_DIR, "08_marker_dotplot.png"), p,
                   width = max(6, length(top_genes) * 0.3), height = 5, dpi = 150, limitsize = FALSE)
}


## =============================================================================
## STEP 09 - (Optional) Reference-based cell type annotation
## =============================================================================
## Tools: celldex::<reference>(), SingleR::SingleR()
## Only runs if annotation$run is TRUE. Requires internet access on first use
## (the reference is downloaded and cached by celldex/ExperimentHub).
log_msg("STEP 09: cell type annotation (run = ", isTRUE(cfg$annotation$run), ")")

if (isTRUE(cfg$annotation$run)) {
  ref_fun_name <- cfg$annotation$reference
  if (!exists(ref_fun_name, where = asNamespace("celldex"), inherits = FALSE)) {
    stop("annotation$reference '", ref_fun_name, "' is not a function in celldex. ",
         "See ?celldex for available references.")
  }
  ref_fun <- get(ref_fun_name, envir = asNamespace("celldex"))
  log_msg("Downloading/loading reference '", ref_fun_name, "'...")
  ref <- ref_fun()

  label_col <- cfg$annotation$label_column
  if (!label_col %in% colnames(colData(ref))) {
    stop("annotation$label_column '", label_col, "' not found in the reference. ",
         "Available: ", paste(colnames(colData(ref)), collapse = ", "))
  }

  pred <- SingleR::SingleR(test = sce, ref = ref, labels = colData(ref)[[label_col]])
  colData(sce)$cell_type <- pred$labels
  write.csv(as.data.frame(pred[, c("labels", "pruned.labels")]),
            file.path(cfg$paths$output_dir, "annotation_scores.csv"))
  log_msg("Annotation complete, ", length(unique(pred$labels)), " label(s) assigned.")

  if ("UMAP" %in% reducedDimNames(sce)) {
    p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "cell_type") +
      ggplot2::ggtitle("UMAP by SingleR cell type")
    ggplot2::ggsave(file.path(FIG_DIR, "09_umap_by_celltype.png"), p, width = 7, height = 5, dpi = 150)
  }
} else {
  log_msg("Annotation disabled in config - colData(sce)$cell_type set to NA. ",
          "Set a species/tissue-appropriate celldex reference and annotation$run: true to enable.")
  colData(sce)$cell_type <- NA_character_
}


## =============================================================================
## STEP 10 - (Optional) Pseudobulk differential expression between conditions
## =============================================================================
## Tools: scuttle::aggregateAcrossCells(), DESeq2
## Collapses counts to one profile per cluster x condition x sample before
## applying standard bulk DE (DESeq2) - avoids treating cells from the same
## sample as independent replicates, which inflates false discovery rates
## (Squair et al. 2021, Nat Commun). Only runs if differential_expression$run
## is TRUE and condition_var has both levels in group_levels present.
condition_var <- cfg$differential_expression$condition_var
group_levels  <- cfg$differential_expression$group_levels
run_de <- isTRUE(cfg$differential_expression$run) &&
  has_coldata_column(sce, condition_var) &&
  all(group_levels %in% unique(colData(sce)[[condition_var]]))

log_msg("STEP 10: pseudobulk differential expression (run = ", run_de, ")")

if (run_de) {
  sample_col <- if (has_coldata_column(sce, "sample_id")) "sample_id" else condition_var
  log_msg("Using '", sample_col, "' as the pseudobulk replicate unit.")

  ids <- S4Vectors::DataFrame(cluster = colData(sce)$cluster,
                               condition = colData(sce)[[condition_var]],
                               sample = colData(sce)[[sample_col]])
  pb <- scuttle::aggregateAcrossCells(sce, ids = ids)
  clusters_present <- unique(colData(pb)$cluster)
  log_msg("Running DESeq2 for ", length(clusters_present), " cluster(s): '",
          group_levels[2], "' vs '", group_levels[1], "'.")

  for (cl in clusters_present) {
    pb_cl <- pb[, colData(pb)$cluster == cl & colData(pb)$condition %in% group_levels]
    if (length(unique(colData(pb_cl)$condition)) < 2 || ncol(pb_cl) < 2) {
      log_msg("  cluster ", cl, ": skipped (not enough replicates per group).")
      next
    }
    col_data <- data.frame(condition = factor(colData(pb_cl)$condition, levels = group_levels))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts(pb_cl)),
                                           colData = col_data, design = ~condition)
    dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("condition", group_levels[2], group_levels[1]))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[order(res_df$padj), ]
    out_path <- file.path(cfg$paths$output_dir, paste0("de_results_cluster_", cl, ".csv"))
    write.csv(res_df, out_path, row.names = FALSE)
    log_msg("  cluster ", cl, ": ", sum(res_df$padj < 0.05, na.rm = TRUE), " gene(s) at padj < 0.05.")
  }
} else if (isTRUE(cfg$differential_expression$run)) {
  log_msg("differential_expression$run is TRUE but column '", condition_var,
          "' is missing or lacks both levels ", paste(group_levels, collapse = " / "), " - skipping.")
} else {
  log_msg("Differential expression disabled in config - skipping.")
}


## =============================================================================
## STEP 11 - Final summary figure panel + markdown run report
## =============================================================================
## Combines the key UMAP views (cluster, cell type, condition, QC) into one
## figure, and writes a short report summarising every parameter and output -
## a sanity check and a methods-section draft in one.
log_msg("STEP 11: final visualization and report")

panels <- list()
if ("UMAP" %in% reducedDimNames(sce)) {
  panels$cluster <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "cluster") + ggplot2::ggtitle("Cluster")
  if (!all(is.na(colData(sce)$cell_type))) {
    panels$cell_type <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "cell_type") +
      ggplot2::ggtitle("Cell type (SingleR)")
  }
  if (has_coldata_column(sce, condition_var)) {
    panels$condition <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = condition_var) +
      ggplot2::ggtitle("Condition")
  }
  panels$qc <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = "sum") + ggplot2::ggtitle("Total UMI count")
}

if (length(panels) > 0) {
  combined <- Reduce(`+`, panels) + patchwork::plot_layout(ncol = 2)
  ggplot2::ggsave(file.path(FIG_DIR, "11_summary_panel.png"), combined,
                   width = 11, height = 4 * ceiling(length(panels) / 2), dpi = 150, limitsize = FALSE)
}

saveRDS(sce, file.path(cfg$paths$output_dir, "final_sce.rds"))

report_path <- file.path(cfg$paths$output_dir, "pipeline_summary.md")
n_clusters <- nlevels(colData(sce)$cluster)
report_lines <- c(
  paste0("# ", cfg$project$name, " - pipeline run summary"),
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Random seed: ", cfg$project$seed),
  "",
  "## Result",
  paste0("- Final object: ", nrow(sce), " genes x ", ncol(sce), " cells"),
  paste0("- Clusters found: ", n_clusters),
  paste0("- Batch correction: ", if (isTRUE(cfg$batch_correction$run)) "yes (fastMNN)" else "not run"),
  paste0("- Annotation: ", if (isTRUE(cfg$annotation$run)) paste0("yes (", cfg$annotation$reference, ")") else "not run"),
  paste0("- Differential expression: ", if (run_de) "yes (pseudobulk DESeq2)" else "not run"),
  "",
  "## Parameters used",
  paste0("- QC: mito_pattern = '", cfg$qc$mito_pattern, "', nmads = ", cfg$qc$nmads,
         ", min_genes_per_cell = ", cfg$qc$min_genes_per_cell,
         ", min_cells_per_gene = ", cfg$qc$min_cells_per_gene),
  paste0("- Normalisation: ", cfg$normalization$method),
  paste0("- HVGs: ", cfg$feature_selection$n_hvgs),
  paste0("- PCA components: ", cfg$dimensionality_reduction$n_pcs),
  paste0("- Clustering: ", cfg$clustering$method, ", k = ", cfg$clustering$k_neighbors),
  "",
  "## Output files",
  "See output/ for final_sce.rds, *.csv tables (QC metrics, gene variance, ",
  "cluster sizes, marker genes, DE results), output/figures/ for every plot, ",
  "logs/pipeline.log for the full run log, and logs/sessionInfo.txt for exact ",
  "package versions."
)
writeLines(report_lines, report_path)
log_msg("Run report written -> ", report_path)


## =============================================================================
## Reproducibility artefacts: copy the config used, capture sessionInfo()
## =============================================================================
file.copy(config_path, file.path(cfg$paths$output_dir, basename(config_path)), overwrite = TRUE)
if (isTRUE(cfg$reproducibility$save_session_info)) {
  writeLines(capture.output(sessionInfo()), file.path(cfg$paths$log_dir, "sessionInfo.txt"))
  log_msg("sessionInfo() saved.")
}

log_msg("=================================================================")
log_msg("Pipeline finished in ", round(difftime(Sys.time(), RUN_START, units = "mins"), 2),
        " min. See ", cfg$paths$output_dir, " for results.")
log_msg("=================================================================")

## -----------------------------------------------------------------------------
## Stronger reproducibility with renv (recommended for real projects)
## -----------------------------------------------------------------------------
## This script pins the Bioconductor RELEASE, one lever for consistent package
## versions across machines. For a fully locked, bit-for-bit reproducible
## environment that you can archive or hand to a collaborator, run once in
## this project folder (not part of this script - a manual, one-time step):
##   install.packages("renv"); renv::init()
##   # ... run this script successfully ...
##   renv::snapshot()            # records the exact version of every package used
## Anyone else can then reproduce your exact environment with renv::restore().
## =============================================================================
