#!/usr/bin/env Rscript

# ==============================================================================
# PhiSpace Deconvolution Pipeline for Benchmarking
# ==============================================================================
# This script runs PhiSpace on prepared benchmark datasets and outputs
# cell type proportions in the standard format used by the benchmarking framework.
#
# Supported normalization methods:
#   - log1p: log(1+x) transformation
#   - libsize: Library size normalization (normalize to median total counts)
#   - scran: scran pooling-based size factor normalization
#   - sctransform: Seurat's SCTransform (variance stabilizing transformation)
#   - spanorm: SpaNorm spatial-aware normalization (for spatial data)
#
# Usage:
#   Rscript 02_run_phispace.R <dataset_dir> <output_dir> [ref_norm] [query_norm]
#
# Arguments:
#   dataset_dir: Directory containing prepared dataset files
#   output_dir: Directory for output results
#   ref_norm: Normalization for reference scRNA-seq (default: sctransform)
#   query_norm: Normalization for query spatial (default: log1p)
# ==============================================================================

suppressPackageStartupMessages({
  library(PhiSpace)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# Increase future globals size limit for SCTransform with large datasets
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB

# ==============================================================================
# Normalization Functions
# ==============================================================================

normalize_log1p <- function(sce, assay_name = "log1p") {
  # Simple log(1+x) transformation
  cat("    Applying log1p normalization...\n")

  counts_mat <- assay(sce, "counts")
  log1p_mat <- log1p(counts_mat)

  assay(sce, assay_name) <- log1p_mat
  return(sce)
}

normalize_libsize <- function(sce, assay_name = "libsize") {
  # Library size normalization: normalize to median library size
  cat("    Applying library size normalization...\n")

  counts_mat <- assay(sce, "counts")
  lib_sizes <- colSums(counts_mat)
  median_lib <- median(lib_sizes)

  # Normalize to median library size, then log1p
  scale_factors <- median_lib / lib_sizes
  scale_factors[!is.finite(scale_factors)] <- 0

  norm_mat <- t(t(counts_mat) * scale_factors)
  norm_mat <- log1p(norm_mat)

  assay(sce, assay_name) <- norm_mat
  return(sce)
}

normalize_scran <- function(sce, assay_name = "scran") {
  # scran pooling-based normalization
  cat("    Applying scran normalization...\n")

  # Load scran if available
  if (!requireNamespace("scran", quietly = TRUE)) {
    stop("scran package is required for scran normalization. Install with: BiocManager::install('scran')")
  }
  if (!requireNamespace("scuttle", quietly = TRUE)) {
    stop("scuttle package is required for scran normalization. Install with: BiocManager::install('scuttle')")
  }

  library(scran)
  library(scuttle)

  # Quick clustering for size factor estimation
  counts_mat <- assay(sce, "counts")

  # For large datasets, use a subset for clustering
  if (ncol(sce) > 5000) {
    set.seed(42)
    subset_idx <- sample(ncol(sce), 5000)
    quick_clust <- quickCluster(counts_mat[, subset_idx])
    # Extend clusters to full dataset (assign to nearest)
    clusters <- rep(1, ncol(sce))
    clusters[subset_idx] <- as.integer(quick_clust)
  } else {
    clusters <- quickCluster(counts_mat)
  }

  # Compute size factors
  sce <- computeSumFactors(sce, clusters = clusters)

  # Handle zero or negative size factors
  sf <- sizeFactors(sce)
  sf[sf <= 0] <- min(sf[sf > 0])
  sizeFactors(sce) <- sf

  # Normalize and log transform
  sce <- logNormCounts(sce, log = TRUE, transform = "log1p")

  # Rename assay to our standard name
  assay(sce, assay_name) <- assay(sce, "logcounts")

  return(sce)
}

normalize_sctransform <- function(counts_mat, col_data, assay_name = "sctransform") {
  # Seurat's SCTransform (variance stabilizing transformation)
  # Returns SCE with SCT normalized data
  cat("    Applying SCTransform normalization...\n")

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = as(counts_mat, "dgCMatrix"),
    meta.data = as.data.frame(col_data)
  )

  # Apply SCTransform
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

  # Extract normalized data
  sct_data <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")

  # Create SCE with SCT data
  sce <- SingleCellExperiment(
    assays = list(sctransform = sct_data),
    colData = col_data
  )

  return(sce)
}

normalize_spanorm <- function(counts_mat, col_data, assay_name = "spanorm") {
  # SpaNorm spatial-aware normalization
  # Returns SCE with SpaNorm normalized data
  cat("    Applying SpaNorm normalization...\n")

  # Load SpaNorm
  if (!requireNamespace("SpaNorm", quietly = TRUE)) {
    stop("SpaNorm package is required. Install with: BiocManager::install('SpaNorm')")
  }
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment package is required. Install with: BiocManager::install('SpatialExperiment')")
  }

  library(SpaNorm)
  library(SpatialExperiment)
  library(matrixStats)

  # Filter out zero-sum genes first
  gene_sums <- rowSums(counts_mat)
  counts_mat <- counts_mat[gene_sums > 0, ]
  cat(sprintf("    After removing zero genes: %d genes\n", nrow(counts_mat)))

  # Create spatial coordinates (grid layout for simulated data without coordinates)
  n_spots <- ncol(counts_mat)
  grid_size <- ceiling(sqrt(n_spots))

  # Create grid coordinates with scaling (SpaNorm expects pixel coordinates)
  x_coords <- rep(seq_len(grid_size), each = grid_size)[seq_len(n_spots)]
  y_coords <- rep(seq_len(grid_size), grid_size)[seq_len(n_spots)]
  coords <- cbind(
    pxl_col_in_fullres = x_coords * 100,
    pxl_row_in_fullres = y_coords * 100
  )
  rownames(coords) <- colnames(counts_mat)

  cat(sprintf("    Created grid coordinates: %d x %d\n", grid_size, grid_size))

  # Create SpatialExperiment object
  spe <- SpatialExperiment(
    assays = list(counts = as(counts_mat, "dgCMatrix")),
    colData = col_data,
    spatialCoords = coords
  )

  # Filter genes expressed in at least 10% of spots
  keep <- filterGenes(spe, prop = 0.10)
  spe <- spe[keep, ]

  # Cap at 5000 genes for computational efficiency
  # Keep most variable genes if exceeds cap
  MAX_GENES <- 5000
  if (nrow(spe) > MAX_GENES) {
    cat(sprintf("    Reducing from %d to %d most variable genes for efficiency\n", nrow(spe), MAX_GENES))
    gene_vars <- rowVars(as.matrix(assay(spe, "counts")))
    top_genes <- order(gene_vars, decreasing = TRUE)[1:MAX_GENES]
    spe <- spe[top_genes, ]
  }
  cat(sprintf("    After gene filtering: %d genes\n", nrow(spe)))

  # Compute size factors first (required for SpaNorm to work properly)
  spe <- fastSizeFactors(spe)

  # Apply SpaNorm with default settings and 1-hour timeout
  SPANORM_TIMEOUT <- 3600  # 1 hour in seconds

  spanorm_result <- tryCatch({
    # Set CPU time limit
    setTimeLimit(cpu = SPANORM_TIMEOUT, elapsed = SPANORM_TIMEOUT, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)

    SpaNorm(
      spe,
      sample.p = 0.1,
      df.tps = 2,
      verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    SpaNorm timeout or error after %d seconds: %s\n", SPANORM_TIMEOUT, e$message))
    return(NULL)
  })

  # Check if SpaNorm succeeded
  if (is.null(spanorm_result)) {
    cat("    Falling back to simple log1p normalization due to SpaNorm timeout\n")
    # Fallback to simple log1p normalization
    norm_data <- log1p(as.matrix(assay(spe, "counts")))
  } else {
    spe <- spanorm_result
    # Extract normalized data (SpaNorm stores in logcounts)
    norm_data <- assay(spe, "logcounts")
  }

  # Create SCE with normalized data
  sce <- SingleCellExperiment(
    assays = list(spanorm = norm_data),
    colData = col_data
  )

  return(sce)
}

apply_normalization <- function(counts_mat, col_data, method, is_reference = TRUE) {
  # Apply the specified normalization method
  # Returns SCE with normalized assay

  data_type <- ifelse(is_reference, "reference", "query")
  cat(sprintf("  Normalizing %s with method: %s\n", data_type, method))

  if (method == "sctransform") {
    # SCTransform is special - creates its own SCE
    sce <- normalize_sctransform(counts_mat, col_data, assay_name = method)
  } else if (method == "spanorm") {
    # SpaNorm is special - needs spatial coordinates
    sce <- normalize_spanorm(counts_mat, col_data, assay_name = method)
  } else {
    # Other methods start with counts SCE
    sce <- SingleCellExperiment(
      assays = list(counts = as(counts_mat, "dgCMatrix")),
      colData = col_data
    )

    # Remove zero-count features
    sce <- zeroFeatQC(sce)

    if (method == "log1p") {
      sce <- normalize_log1p(sce, assay_name = method)
    } else if (method == "libsize") {
      sce <- normalize_libsize(sce, assay_name = method)
    } else if (method == "scran") {
      sce <- normalize_scran(sce, assay_name = method)
    } else {
      stop(sprintf("Unknown normalization method: %s", method))
    }
  }

  cat(sprintf("    Completed. Features: %d, Samples: %d\n", nrow(sce), ncol(sce)))
  return(sce)
}

# ==============================================================================
# Parse command line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 02_run_phispace.R <dataset_dir> <output_dir> [ref_norm] [query_norm]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  dataset_dir: Directory containing prepared dataset files\n")
  cat("  output_dir: Directory for output results\n")
  cat("  ref_norm: Normalization for reference scRNA-seq (default: sctransform)\n")
  cat("  query_norm: Normalization for query spatial (default: log1p)\n")
  cat("\n")
  cat("Normalization options:\n")
  cat("  log1p      - Simple log(1+x) transformation\n")
  cat("  libsize    - Library size normalization (normalize to median, then log1p)\n")
  cat("  scran      - scran pooling-based size factor normalization\n")
  cat("  sctransform - Seurat's SCTransform (variance stabilizing transformation)\n")
  cat("  spanorm    - SpaNorm spatial-aware normalization (recommended for spatial query)\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript 02_run_phispace.R data/dataset1 results/dataset1 sctransform log1p\n")
  quit(status = 1)
}

dataset_dir <- args[1]
output_dir <- args[2]
ref_norm <- ifelse(length(args) >= 3, args[3], "sctransform")
query_norm <- ifelse(length(args) >= 4, args[4], "log1p")
n_hvg <- ifelse(length(args) >= 5, as.integer(args[5]), 0L)  # 0 = use all genes

# Validate normalization methods
valid_methods <- c("log1p", "libsize", "scran", "sctransform", "spanorm")
if (!ref_norm %in% valid_methods) {
  stop(sprintf("Invalid reference normalization method: %s. Must be one of: %s",
               ref_norm, paste(valid_methods, collapse = ", ")))
}
if (!query_norm %in% valid_methods) {
  stop(sprintf("Invalid query normalization method: %s. Must be one of: %s",
               query_norm, paste(valid_methods, collapse = ", ")))
}

cat(paste(rep("=", 60), collapse=""), "\n")
cat("PhiSpace Benchmarking Pipeline\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Reference normalization:", ref_norm, "\n")
cat("Query normalization:", query_norm, "\n")
cat("HVG selection:", ifelse(n_hvg > 0, paste0(n_hvg, " genes"), "all genes"), "\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load data
# ==============================================================================

cat("\nLoading data...\n")

# Load scRNA-seq counts
scrna_counts_file <- file.path(dataset_dir, "scRNA_counts.csv")
if (!file.exists(scrna_counts_file)) {
  stop("scRNA counts file not found: ", scrna_counts_file)
}
scrna_counts <- fread(scrna_counts_file, header = TRUE)
gene_names <- scrna_counts[[1]]
scrna_counts <- as.matrix(scrna_counts[, -1, with = FALSE])
rownames(scrna_counts) <- gene_names
cat("  scRNA-seq counts loaded:", nrow(scrna_counts), "genes x", ncol(scrna_counts), "cells\n")

# Load cell type annotations
celltype_file <- file.path(dataset_dir, "scRNA_celltype.csv")
if (!file.exists(celltype_file)) {
  stop("Cell type file not found: ", celltype_file)
}
celltype_df <- fread(celltype_file)
celltypes <- celltype_df$celltype
names(celltypes) <- as.character(celltype_df$cell)
cat("  Cell types loaded:", length(unique(celltypes)), "unique types\n")

# Filter cell types with < 25 cells (match RCTD filtering criterion)
celltype_counts <- table(celltypes)
rare_types <- names(celltype_counts[celltype_counts < 25])
if (length(rare_types) > 0) {
  cat("  Removing cell types with < 25 cells:", paste(rare_types, collapse = ", "), "\n")
  keep_cells <- !celltypes %in% rare_types
  celltypes <- celltypes[keep_cells]
  scrna_counts <- scrna_counts[, names(celltypes)]
  cat("  After filtering:", ncol(scrna_counts), "cells,",
      length(unique(celltypes)), "cell types\n")
} else {
  cat("  No rare cell types to remove\n")
}

# Drop unused factor levels
cell_names <- names(celltypes)
celltypes <- as.character(celltypes)
names(celltypes) <- cell_names

# Load spatial counts
spatial_counts_file <- file.path(dataset_dir, "spatial_counts.csv")
if (!file.exists(spatial_counts_file)) {
  stop("Spatial counts file not found: ", spatial_counts_file)
}
spatial_counts <- fread(spatial_counts_file)
spatial_gene_names <- spatial_counts[[1]]
spatial_counts <- as.matrix(spatial_counts[, -1, with = FALSE])
rownames(spatial_counts) <- spatial_gene_names
cat("  Spatial counts loaded:", nrow(spatial_counts), "genes x", ncol(spatial_counts), "spots\n")

# ==============================================================================
# Apply normalization
# ==============================================================================

cat("\nApplying normalization...\n")

# Prepare column data for reference
ref_col_data <- DataFrame(
  celltype = celltypes[colnames(scrna_counts)],
  row.names = colnames(scrna_counts)
)

# Prepare column data for query (empty for spatial)
query_col_data <- DataFrame(row.names = colnames(spatial_counts))

# Check for cached normalized counts
ref_cache_file <- file.path(dataset_dir, paste0("scRNA_counts_", ref_norm, ".csv"))
query_cache_file <- file.path(dataset_dir, paste0("spatial_counts_", query_norm, ".csv"))

# Normalize reference (or load from cache)
if (file.exists(ref_cache_file)) {
  cat(sprintf("  Loading cached reference normalization from %s\n", ref_cache_file))
  ref_cached <- fread(ref_cache_file, header = TRUE)
  ref_gene_names <- ref_cached[[1]]
  ref_mat <- as.matrix(ref_cached[, -1, with = FALSE])
  rownames(ref_mat) <- ref_gene_names
  # Subset to filtered cells (in case cache has pre-filtering data)
  keep_cols <- intersect(colnames(ref_mat), names(celltypes))
  ref_mat <- ref_mat[, keep_cols]
  reference <- SingleCellExperiment(
    assays = setNames(list(ref_mat), ref_norm),
    colData = ref_col_data[colnames(ref_mat), , drop = FALSE]
  )
  cat(sprintf("    Loaded. Features: %d, Samples: %d\n", nrow(reference), ncol(reference)))
} else {
  reference <- apply_normalization(scrna_counts, ref_col_data, ref_norm, is_reference = TRUE)
  # Save normalized counts to cache
  cat(sprintf("  Saving reference normalization to %s\n", ref_cache_file))
  ref_mat <- as.matrix(assay(reference, ref_norm))
  ref_out <- data.table(gene = rownames(ref_mat), as.data.table(ref_mat))
  fwrite(ref_out, ref_cache_file)
}

# Normalize query (or load from cache)
if (file.exists(query_cache_file)) {
  cat(sprintf("  Loading cached query normalization from %s\n", query_cache_file))
  query_cached <- fread(query_cache_file)
  query_gene_names <- query_cached[[1]]
  query_mat <- as.matrix(query_cached[, -1, with = FALSE])
  rownames(query_mat) <- query_gene_names
  query <- SingleCellExperiment(
    assays = setNames(list(query_mat), query_norm),
    colData = query_col_data[colnames(query_mat), , drop = FALSE]
  )
  cat(sprintf("    Loaded. Features: %d, Samples: %d\n", nrow(query), ncol(query)))
} else {
  query <- apply_normalization(spatial_counts, query_col_data, query_norm, is_reference = FALSE)
  # Save normalized counts to cache
  cat(sprintf("  Saving query normalization to %s\n", query_cache_file))
  query_mat <- as.matrix(assay(query, query_norm))
  query_out <- data.table(gene = rownames(query_mat), as.data.table(query_mat))
  fwrite(query_out, query_cache_file)
}

# ==============================================================================
# Find common genes and subset
# ==============================================================================

cat("\nFinding common genes...\n")

common_genes <- intersect(rownames(reference), rownames(query))
cat("  Common genes:", length(common_genes), "\n")

if (length(common_genes) < 50) {
  warning("Very few common genes found. Results may be unreliable.")
}

# HVG selection: pick top n_hvg by variance in the reference
if (n_hvg > 0 && length(common_genes) > n_hvg) {
  cat(sprintf("  Selecting top %d HVGs from reference...\n", n_hvg))
  ref_mat_tmp <- as.matrix(assay(reference, ref_norm)[common_genes, ])
  gene_vars <- apply(ref_mat_tmp, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_hvg)]
  common_genes <- top_genes
  cat(sprintf("  After HVG selection: %d genes\n", length(common_genes)))
  rm(ref_mat_tmp, gene_vars, top_genes)
}

# Subset to common genes
reference <- reference[common_genes, ]
query <- query[common_genes, ]

# ==============================================================================
# Run PhiSpace
# ==============================================================================

cat("\nRunning PhiSpace...\n")

# Start timing
start_time <- Sys.time()

# Run PhiSpace with the normalized assays
PhiRes <- tryCatch({
  PhiSpaceR_1ref(
    reference = reference,
    query = query,
    phenotypes = "celltype",
    refAssay = ref_norm,
    queryAssay = query_norm,
    regMethod = "PLS",
    center = TRUE,
    scale = FALSE
  )
}, error = function(e) {
  cat("Error in PhiSpace:", conditionMessage(e), "\n")
  return(NULL)
})

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (is.null(PhiRes)) {
  cat("PhiSpace failed. Exiting.\n")
  quit(status = 1)
}

cat("  PhiSpace completed in", round(runtime, 2), "seconds\n")

# ==============================================================================
# Extract and process results using Score2Prob
# ==============================================================================

cat("\nExtracting results...\n")

raw_scores <- PhiRes$PhiSpaceScore
cat("  Raw PhiSpace scores shape:", nrow(raw_scores), "x", ncol(raw_scores), "\n")

# --- Normalised scores ---
phi_scores <- normPhiScores(raw_scores)
phi_proportions <- Score2Prob(phi_scores)
result_df <- as.data.frame(phi_proportions)
rownames(result_df) <- colnames(query)
cat("  Normalised proportions computed:", nrow(result_df), "x", ncol(result_df), "\n")

# --- Unnormalised scores ---
phi_proportions_unnorm <- Score2Prob(raw_scores)
result_df_unnorm <- as.data.frame(phi_proportions_unnorm)
rownames(result_df_unnorm) <- colnames(query)
cat("  Unnormalised proportions computed:", nrow(result_df_unnorm), "x", ncol(result_df_unnorm), "\n")

# --- Min-max rescaled normalised scores (linear, no Score2Prob) ---
phi_minmax <- t(apply(phi_scores, 1, function(x) {
  x <- x - min(x)
  if (max(x) > 0) x <- x / max(x)
  x / sum(x)
}))
result_df_minmax <- as.data.frame(phi_minmax)
rownames(result_df_minmax) <- colnames(query)
cat("  Min-max proportions computed:", nrow(result_df_minmax), "x", ncol(result_df_minmax), "\n")

# --- Raw normalised scores (no Score2Prob, no rescaling) ---
result_df_rawscore <- as.data.frame(phi_scores)
rownames(result_df_rawscore) <- colnames(query)
cat("  Raw normalised scores:", nrow(result_df_rawscore), "x", ncol(result_df_rawscore),
    "range [", round(min(phi_scores), 3), ",", round(max(phi_scores), 3), "]\n")

# ==============================================================================
# Save results
# ==============================================================================

cat("\nSaving results...\n")

# Save normalised result
output_file <- file.path(output_dir, "PhiSpace_result.txt")
write.csv(result_df, output_file, row.names = TRUE)
cat("  Normalised results saved to:", output_file, "\n")

# Save unnormalised result
output_file_unnorm <- file.path(output_dir, "PhiSpace_unnorm_result.txt")
write.csv(result_df_unnorm, output_file_unnorm, row.names = TRUE)
cat("  Unnormalised results saved to:", output_file_unnorm, "\n")

# Save min-max rescaled result
output_file_minmax <- file.path(output_dir, "PhiSpace_minmax_result.txt")
write.csv(result_df_minmax, output_file_minmax, row.names = TRUE)
cat("  Min-max results saved to:", output_file_minmax, "\n")

# Save raw normalised scores (before Score2Prob)
output_file_rawscore <- file.path(output_dir, "PhiSpace_rawscore_result.txt")
write.csv(result_df_rawscore, output_file_rawscore, row.names = TRUE)
cat("  Raw normalised scores saved to:", output_file_rawscore, "\n")

# Save runtime and configuration info
runtime_file <- file.path(output_dir, "PhiSpace_runtime.txt")
write.csv(
  data.frame(
    method = "PhiSpace",
    ref_normalization = ref_norm,
    query_normalization = query_norm,
    runtime_seconds = runtime,
    n_spots = nrow(result_df),
    n_celltypes = ncol(result_df),
    n_common_genes = length(common_genes)
  ),
  runtime_file,
  row.names = FALSE
)
cat("  Runtime info saved to:", runtime_file, "\n")

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("PhiSpace benchmarking completed successfully!\n")
cat(sprintf("  Reference normalization: %s\n", ref_norm))
cat(sprintf("  Query normalization: %s\n", query_norm))
cat(sprintf("  Runtime: %.2f seconds\n", runtime))
cat(paste(rep("=", 60), collapse=""), "\n")
