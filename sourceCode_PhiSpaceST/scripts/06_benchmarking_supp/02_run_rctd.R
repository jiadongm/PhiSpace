#!/usr/bin/env Rscript

# ==============================================================================
# RCTD Deconvolution Pipeline for Benchmarking
# ==============================================================================
# This script runs RCTD (spacexr) on prepared benchmark datasets and outputs
# cell type proportions in the standard format used by the benchmarking framework.
#
# Usage:
#   Rscript 02_run_rctd.R <dataset_dir> <output_dir>
#
# Arguments:
#   dataset_dir: Directory containing prepared dataset files
#                (scRNA_counts.csv, scRNA_celltype.csv, spatial_counts.csv)
#   output_dir: Directory for output results
# ==============================================================================

suppressPackageStartupMessages({
  library(spacexr)
  library(Matrix)
  library(data.table)
})

# ==============================================================================
# Parse command line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 02_run_rctd.R <dataset_dir> <output_dir>\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  dataset_dir: Directory containing prepared dataset files\n")
  cat("  output_dir: Directory for output results\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript 02_run_rctd.R prepared_data/dataset1 results_rctd/dataset1\n")
  quit(status = 1)
}

dataset_dir <- args[1]
output_dir <- args[2]

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("RCTD Benchmarking Pipeline\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

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
# Filter cell types with < 25 cells
# ==============================================================================

cat("\nFiltering rare cell types...\n")

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

# Drop unused factor levels and preserve names
cell_names <- names(celltypes)
celltypes <- as.factor(as.character(celltypes))
names(celltypes) <- cell_names

# ==============================================================================
# Create RCTD reference and spatial objects
# ==============================================================================

cat("\nPreparing RCTD inputs...\n")

# Convert to data.frame for spacexr (expects genes x cells data.frame)
sc_counts_df <- as.data.frame(scrna_counts)
nUMI_sc <- colSums(scrna_counts)
names(nUMI_sc) <- colnames(scrna_counts)

# Create Reference object
cat("  Creating Reference object...\n")
reference <- Reference(sc_counts_df, celltypes, nUMI_sc)

# Create SpatialRNA object with dummy coordinates
spatial_counts_df <- as.data.frame(spatial_counts)
coords <- data.frame(
  xcoord = seq_len(ncol(spatial_counts)),
  ycoord = seq_len(ncol(spatial_counts)),
  row.names = colnames(spatial_counts)
)
nUMI_spatial <- colSums(spatial_counts)

cat("  Creating SpatialRNA object...\n")
puck <- SpatialRNA(coords, spatial_counts_df, nUMI_spatial)

# ==============================================================================
# Run RCTD
# ==============================================================================

cat("\nRunning RCTD...\n")

start_time <- Sys.time()

myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("  RCTD completed in", round(runtime, 2), "seconds\n")

# ==============================================================================
# Extract and normalize results
# ==============================================================================

cat("\nExtracting results...\n")

results <- myRCTD@results
# Convert to dense matrix if needed, then normalize to sum to 1
weights_mat <- as.matrix(results$weights)
norm_weights <- sweep(weights_mat, 1, rowSums(weights_mat), "/")

cat("  Result shape:", nrow(norm_weights), "x", ncol(norm_weights), "\n")

# ==============================================================================
# Save results
# ==============================================================================

cat("\nSaving results...\n")

output_file <- file.path(output_dir, "RCTD_result.txt")
write.csv(norm_weights, output_file)
cat("  Results saved to:", output_file, "\n")

# Save runtime info
runtime_file <- file.path(output_dir, "RCTD_runtime.txt")
write.csv(
  data.frame(
    method = "RCTD",
    runtime_seconds = runtime,
    n_spots = nrow(norm_weights),
    n_celltypes = ncol(norm_weights)
  ),
  runtime_file,
  row.names = FALSE
)
cat("  Runtime info saved to:", runtime_file, "\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RCTD benchmarking completed successfully!\n")
cat(sprintf("  Runtime: %.2f seconds\n", runtime))
cat(paste(rep("=", 60), collapse = ""), "\n")
