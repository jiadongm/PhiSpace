#!/usr/bin/env Rscript

# ==============================================================================
# Seurat Label Transfer Deconvolution Pipeline for Benchmarking
# ==============================================================================
# Runs Seurat's FindTransferAnchors + TransferData on prepared benchmark
# datasets and outputs cell type proportions in the standard format.
#
# Usage:
#   Rscript 02_run_seurat.R <dataset_dir> <output_dir>
#
# Arguments:
#   dataset_dir: Directory containing prepared dataset files
#                (scRNA_counts.csv, scRNA_celltype.csv, spatial_counts.csv)
#   output_dir: Directory for output results
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# Allow large globals for SCTransform parallelisation
options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB

# ==============================================================================
# Parse command line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 02_run_seurat.R <dataset_dir> <output_dir>\n")
  cat("\nArguments:\n")
  cat("  dataset_dir: Directory containing prepared dataset files\n")
  cat("  output_dir: Directory for output results\n")
  cat("\nExample:\n")
  cat("  Rscript 02_run_seurat.R prepared_data/dataset1 results_seurat/dataset1\n")
  quit(status = 1)
}

dataset_dir <- args[1]
output_dir <- args[2]

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Seurat Benchmarking Pipeline\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load data
# ==============================================================================

cat("\nLoading data...\n")

# Load scRNA-seq counts (genes x cells)
scrna_counts_file <- file.path(dataset_dir, "scRNA_counts.csv")
if (!file.exists(scrna_counts_file)) {
  stop("scRNA counts file not found: ", scrna_counts_file)
}
scrna_counts <- fread(scrna_counts_file, header = TRUE)
gene_names <- scrna_counts[[1]]
scrna_counts <- as.matrix(scrna_counts[, -1, with = FALSE])
rownames(scrna_counts) <- gene_names
cat("  scRNA-seq counts loaded:", nrow(scrna_counts), "genes x",
    ncol(scrna_counts), "cells\n")

# Load cell type annotations
celltype_file <- file.path(dataset_dir, "scRNA_celltype.csv")
if (!file.exists(celltype_file)) {
  stop("Cell type file not found: ", celltype_file)
}
celltype_df <- fread(celltype_file)
celltypes <- celltype_df$celltype
names(celltypes) <- as.character(celltype_df$cell)
cat("  Cell types loaded:", length(unique(celltypes)), "unique types\n")

# Load spatial counts (genes x spots)
spatial_counts_file <- file.path(dataset_dir, "spatial_counts.csv")
if (!file.exists(spatial_counts_file)) {
  stop("Spatial counts file not found: ", spatial_counts_file)
}
spatial_counts <- fread(spatial_counts_file)
spatial_gene_names <- spatial_counts[[1]]
spatial_counts <- as.matrix(spatial_counts[, -1, with = FALSE])
rownames(spatial_counts) <- spatial_gene_names
cat("  Spatial counts loaded:", nrow(spatial_counts), "genes x",
    ncol(spatial_counts), "spots\n")

# ==============================================================================
# Filter cell types with < 25 cells (matching RCTD/PhiSpace)
# ==============================================================================

cat("\nFiltering rare cell types...\n")

celltype_counts <- table(celltypes)
rare_types <- names(celltype_counts[celltype_counts < 25])

if (length(rare_types) > 0) {
  cat("  Removing cell types with < 25 cells:",
      paste(rare_types, collapse = ", "), "\n")
  keep_cells <- !celltypes %in% rare_types
  celltypes <- celltypes[keep_cells]
  scrna_counts <- scrna_counts[, names(celltypes)]
  cat("  After filtering:", ncol(scrna_counts), "cells,",
      length(unique(celltypes)), "cell types\n")
} else {
  cat("  No rare cell types to remove\n")
}

# ==============================================================================
# Create Seurat objects
# ==============================================================================

cat("\nCreating Seurat objects...\n")

# Convert to sparse matrices (avoids Seurat v5 LogMap issue with dense matrices)
scrna_counts <- Matrix::Matrix(scrna_counts, sparse = TRUE)
spatial_counts <- Matrix::Matrix(spatial_counts, sparse = TRUE)

# Remove rows with empty gene names (can arise from CSV index columns)
keep_sc <- nchar(rownames(scrna_counts)) > 0
if (any(!keep_sc)) {
  cat("  Removing", sum(!keep_sc), "empty gene names from scRNA\n")
  scrna_counts <- scrna_counts[keep_sc, ]
}
keep_sp <- nchar(rownames(spatial_counts)) > 0
if (any(!keep_sp)) {
  cat("  Removing", sum(!keep_sp), "empty gene names from spatial\n")
  spatial_counts <- spatial_counts[keep_sp, ]
}

# Reference (scRNA-seq)
sc_rna <- CreateSeuratObject(counts = scrna_counts, project = "reference")
sc_rna$celltype <- celltypes[colnames(sc_rna)]

# Query (spatial)
spatial <- CreateSeuratObject(counts = spatial_counts, project = "spatial")

# ==============================================================================
# Run Seurat label transfer
# ==============================================================================

cat("\nRunning Seurat SCTransform + label transfer...\n")

start_time <- Sys.time()

# SCTransform both
sc_rna <- SCTransform(sc_rna, verbose = FALSE)
spatial <- SCTransform(spatial, verbose = FALSE)

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = sc_rna,
  query = spatial,
  dims = 1:30,
  normalization.method = "SCT"
)

# Transfer cell type labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = sc_rna$celltype,
  dims = 1:30
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("  Seurat completed in", round(runtime, 2), "seconds\n")

# ==============================================================================
# Extract and save results
# ==============================================================================

cat("\nExtracting results...\n")

# predictions contains: predicted.id, prediction.score.{celltype}, prediction.score.max
# We need the prediction.score.{celltype} columns as proportions
score_cols <- grep("^prediction\\.score\\.", colnames(predictions), value = TRUE)
score_cols <- setdiff(score_cols, "prediction.score.max")

result <- predictions[, score_cols]
# Column names keep the prediction.score. prefix — 03_calculate_metrics.py
# strips this prefix when method_name == "Seurat"

cat("  Result shape:", nrow(result), "x", ncol(result), "\n")

output_file <- file.path(output_dir, "Seurat_result.txt")
write.csv(result, output_file)
cat("  Results saved to:", output_file, "\n")

# Save runtime info
runtime_file <- file.path(output_dir, "Seurat_runtime.txt")
write.csv(
  data.frame(
    method = "Seurat",
    runtime_seconds = runtime,
    n_spots = nrow(result),
    n_celltypes = ncol(result)
  ),
  runtime_file,
  row.names = FALSE
)
cat("  Runtime info saved to:", runtime_file, "\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Seurat benchmarking completed successfully!\n")
cat(sprintf("  Runtime: %.2f seconds\n", runtime))
cat(paste(rep("=", 60), collapse = ""), "\n")
