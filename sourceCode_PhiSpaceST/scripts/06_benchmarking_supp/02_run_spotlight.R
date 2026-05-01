#!/usr/bin/env Rscript

# ==============================================================================
# SPOTlight Deconvolution Pipeline for Benchmarking
# ==============================================================================
# Runs SPOTlight (NMF-based deconvolution) on prepared benchmark datasets and
# outputs cell type proportions in the standard format.
#
# Uses SPOTlight Bioconductor API (v1.x): SPOTlight(), getMGS()
#
# Usage:
#   Rscript 02_run_spotlight.R <dataset_dir> <output_dir>
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SPOTlight)
  library(SingleCellExperiment)
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
  cat("Usage: Rscript 02_run_spotlight.R <dataset_dir> <output_dir>\n")
  quit(status = 1)
}

dataset_dir <- args[1]
output_dir <- args[2]

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SPOTlight Benchmarking Pipeline\n")
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
scrna_counts <- fread(file.path(dataset_dir, "scRNA_counts.csv"), header = TRUE)
gene_names <- scrna_counts[[1]]
scrna_counts <- as.matrix(scrna_counts[, -1, with = FALSE])
rownames(scrna_counts) <- gene_names
cat("  scRNA-seq counts loaded:", nrow(scrna_counts), "genes x",
    ncol(scrna_counts), "cells\n")

# Load cell type annotations
celltype_df <- fread(file.path(dataset_dir, "scRNA_celltype.csv"))
celltypes <- celltype_df$celltype
names(celltypes) <- as.character(celltype_df$cell)
cat("  Cell types loaded:", length(unique(celltypes)), "unique types\n")

# Load spatial counts (genes x spots)
spatial_counts <- fread(file.path(dataset_dir, "spatial_counts.csv"), header = TRUE)
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
      length(rare_types), "types\n")
  keep_cells <- !celltypes %in% rare_types
  celltypes <- celltypes[keep_cells]
  scrna_counts <- scrna_counts[, names(celltypes)]
  cat("  After filtering:", ncol(scrna_counts), "cells,",
      length(unique(celltypes)), "cell types\n")
} else {
  cat("  No rare cell types to remove\n")
}

# ==============================================================================
# Prepare objects
# ==============================================================================

cat("\nPreparing objects...\n")

# Convert to sparse
scrna_counts <- Matrix::Matrix(scrna_counts, sparse = TRUE)
spatial_counts <- Matrix::Matrix(spatial_counts, sparse = TRUE)

# Remove empty gene names
keep_sc <- nchar(rownames(scrna_counts)) > 0
if (any(!keep_sc)) {
  scrna_counts <- scrna_counts[keep_sc, ]
}
keep_sp <- nchar(rownames(spatial_counts)) > 0
if (any(!keep_sp)) {
  spatial_counts <- spatial_counts[keep_sp, ]
}

# Create SingleCellExperiment for scRNA-seq
sce <- SingleCellExperiment(
  assays = list(counts = scrna_counts),
  colData = data.frame(celltype = celltypes[colnames(scrna_counts)],
                       row.names = colnames(scrna_counts))
)
colLabels(sce) <- sce$celltype

# Spatial data passed as sparse matrix directly to SPOTlight
# (avoids Seurat v5 extraction issues in runDeconvolution)

# ==============================================================================
# Run SPOTlight
# ==============================================================================

cat("\nRunning SPOTlight...\n")

start_time <- Sys.time()

# Find marker genes using Seurat on the scRNA-seq data
sc_seurat <- CreateSeuratObject(counts = scrna_counts, project = "reference")
sc_seurat$celltype <- celltypes[colnames(sc_seurat)]
Idents(sc_seurat) <- sc_seurat$celltype

set.seed(123)
sc_seurat <- SCTransform(sc_seurat, verbose = FALSE)

cluster_markers <- FindAllMarkers(
  object = sc_seurat,
  assay = "SCT",
  slot = "data",
  verbose = FALSE,
  only.pos = TRUE
)

# Format marker genes for SPOTlight's mgs parameter
mgs <- cluster_markers[, c("gene", "cluster", "avg_log2FC")]
colnames(mgs) <- c("gene", "cluster", "weight")
mgs$cluster <- as.character(mgs$cluster)

# Run SPOTlight deconvolution (Bioconductor API)
# Pass spatial as matrix (Seurat v5 object causes issues with runDeconvolution)
res <- SPOTlight(
  x = sce,
  y = spatial_counts,
  groups = as.character(sce$celltype),
  mgs = mgs,
  gene_id = "gene",
  group_id = "cluster",
  weight_id = "weight",
  n_top = NULL,
  scale = TRUE,
  min_prop = 0,
  verbose = TRUE
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("  SPOTlight completed in", round(runtime, 2), "seconds\n")

# ==============================================================================
# Extract and save results
# ==============================================================================

cat("\nExtracting results...\n")

# SPOTlight returns a list with 'mat' (proportions) and 'NMF' (model)
decon_mtrx <- res$mat
# Rows are spots, columns are cell types

cat("  Result shape:", nrow(decon_mtrx), "x", ncol(decon_mtrx), "\n")

output_file <- file.path(output_dir, "SPOTlight_result.txt")
write.csv(decon_mtrx, output_file)
cat("  Results saved to:", output_file, "\n")

# Save runtime info
runtime_file <- file.path(output_dir, "SPOTlight_runtime.txt")
write.csv(
  data.frame(
    method = "SPOTlight",
    runtime_seconds = runtime,
    n_spots = nrow(decon_mtrx),
    n_celltypes = ncol(decon_mtrx)
  ),
  runtime_file,
  row.names = FALSE
)
cat("  Runtime info saved to:", runtime_file, "\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SPOTlight benchmarking completed successfully!\n")
cat(sprintf("  Runtime: %.2f seconds\n", runtime))
cat(paste(rep("=", 60), collapse = ""), "\n")
