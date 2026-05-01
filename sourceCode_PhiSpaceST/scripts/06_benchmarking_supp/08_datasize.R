## 08_datasize.R
## Collect and store dataset size information (n_cells, n_spots, n_common_genes, etc.)
## Skips computation if figures/dataset_sizes.csv already exists.

library(data.table)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  base_dir <- dirname(normalizePath(script_path))
} else {
  base_dir <- getwd()
}

outdir <- file.path(base_dir, "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(outdir, "dataset_sizes.csv")

# Skip if results already exist
if (file.exists(out_file)) {
  cat("Dataset sizes already computed:", out_file, "\n")
  dataset_info <- fread(out_file)
  print(dataset_info)
  cat("\nTo recompute, delete", out_file, "and re-run this script.\n")
  q(save = "no")
}

datasets <- paste0("dataset", 1:32)

collect_dataset_info <- function() {
  rows <- list()

  for (ds in datasets) {
    ds_num <- as.integer(sub("dataset", "", ds))
    ds_dir <- file.path(base_dir, "prepared_data", ds)

    # Read metadata for n_cells and n_spots
    meta_file <- file.path(ds_dir, "metadata.csv")
    if (!file.exists(meta_file)) next
    meta <- read.csv(meta_file)

    # Compute common genes from row names of count CSVs
    scrna_file <- file.path(ds_dir, "scRNA_counts.csv")
    spatial_file <- file.path(ds_dir, "spatial_counts.csv")

    if (file.exists(scrna_file) && file.exists(spatial_file)) {
      scrna_genes <- fread(scrna_file, select = 1, header = TRUE)[[1]]
      spatial_genes <- fread(spatial_file, select = 1, header = TRUE)[[1]]
      n_common <- length(intersect(scrna_genes, spatial_genes))
    } else {
      n_common <- NA
    }

    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds_num,
      n_cells = meta$n_cells,
      n_spots = meta$n_spots,
      n_genes_scrna = meta$n_genes_scrna,
      n_genes_spatial = meta$n_genes_spatial,
      n_celltypes = meta$n_celltypes,
      n_common_genes = n_common
    )
  }

  rbindlist(rows)
}

dataset_info <- collect_dataset_info()

cat("=== Dataset sizes ===\n")
print(dataset_info)
cat("\nSummary:\n")
cat(sprintf("  n_cells:        %d (constant)\n", dataset_info$n_cells[1]))
cat(sprintf("  n_spots:        %d (constant)\n", dataset_info$n_spots[1]))
cat(sprintf("  n_common_genes: %d - %d (median %d)\n",
            min(dataset_info$n_common_genes, na.rm = TRUE),
            max(dataset_info$n_common_genes, na.rm = TRUE),
            median(dataset_info$n_common_genes, na.rm = TRUE)))
cat(sprintf("  n_celltypes:    %d - %d\n",
            min(dataset_info$n_celltypes), max(dataset_info$n_celltypes)))

# Save
write.csv(dataset_info, out_file, row.names = FALSE)
cat("\nDataset sizes saved to", out_file, "\n")
