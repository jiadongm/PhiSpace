## ============================================================================
## Prepare all 40 Visium NSCLC samples for PhiSpace analysis
## Creates per-sample directory structure, loads with Seurat, converts to SCE
## ============================================================================

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(PhiSpace)
  library(Seurat)
  library(SingleCellExperiment)
  library(qs)
  library(dplyr)
})

raw_dir <- file.path(paths$data_root, "Visium_NSCLC", "allSamples")
out_dir <- file.path(paths$output_root, "pDC", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# All 40 sample names
sampleNames <- c(
  "D1_1", "D1_2", "D2_1", "D2_2",
  "P10_B1", "P10_B2", "P10_T1", "P10_T2", "P10_T3", "P10_T4",
  "P11_B1", "P11_B2", "P11_T1", "P11_T2", "P11_T3", "P11_T4",
  "P15_B1", "P15_B2", "P15_T1", "P15_T2",
  "P16_B1", "P16_B2", "P16_T1", "P16_T2",
  "P17_B1", "P17_B2", "P17_T1", "P17_T2",
  "P19_B1", "P19_B2", "P19_T1", "P19_T2",
  "P24_B1", "P24_B2", "P24_T1", "P24_T2",
  "P25_B1", "P25_B2", "P25_T1", "P25_T2"
)

# Condition labels
cancerTypes <- c(
  "Healthy", "Healthy", "Healthy", "Healthy",  # D1, D2
  "Background", "Background", "LUAD", "LUAD", "LUAD", "LUAD",  # P10
  "Background", "Background", "LUSC", "LUSC", "LUSC", "LUSC",  # P11
  "Background", "Background", "LUAD", "LUAD",  # P15
  "Background", "Background", "LUAD", "LUAD",  # P16
  "Background", "Background", "LUSC", "LUSC",  # P17
  "Background", "Background", "LUSC", "LUSC",  # P19
  "Background", "Background", "LUAD", "LUAD",  # P24
  "Background", "Background", "LUAD", "LUAD"   # P25
)
names(cancerTypes) <- sampleNames

# Save sample metadata
write.csv(
  data.frame(sample = sampleNames, condition = cancerTypes),
  file.path(out_dir, "sample_metadata.csv"), row.names = FALSE
)

# Process each sample
query_list <- vector("list", length(sampleNames))
names(query_list) <- sampleNames

for (samp in sampleNames) {

  sce_path <- file.path(out_dir, paste0(samp, "_sce.qs"))

  if (file.exists(sce_path)) {
    cat(sprintf("  [%s] Loading existing SCE\n", samp))
    query_list[[samp]] <- qread(sce_path)
    next
  }

  cat(sprintf("  [%s] Processing from raw data...\n", samp))

  # Create per-sample directory for Seurat::Load10X_Spatial
  samp_dir <- file.path(out_dir, samp)
  dir.create(samp_dir, showWarnings = FALSE)

  # Symlink h5 file
  h5_src <- file.path(raw_dir, paste0(samp, "-filtered_feature_bc_matrix.h5"))
  h5_dst <- file.path(samp_dir, paste0(samp, "-filtered_feature_bc_matrix.h5"))
  if (!file.exists(h5_dst)) file.symlink(h5_src, h5_dst)

  # Set up spatial directory
  spatial_dst <- file.path(samp_dir, "spatial")
  if (!dir.exists(spatial_dst)) {
    # Check if already extracted in raw dir
    spatial_extracted <- file.path(raw_dir, paste0(samp, "-spatial"))
    if (dir.exists(spatial_extracted)) {
      file.symlink(spatial_extracted, spatial_dst)
    } else {
      # Extract from tar
      tar_file <- file.path(raw_dir, paste0(samp, "-spatial.tar"))
      dir.create(spatial_dst, showWarnings = FALSE)
      system2("tar", c("-xf", tar_file, "-C", spatial_dst, "--strip-components=0"))
    }
  }

  # Load with Seurat
  qu_vis <- Load10X_Spatial(
    samp_dir,
    filename = paste0(samp, "-filtered_feature_bc_matrix.h5")
  )

  # Extract spatial coordinates
  spat_info <- data.frame(
    x = qu_vis@images$slice1$centroids@coords[, 2],
    y = -qu_vis@images$slice1$centroids@coords[, 1]
  )
  qu_vis <- AddMetaData(qu_vis, spat_info)

  # Convert to SCE
  sce <- SingleCellExperiment(
    list(counts = qu_vis[["Spatial"]]$counts),
    colData = qu_vis@meta.data %>% as.data.frame()
  )
  assay(sce, "log1p") <- log1p(assay(sce, "counts"))

  # Save
  qsave(sce, sce_path)
  query_list[[samp]] <- sce
  cat(sprintf("  [%s] Done: %d genes x %d spots\n", samp, nrow(sce), ncol(sce)))
}

# Save combined list
allPath <- file.path(out_dir, "allSamples.qs")
qsave(query_list, allPath)
cat(sprintf("\nSaved combined list of %d samples to %s\n", length(query_list), allPath))
