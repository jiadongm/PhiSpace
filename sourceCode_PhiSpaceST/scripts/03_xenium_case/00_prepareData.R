# prepareData.R
# Prepare Xenium FFPE Human Lung Cancer data for PhiSpace analysis
# Download source: https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(qs)
  library(Matrix)
  library(DropletUtils) # for read10xCounts (reads h5 format)
})

# ---- Paths ----
xenium_dir <- "data/Xenium_NSCLC"
dat_dir <- paste0(paths$data_root, "/")

# ---- 1. Load Xenium cell-feature matrix ----
cat("Loading Xenium cell-feature matrix...\n")
h5_path <- file.path(xenium_dir, "Xenium_V1_humanLung_Cancer_FFPE_cell_feature_matrix.h5")
sce_raw <- read10xCounts(h5_path, col.names = TRUE)
cat(sprintf("  Raw: %d genes x %d cells\n", nrow(sce_raw), ncol(sce_raw)))

# ---- 2. Load cell metadata (spatial coordinates) ----
cat("Loading cell metadata...\n")
cells_path <- file.path(xenium_dir, "Xenium_V1_humanLung_Cancer_FFPE_cells.csv.gz")
cell_meta <- read.csv(cells_path)
cat(sprintf("  Cell metadata: %d cells x %d columns\n", nrow(cell_meta), ncol(cell_meta)))
cat(sprintf("  Columns: %s\n", paste(colnames(cell_meta), collapse = ", ")))

# ---- 3. Match cells between h5 and metadata ----
# The h5 barcodes should match cell_id in metadata
h5_barcodes <- colnames(sce_raw)
meta_ids <- cell_meta$cell_id
cat(sprintf("  h5 barcodes: %d, metadata cell_ids: %d\n", length(h5_barcodes), length(meta_ids)))

# Match ordering
common_cells <- intersect(h5_barcodes, meta_ids)
cat(sprintf("  Common cells: %d\n", length(common_cells)))

# Subset and reorder
sce_raw <- sce_raw[, common_cells]
cell_meta <- cell_meta[match(common_cells, cell_meta$cell_id), ]
rownames(cell_meta) <- cell_meta$cell_id

# ---- 4. Add spatial coords and metadata to colData ----
colData(sce_raw) <- DataFrame(cell_meta)

# ---- 5. Filter genes: keep only gene targets (exclude controls) ----
gene_info <- rowData(sce_raw)
cat(sprintf("\n  Gene types in the data:\n"))
print(table(gene_info$Type))

# Keep only "Gene Expression" features
gene_idx <- gene_info$Type == "Gene Expression"
sce <- sce_raw[gene_idx, ]
cat(sprintf("\n  After filtering controls: %d genes x %d cells\n", nrow(sce), ncol(sce)))

# Use gene symbols as rownames
rownames(sce) <- rowData(sce)$Symbol

# ---- 6. Basic QC and filtering ----
cat("\nQC filtering...\n")
# Per-cell QC
qc <- perCellQCMetrics(sce)
colData(sce)$detected <- qc$detected
colData(sce)$sum <- qc$sum

cat(sprintf("  Median genes detected per cell: %.0f\n", median(qc$detected)))
cat(sprintf("  Median total counts per cell: %.0f\n", median(qc$sum)))

# Filter cells with 0 transcript counts
keep_cells <- qc$sum > 0
cat(sprintf("  Cells with >0 counts: %d / %d\n", sum(keep_cells), length(keep_cells)))
sce <- sce[, keep_cells]

# Filter genes detected in at least 10 cells
gene_detected <- rowSums(assay(sce, "counts") > 0)
keep_genes <- gene_detected >= 10
cat(sprintf("  Genes detected in >=10 cells: %d / %d\n", sum(keep_genes), length(keep_genes)))
sce <- sce[keep_genes, ]

cat(sprintf("\n  Final SCE: %d genes x %d cells\n", nrow(sce), ncol(sce)))

# ---- 7. Normalisation: log1p ----
cat("Computing log1p normalisation...\n")
assay(sce, "log1p") <- log1p(assay(sce, "counts"))

# ---- 8. Save SCE object ----
out_path <- file.path(xenium_dir, "xenium_lung_sce.qs")
cat(sprintf("Saving SCE to %s...\n", out_path))
qsave(sce, out_path)

# ---- 9. Check gene overlap with 4 references ----
cat("\n=== Gene overlap with scRNA-seq references ===\n")
xenium_genes <- rownames(sce)
lineages <- c("immune", "epithelial", "endothelial", "mesenchymal")
for (lineage in lineages) {
  ref_path <- paste0(dat_dir, "LungFibrosis/scRNA-seq/", lineage, "_sce.qs")
  if (file.exists(ref_path)) {
    ref <- qread(ref_path)
    ref_genes <- rownames(ref)
    overlap <- intersect(xenium_genes, ref_genes)
    cat(sprintf("  %s: %d ref genes, %d overlap with Xenium (%d Xenium genes)\n",
                lineage, length(ref_genes), length(overlap), length(xenium_genes)))
  } else {
    cat(sprintf("  %s: reference file not found at %s\n", lineage, ref_path))
  }
}

# Also check Azimuth lung reference
azimuth_path <- paste0(dat_dir, "LungRef/AzimuthLung2.0_sce_0.1sub.qs")
if (file.exists(azimuth_path)) {
  ref <- qread(azimuth_path)
  ref_genes <- rownames(ref)
  overlap <- intersect(xenium_genes, ref_genes)
  cat(sprintf("  Azimuth Lung: %d ref genes, %d overlap with Xenium\n",
              length(ref_genes), length(overlap)))
}


# ---- 10. Load pathology annotations from Xenium Explorer coordinates ----
# Polygon boundaries exported from Xenium Explorer (in Xenium micron space).
# Each row is a vertex; polygons grouped by Selection column; Class gives the
# annotation type (Tumor, Blood Vessels, Lymphoid, Unclassified).
# "Unclassified" covers Normal lung and Benign bronchial epithelium selections.
cat("\n=== Pathology annotations from Xenium Explorer ===\n")
library(readr)
library(sf)

annot_path <- file.path(xenium_dir, "patho_coordinates.csv")
patho_coordinates <- read_csv(annot_path, skip = 2, show_col_types = FALSE)
cat(sprintf("  Loaded %d polygon vertices across %d selections\n",
            nrow(patho_coordinates), length(unique(patho_coordinates$Selection))))
cat("  Classes:", paste(unique(patho_coordinates$Class), collapse = ", "), "\n")

# Build sf polygons per Selection
selections <- split(patho_coordinates, patho_coordinates$Selection)
poly_list <- lapply(names(selections), function(sel) {
  df <- selections[[sel]]
  coords <- cbind(df$X, df$Y)
  # Close the ring if not already closed
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, , drop = FALSE])
  }
  st_polygon(list(coords))
})

poly_sf <- st_sf(
  Selection = names(selections),
  Class = sapply(selections, function(df) df$Class[1]),
  geometry = st_sfc(poly_list)
)
cat(sprintf("  Built %d polygons\n", nrow(poly_sf)))

# Create sf points from cell coordinates
cell_df <- data.frame(
  cell_id = colnames(sce),
  x = colData(sce)$x_centroid,
  y = colData(sce)$y_centroid
)
cell_pts <- st_as_sf(cell_df, coords = c("x", "y"))

# Point-in-polygon: assign each cell to the smallest containing polygon
# (smallest = most specific, e.g. Lymphoid inside Tumor)
poly_sf$area <- st_area(poly_sf)
poly_sf <- poly_sf[order(poly_sf$area), ]  # smallest first

cat("  Running point-in-polygon...\n")
# st_join with largest = FALSE keeps only the first match per point
# Since polygons are sorted smallest-first, this gives the most specific label
joined <- st_join(cell_pts, poly_sf, join = st_within, largest = FALSE)

# Handle cells in multiple polygons: keep smallest (first after sort)
joined <- joined[!duplicated(joined$cell_id), ]

# Match back to sce column order
m <- match(colnames(sce), joined$cell_id)
sce$pathology <- joined$Class[m]
sce$pathology[is.na(sce$pathology)] <- "Unannotated"
sce$pathology_detail <- joined$Selection[m]
sce$pathology_detail[is.na(sce$pathology_detail)] <- "Unannotated"

cat("  Pathology distribution:\n")
print(table(sce$pathology))
cat("\n  Detailed:\n")
print(table(sce$pathology_detail))

# Re-save SCE with pathology annotations
qsave(sce, out_path)
cat("  SCE updated with pathology annotations.\n")



cat("\nDone!\n")















