# runPhiSpace4Refs.R
# Run PhiSpace annotation on Xenium lung cancer data using 4 scRNA-seq references
# Analogous to CaseCosMx/runPhiSpace4Refs.R

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(PhiSpace)
  library(SingleCellExperiment)
  library(scuttle)
  library(scran)
  library(qs)
  library(ggplot2)
})

# ---- Paths ----
ref_dir <- file.path(paths$data_root, "LungFibrosis", "scRNA-seq/")
xenium_dir <- "data/Xenium_NSCLC"
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Parameters ----
YtrainName <- "manual_annotation_1"
refAssay <- "logcounts"  # scran-normalised (already in references)
queryAssay <- "log1p"    # simple log1p for Xenium query
minCellsPerType <- 25    # filter rare cell types

# ---- 1. Load Xenium query ----
cat("Loading Xenium query SCE...\n")
query <- qread(file.path(xenium_dir, "xenium_lung_sce.qs"))
cat(sprintf("  Query: %d genes x %d cells\n", nrow(query), ncol(query)))

# Ensure log1p assay exists
if (!queryAssay %in% assayNames(query)) {
  cat(sprintf("  Computing %s assay for query...\n", queryAssay))
  assay(query, queryAssay) <- log1p(assay(query, "counts"))
}

# ---- 2. Load references ----
cat("\nLoading 4 lineage-specific references...\n")
lineages <- c("immune", "epithelial", "endothelial", "mesenchymal")
ref_list <- lapply(lineages, function(lineage) {
  ref <- qread(paste0(ref_dir, lineage, "_sce.qs"))
  cat(sprintf("  %s: %d genes x %d cells\n", lineage, nrow(ref), ncol(ref)))
  ref
})
names(ref_list) <- lineages

# ---- 3. Process references: scran normalisation + cell type filtering ----
cat("\nProcessing references...\n")
for (ii in seq_along(ref_list)) {
  lineage <- lineages[ii]
  ref <- ref_list[[lineage]]

  # Check if scran normalisation (logcounts) already exists
  if (refAssay %in% assayNames(ref)) {
    cat(sprintf("  %s: '%s' assay already present, skipping normalisation\n",
                lineage, refAssay))
  } else {
    cat(sprintf("  %s: computing scran normalisation...\n", lineage))
    # Compute scran size factors
    clusters <- quickCluster(ref)
    ref <- computeSumFactors(ref, clusters = clusters)
    ref <- logNormCounts(ref)  # creates 'logcounts' assay
    cat(sprintf("    Done. Size factor range: [%.3f, %.3f]\n",
                min(sizeFactors(ref)), max(sizeFactors(ref))))
  }

  # Filter cell types with fewer than minCellsPerType cells
  ct_counts <- table(ref[[YtrainName]])
  rare_types <- names(ct_counts[ct_counts < minCellsPerType])
  if (length(rare_types) > 0) {
    cat(sprintf("  %s: removing %d rare cell types (<%d cells): %s\n",
                lineage, length(rare_types), minCellsPerType,
                paste(rare_types, collapse = ", ")))
    keep <- !ref[[YtrainName]] %in% rare_types
    ref <- ref[, keep]
    # Drop unused factor levels
    ref[[YtrainName]] <- droplevels(as.factor(ref[[YtrainName]]))
    cat(sprintf("    After filtering: %d cells, %d cell types\n",
                ncol(ref), length(unique(ref[[YtrainName]]))))
  } else {
    cat(sprintf("  %s: all cell types have >= %d cells\n", lineage, minCellsPerType))
  }

  ref_list[[lineage]] <- ref
}

# ---- 4. Run PhiSpace for each reference ----
cat("\n=== Running PhiSpace ===\n")
sc_list <- vector("list", length(ref_list))
names(sc_list) <- lineages

for (ii in seq_along(ref_list)) {
  lineage <- lineages[ii]
  reference <- ref_list[[lineage]]

  cat(sprintf("\n--- %s (%d cell types) ---\n", lineage,
              length(unique(reference[[YtrainName]]))))

  # Find overlapping genes
  common_genes <- intersect(rownames(reference), rownames(query))
  cat(sprintf("  Overlapping genes: %d\n", length(common_genes)))

  t0 <- proc.time()

  PhiRes <- PhiSpaceR_1ref(
    reference = reference[common_genes, ],
    query = query[common_genes, ],
    phenotypes = YtrainName,
    refAssay = refAssay,
    queryAssay = queryAssay,
    regMethod = "PLS",
    center = TRUE,
    scale = FALSE
  )

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  PhiSpace completed in %.1f seconds\n", elapsed))

  # Normalise scores and add lineage suffix to column names
  sc <- PhiRes$PhiSpaceScore
  sc_norm <- normPhiScores(sc)
  colnames(sc_norm) <- paste0(colnames(sc_norm), "(", lineage, ")")

  cat(sprintf("  Output: %d cells x %d cell types\n", nrow(sc_norm), ncol(sc_norm)))
  sc_list[[lineage]] <- sc_norm
}

# ---- 5. Combine scores from all references ----
cat("\nCombining PhiSpace scores from all references...\n")
combined_scores <- Reduce(cbind, sc_list)
cat(sprintf("  Combined: %d cells x %d cell types\n",
            nrow(combined_scores), ncol(combined_scores)))

# Store in query SCE
reducedDim(query, "PhiSpace") <- combined_scores

# ---- 6. Save results ----
# Save individual scores per reference
phi_res_path <- file.path(output_dir, "PhiRes4Refs.qs")
cat(sprintf("\nSaving per-reference scores to %s\n", phi_res_path))
qsave(sc_list, phi_res_path)


# Define pathological regions
pathology_combined <- query$pathology_detail
pathology_combined[query$pathology_detail %in% c("Blood vessels", "Blood vessels (1)", "Blood vessels (2)")] <- "BloodVessels"
pathology_combined[query$pathology_detail %in% c("Lymphoid aggregate", "Lymphoid aggregate (1)", "Lymphoid aggregate (2)", "Lymphoid aggregate (3)")] <- "Lymphoid"
pathology_combined[query$pathology_detail %in% c("Tumor", "Tumor (1)")] <- "Tumour"
pathology_combined[query$pathology_detail == "Benign bronchial epithelium"] <- "Normal"
pathology_combined[query$pathology_detail == "Normal lung"] <- "Normal"
query$pathology_combined <- pathology_combined

# Save query SCE with PhiSpace scores
sce_out_path <- file.path(xenium_dir, "xenium_lung_sce.qs")
cat(sprintf("Saving annotated SCE to %s\n", sce_out_path))
qsave(query, sce_out_path)

# ---- 7. Save spatial cell type heatmaps ----
cat("\nSaving spatial cell type heatmaps...\n")
# Flip y to match Xenium Explorer / image convention (y increases downward)
query$y_centroid_flipped <- -query$y_centroid
fig_dir <- "figs/spatialCellTypeHeatmaps"
saveCellTypeMaps(
  query,
  methodName = "PhiSpace",
  tissueName = "XeniumLungCancer",
  coordNames = c("x_centroid", "y_centroid_flipped"),
  freeColScale = F,
  psize = 0.4,
  outputDir = fig_dir,
  fignrow = 2,
  figncol = 2,
  width = 12,
  height = 12
)
cat(sprintf("  Heatmaps saved to %s\n", fig_dir))

# ---- 8. Summary ----
cat("\n=== Summary ===\n")
cat(sprintf("  Query: %d genes x %d cells\n", nrow(query), ncol(query)))
cat(sprintf("  PhiSpace scores: %d cell types across 4 references\n", ncol(combined_scores)))
cat(sprintf("  Cell types per reference:\n"))
for (lineage in lineages) {
  ct_names <- colnames(sc_list[[lineage]])
  cat(sprintf("    %s: %d types\n", lineage, length(ct_names)))
  for (ct in ct_names) cat(sprintf("      - %s\n", ct))
}

cat("\nDone!\n")
