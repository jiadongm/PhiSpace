# clusterEnrichment.R
# Differentially enriched cell types per cluster
# For each k (2..10): heatmap of mean PhiSpace scores + PLS-DA top 5 enriched cell types

suppressPackageStartupMessages({
  library(PhiSpace)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---- Paths ----
xenium_dir <- "data/Xenium_NSCLC"
output_dir <- "output"
fig_dir <- "figs/clusterEnrichment"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- 1. Load data ----
cat("Loading Xenium SCE...\n")
query <- qread(file.path(xenium_dir, "xenium_lung_sce.qs"))
cat(sprintf("  %d cells, %d cell types\n", ncol(query),
            ncol(reducedDim(query, "PhiSpace"))))

# Use unsmoothed PhiSpace scores
PhiScores <- reducedDim(query, "PhiSpace")
cat(sprintf("  Using unsmoothed PhiSpace scores\n"))

# Get pathology labels
pathology <- query$pathology_combined

# ---- 2. Define color schemes ----
# Cluster colors (up to 10)
clust_cols <- c(
"1" = "#E5D8BD", "2" = "#CCEBC5", "3" = "#FED9A6", "4" = "#B3CDE3",
  "5" = "gray",    "6" = "#FBB4AE", "7" = "#DECBE4", "8" = "#FFFFCC",
  "9" = "gray50",  "10" = "#A6CEE3"
)

# Pathology colors (consistent with coPresenceAnalysis_multiK.R)
patho_cols <- c(Tumour = "black", Lymphoid = "#E41A1C",
                `BloodVessels` = "#984EA3", Normal = "#377EB8",
                Unannotated = "gray")

# Heatmap color scale for raw PhiSpace scores (white to red)
# Scores are typically in [-1, 1] range after normPhiScores
heatmap_col <- colorRamp2(c(-0.5, 0, 0.5, 1), c("blue", "white", "orange", "red"))

# ---- 3. Analysis for each k ----
cat("\n=== Cluster Enrichment Analysis ===\n")

plsda_results <- list()

for (k in 10:2) {

  cat(sprintf("\n--- k = %d ---\n", k))

  clust_col <- paste0("PhiClust_k", k)
  clusters <- query[[clust_col]]

  # ----- 3a. Compute mean PhiSpace scores per cluster -----
  cat("  Computing mean scores per cluster...\n")
  cluster_means <- aggregate(
    PhiScores,
    by = list(cluster = clusters),
    FUN = mean
  )
  rownames(cluster_means) <- cluster_means$cluster
  cluster_means <- cluster_means[, -1] %>% as.matrix() %>% t()

  # Order clusters numerically
  cluster_means <- cluster_means[, order(as.numeric(colnames(cluster_means)))]

  # ----- 3b. Compute pathology composition per cluster -----
  path_comp <- table(clusters, pathology)
  path_prop <- prop.table(path_comp, margin = 1)  # row proportions
  ordered_clusters <- as.character(sort(as.numeric(rownames(path_prop))))
  path_prop_mat <- path_prop[ordered_clusters, , drop = FALSE]
  # Ensure all patho categories are present (add zero columns if missing)
  for (cat_name in names(patho_cols)) {
    if (!cat_name %in% colnames(path_prop_mat)) {
      path_prop_mat <- cbind(path_prop_mat, 0)
      colnames(path_prop_mat)[ncol(path_prop_mat)] <- cat_name
    }
  }
  path_prop_mat <- path_prop_mat[, names(patho_cols), drop = FALSE]

  # Cluster sizes
  clust_sizes <- table(clusters)
  clust_sizes <- clust_sizes[order(as.numeric(names(clust_sizes)))]

  # ----- 3c. Create heatmap -----
  cat("  Creating heatmap...\n")

  # Use raw mean scores (not z-scored) to show absolute enrichment levels
  # PhiSpace scores are already normalised, so comparable across cell types

  # Column annotation: stacked bar of pathology proportions
  col_anno <- HeatmapAnnotation(
    Pathology = anno_barplot(
      path_prop_mat,
      gp = gpar(fill = patho_cols[colnames(path_prop_mat)]),
      height = unit(2, "cm"),
      bar_width = 0.8
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 8)
  )

  # Row annotation: lineage
  lineages <- sapply(strsplit(rownames(cluster_means), "\\("), function(x) {
    gsub("\\)", "", x[length(x)])
  })
  lineage_cols <- c(
    "immune" = "#66C2A5", "epithelial" = "#FC8D62",
    "endothelial" = "#8DA0CB", "mesenchymal" = "#E78AC3"
  )
  row_anno <- rowAnnotation(
    Lineage = lineages,
    col = list(Lineage = lineage_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Lineage = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
    )
  )

  # Create heatmap
  ht <- Heatmap(
    cluster_means,
    name = "Mean\nPhiScore",
    col = heatmap_col,
    top_annotation = col_anno,
    left_annotation = row_anno,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title = sprintf("k = %d clusters", k),
    column_title_gp = gpar(fontsize = 6, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    ),
    width = unit(k * 0.8, "cm"),
    row_names_max_width = unit(6, "cm")
  )

  # Save heatmap
  png(file.path(fig_dir, sprintf("heatmap_k%02d.png", k)),
      width = 4, height = 5, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  cat(sprintf("  Saved heatmap to %s/heatmap_k%02d.png\n", fig_dir, k))

  # ----- 3d. PLS-DA for enriched cell types -----
  cat("  Running PLS-DA...\n")

  # Code cluster labels as Y matrix
  Y <- codeY(query, clust_col)
  X <- PhiScores
  ncomp <- min(ncol(Y) - 1, 10)  # cap at 10 components

  # Fit PLS-DA
  pls_res <- mvr(X, Y, ncomp = ncomp, method = "PLS")
  regCoef <- pls_res$coefficients[, , ncomp]

  # Get top 5 enriched cell types per cluster
  top5 <- selectFeat(regCoef, nfeat = 5, absVal = FALSE)
  top5_table <- top5$orderedFeatMat

  # Store results
  plsda_results[[as.character(k)]] <- list(
    regCoef = regCoef,
    top5 = top5_table,
    cluster_means = cluster_means,
    path_prop_mat = path_prop_mat,
    clust_sizes = clust_sizes
  )

  # Save table
  write.csv(top5_table,
            file.path(fig_dir, sprintf("top5_enriched_k%02d.csv", k)),
            row.names = TRUE)
  cat(sprintf("  Saved top 5 table to %s/top5_enriched_k%02d.csv\n", fig_dir, k))

  # Print summary
  cat(sprintf("  Cluster composition:\n"))
  for (cl in colnames(top5_table)) {
    cl_num <- gsub("PhiClust_k\\d+_", "", cl)
    if (cl_num %in% rownames(path_prop_mat)) {
      props <- path_prop_mat[cl_num, ]
      prop_str <- paste(sprintf("%s=%.0f%%", names(props), props * 100), collapse = ", ")
      cat(sprintf("    Cluster %s: %d cells (%s)\n",
                  cl_num, clust_sizes[cl_num], prop_str))
    }
  }
}


# Output k=8 heatmap (reverse x and y)
if(T){
  k = 8
  
  clust_col <- paste0("PhiClust_k", k)
  clusters <- query[[clust_col]]

  cat("  Computing mean scores per cluster...\n")
  cluster_means <- aggregate(
    PhiScores,
    by = list(cluster = clusters),
    FUN = mean
  )
  rownames(cluster_means) <- cluster_means$cluster
  cluster_means <- cluster_means[, -1] %>% as.matrix() %>% t()
  
  # Order clusters numerically
  cluster_means <- cluster_means[, order(as.numeric(colnames(cluster_means)))]
  
  path_comp <- table(clusters, pathology)
  path_prop <- prop.table(path_comp, margin = 1)  # row proportions
  ordered_clusters <- as.character(sort(as.numeric(rownames(path_prop))))
  path_prop_mat <- path_prop[ordered_clusters, , drop = FALSE]
  # Ensure all patho categories are present
  for (cat_name in names(patho_cols)) {
    if (!cat_name %in% colnames(path_prop_mat)) {
      path_prop_mat <- cbind(path_prop_mat, 0)
      colnames(path_prop_mat)[ncol(path_prop_mat)] <- cat_name
    }
  }
  path_prop_mat <- path_prop_mat[, names(patho_cols), drop = FALSE]

  # Cluster sizes
  clust_sizes <- table(clusters)
  clust_sizes <- clust_sizes[order(as.numeric(names(clust_sizes)))]


  # Z-score per cell type (column-scale) so heatmap shows relative enrichment
  # Raw means are small because averaging over thousands of cells washes out signal
  mat <- t(cluster_means)  # rows = clusters, cols = cell types
  mat_scaled <- scale(mat)  # z-score per column (cell type)

  # Strip lineage suffix from cell type names for display
  colnames(mat_scaled) <- gsub("\\(.*\\)$", "", colnames(mat_scaled))

  # Lineage annotation for columns (cell types)
  lineages <- sapply(strsplit(rownames(cluster_means), "\\("), function(x) {
    gsub("\\)", "", x[length(x)])
  })
  lineage_cols <- c(
    "immune" = "#66C2A5", "epithelial" = "#FC8D62",
    "endothelial" = "#8DA0CB", "mesenchymal" = "#E78AC3"
  )
  col_anno <- HeatmapAnnotation(
    Lineage = lineages,
    col = list(Lineage = lineage_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Lineage = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6))
    )
  )

  # Row annotation: stacked bar of pathology proportions
  row_anno <- rowAnnotation(
    Pathology = anno_barplot(
      path_prop_mat,
      gp = gpar(fill = patho_cols[colnames(path_prop_mat)]),
      width = unit(2, "cm"),
      bar_width = 0.8
    ),
    annotation_name_gp = gpar(fontsize = 5)
  )

  # Symmetric color scale for z-scores
  max_z <- max(abs(mat_scaled), na.rm = TRUE)
  zscore_col <- colorRamp2(c(-max_z, 0, max_z), c("blue", "white", "red"))

  # Create heatmap
  ht <- Heatmap(
    mat_scaled,
    name = "Z-scored\nmean PhiScore",
    col = zscore_col,
    top_annotation = col_anno,
    right_annotation = row_anno,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    )
  )

  # Save heatmap
  png(file.path(fig_dir, sprintf("maintext_heatmap_k%02d.png", k)),
      width = 7, height = 3, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
}



# ---- 4. Save all PLS-DA results ----
qsave(plsda_results, file.path(output_dir, "clusterEnrichment_PLSDA.qs"))
cat(sprintf("\nSaved all PLS-DA results to %s/clusterEnrichment_PLSDA.qs\n", output_dir))

# ---- 5. Create summary table (all k values) ----
cat("\nCreating summary table...\n")

summary_list <- lapply(names(plsda_results), function(k) {
  res <- plsda_results[[k]]
  top5 <- res$top5

  # For each cluster, get top enriched cell type
  top1 <- apply(top5, 2, function(col) {
    col[1]  # top 1
  })

  cl_nums <- gsub(".*_", "", names(top1))
  df <- data.frame(
    k = as.integer(k),
    cluster = names(top1),
    n_cells = as.numeric(res$clust_sizes[cl_nums]),
    top1_celltype = top1,
    stringsAsFactors = FALSE
  )
  # Add pathology proportions
  for (cat_name in colnames(res$path_prop_mat)) {
    df[[cat_name]] <- res$path_prop_mat[cl_nums, cat_name]
  }
  df
})

summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df[order(-summary_df$k, summary_df$cluster), ]
rownames(summary_df) <- NULL

write.csv(summary_df, file.path(fig_dir, "summary_all_k.csv"), row.names = FALSE)
cat(sprintf("Saved summary to %s/summary_all_k.csv\n", fig_dir))

# ---- Summary ----
cat("\n=== Cluster Enrichment Analysis Complete ===\n")
cat(sprintf("Figures saved: %s/heatmap_k{02..10}.png\n", fig_dir))
cat(sprintf("Tables saved: %s/top5_enriched_k{02..10}.csv\n", fig_dir))
cat(sprintf("Summary: %s/summary_all_k.csv\n", fig_dir))
cat("\nDone!\n")
