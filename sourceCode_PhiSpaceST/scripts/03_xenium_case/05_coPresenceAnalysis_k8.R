# coPresenceAnalysis.R
# Cell type co-presence analysis for k=8 clusters
# Uses Spearman correlation to measure co-occurrence of cell types within each niche

suppressPackageStartupMessages({
  library(PhiSpace)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
})

# ---- Paths ----
xenium_dir <- "data/Xenium_NSCLC"
output_dir <- "output"
fig_dir <- "figs/coPresence"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- 1. Load data ----
cat("Loading Xenium SCE...\n")
query <- qread(file.path(xenium_dir, "xenium_lung_sce.qs"))
cat(sprintf("  %d cells, %d cell types\n", ncol(query),
            ncol(reducedDim(query, "PhiSpace"))))

# Use unsmoothed PhiSpace scores
PhiScores <- reducedDim(query, "PhiSpace")
cell_types <- colnames(PhiScores)
n_types <- length(cell_types)

# Get cluster labels (k=8)
clusters <- query$PhiClust_k8
k <- 8
cluster_ids <- sort(unique(clusters))
cat(sprintf("  Using k=%d clusters: %s\n", k, paste(cluster_ids, collapse = ", ")))

# Get pathology for annotation
pathology <- query$pathology_combined

# ---- 2. Compute co-presence matrices (Spearman correlation) ----
cat("\nComputing Spearman co-presence matrices per cluster...\n")

cor_mat_list <- list()
cluster_info <- data.frame(
  cluster = character(),
  n_cells = integer(),
  tumor_pct = numeric(),
  stringsAsFactors = FALSE
)

for (cl in cluster_ids) {

  cells_in_cluster <- which(clusters == cl)
  n_cells <- length(cells_in_cluster)
  tumor_pct <- mean(pathology[cells_in_cluster] == "Tumour") * 100

  cat(sprintf("  Cluster %s: %d cells, %.1f%% Tumor\n", cl, n_cells, tumor_pct))

  # Extract scores for this cluster
  scores <- PhiScores[cells_in_cluster, ]

  # Spearman correlation
  cor_mat <- cor(scores, method = "spearman")
  cor_mat_list[[cl]] <- cor_mat

  cluster_info <- rbind(cluster_info, data.frame(
    cluster = cl,
    n_cells = n_cells,
    tumor_pct = tumor_pct
  ))
}

# ---- 3. Visualize co-presence heatmaps ----
cat("\nGenerating co-presence heatmaps...\n")

# Color scale for correlations
cor_col <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

# Lineage annotation
lineages <- sapply(strsplit(cell_types, "\\("), function(x) {
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
  show_legend = FALSE
)

# Generate heatmap for each cluster
for (cl in cluster_ids) {

  info <- cluster_info[cluster_info$cluster == cl, ]

  ht <- Heatmap(
    cor_mat_list[[cl]],
    name = "Spearman\nrho",
    col = cor_col,
    left_annotation = row_anno,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 5),
    column_title = sprintf("Cluster %s (n=%d, %.0f%% Tumor)",
                           cl, info$n_cells, info$tumor_pct),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 7)
    ),
    width = unit(8, "cm"),
    height = unit(8, "cm")
  )

  png(file.path(fig_dir, sprintf("copresence_cluster%s.png", cl)),
      width = 6, height = 6, units = "in", res = 300)
  draw(ht)
  dev.off()
}
cat(sprintf("  Saved 8 heatmaps to %s/copresence_cluster*.png\n", fig_dir))

# ---- 4. Combined panel of all clusters ----
cat("\nGenerating combined panel...\n")

# Create list of heatmaps without legends for compact display
ht_list <- lapply(cluster_ids, function(cl) {
  info <- cluster_info[cluster_info$cluster == cl, ]
  Heatmap(
    cor_mat_list[[cl]],
    name = sprintf("Cl%s", cl),
    col = cor_col,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    column_title = sprintf("%s (%.0f%%)", cl, info$tumor_pct),
    column_title_gp = gpar(fontsize = 8),
    width = unit(3, "cm"),
    height = unit(3, "cm")
  )
})

# Arrange in 2x4 grid using gridExtra
png(file.path(fig_dir, "copresence_all_clusters.png"),
    width = 14, height = 8, units = "in", res = 300)
pushViewport(viewport(layout = grid.layout(2, 4)))
for (i in seq_along(cluster_ids)) {
  row <- (i - 1) %/% 4 + 1
  col <- (i - 1) %% 4 + 1
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
  draw(ht_list[[i]], newpage = FALSE)
  popViewport()
}
popViewport()
dev.off()
cat(sprintf("  Saved combined panel to %s/copresence_all_clusters.png\n", fig_dir))

# ---- 5. Vectorize and compare across clusters ----
cat("\nVectorizing co-presence matrices for cross-cluster comparison...\n")

# Get upper triangle indices (excluding diagonal)
ut_idx <- which(upper.tri(cor_mat_list[[1]]), arr.ind = TRUE)
n_pairs <- nrow(ut_idx)
cat(sprintf("  %d cell type pairs\n", n_pairs))

# Create pair names
pair_names <- paste(cell_types[ut_idx[, 1]], cell_types[ut_idx[, 2]], sep = " <-> ")

# Vectorize each correlation matrix
copresence_mat <- matrix(NA, nrow = length(cluster_ids), ncol = n_pairs)
rownames(copresence_mat) <- cluster_ids
colnames(copresence_mat) <- pair_names

for (i in seq_along(cluster_ids)) {
  cl <- cluster_ids[i]
  cor_mat <- cor_mat_list[[cl]]
  copresence_mat[i, ] <- cor_mat[upper.tri(cor_mat)]
}

cat(sprintf("  Co-presence matrix: %d clusters x %d pairs\n",
            nrow(copresence_mat), ncol(copresence_mat)))

# ---- 6. PCA on co-presence patterns ----
cat("\nPCA on co-presence patterns...\n")

pca_res <- prcomp(copresence_mat, center = TRUE, scale. = FALSE)

# Create PCA plot
pca_df <- data.frame(
  cluster = rownames(copresence_mat),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  tumor_pct = cluster_info$tumor_pct[match(rownames(copresence_mat), cluster_info$cluster)],
  n_cells = cluster_info$n_cells[match(rownames(copresence_mat), cluster_info$cluster)]
)

# Variance explained
var_explained <- summary(pca_res)$importance[2, ]
var_pc1 <- round(var_explained[1] * 100, 1)
var_pc2 <- round(var_explained[2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = tumor_pct, size = n_cells)) +
  geom_text(aes(label = cluster), vjust = 0, size = 3, fontface = "bold") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 50, name = "Tumor %") +
  scale_size_continuous(name = "N cells", range = c(3, 10)) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_pc1),
    y = sprintf("PC2 (%.1f%%)", var_pc2),
    title = "PCA of niche-specific co-presence patterns (k=8)"
  ) +
  theme_bw(base_size = 6) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(fig_dir, "copresence_PCA.png"), p_pca,
       width = 4, height = 3.5)
cat(sprintf("  Saved PCA plot to %s/copresence_PCA.png\n", fig_dir))

# ---- 7. PCA Loading plots ----
cat("\nGenerating PCA loading plots...\n")

loadings <- pca_res$rotation

for (pc_idx in 1:2) {
  pc_name <- paste0("PC", pc_idx)
  pc_loadings <- loadings[, pc_idx]
  names(pc_loadings) <- pair_names

  # Top 30 by absolute value
  top_idx <- order(abs(pc_loadings), decreasing = TRUE)[1:15]
  top_loadings <- pc_loadings[top_idx]

  load_df <- data.frame(
    pair = gsub("\\([^)]+\\)", "", names(top_loadings)),
    loading = as.numeric(top_loadings),
    stringsAsFactors = FALSE
  )
  load_df$pair <- factor(load_df$pair, levels = load_df$pair[order(load_df$loading)])

  p_load <- ggplot(load_df, aes(x = loading, y = pair, fill = loading)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf("%s Loadings", pc_name)
      # title = sprintf("%s Loadings (%.1f%% var)", pc_name, var_explained[pc_idx] * 100)
    ) +
    theme_bw(base_size = 6) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, sprintf("copresence_%s_loadings.png", pc_name)),
         p_load, width = 3.5, height = 4)
  cat(sprintf("  Saved %s loading plot to %s/copresence_%s_loadings.png\n", pc_name, fig_dir, pc_name))
}

# Save loadings table
loadings_df <- data.frame(
  pair = pair_names,
  PC1 = loadings[, 1],
  PC2 = loadings[, 2]
)
loadings_df <- loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ]
write.csv(loadings_df, file.path(fig_dir, "copresence_PCA_loadings.csv"), row.names = FALSE)
cat(sprintf("  Saved loadings table to %s/copresence_PCA_loadings.csv\n", fig_dir))

# ---- 8. Save results ----
cat("\nSaving results...\n")

results <- list(
  cor_mat_list = cor_mat_list,
  copresence_mat = copresence_mat,
  cluster_info = cluster_info,
  pca_res = pca_res,
  pair_names = pair_names
)

qsave(results, file.path(output_dir, "coPresenceAnalysis.qs"))
cat(sprintf("  Saved to %s/coPresenceAnalysis.qs\n", output_dir))

# ---- Summary ----
cat("\n=== Co-Presence Analysis Complete ===\n")
cat(sprintf("Cluster info:\n"))
print(cluster_info)
cat(sprintf("\nOutputs:\n"))
cat(sprintf("  %s/copresence_cluster*.png (8 individual heatmaps)\n", fig_dir))
cat(sprintf("  %s/copresence_all_clusters.png (combined panel)\n", fig_dir))
cat(sprintf("  %s/copresence_PCA.png\n", fig_dir))
cat(sprintf("  %s/copresence_PC1_loadings.png\n", fig_dir))
cat(sprintf("  %s/copresence_PC2_loadings.png\n", fig_dir))
cat(sprintf("  %s/copresence_PCA_loadings.csv\n", fig_dir))
cat(sprintf("  %s/coPresenceAnalysis.qs\n", output_dir))
cat("\nDone!\n")

