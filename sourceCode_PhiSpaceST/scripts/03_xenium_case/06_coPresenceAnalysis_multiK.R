# coPresenceAnalysis_multiK.R
# Cell type co-presence analysis for k=5, 6, and 8 clusters
# Uses Spearman correlation to measure co-occurrence of cell types within each niche

suppressPackageStartupMessages({
  library(PhiSpace)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(scatterpie)
})

# ---- Paths ----
xenium_dir <- "data/Xenium_NSCLC"
output_dir <- "output"
fig_dir_base <- "figs/coPresence"
if (!dir.exists(fig_dir_base)) dir.create(fig_dir_base, recursive = TRUE)

# ---- 1. Load data ----
cat("Loading Xenium SCE...\n")
query <- qread(file.path(xenium_dir, "xenium_lung_sce.qs"))
cat(sprintf("  %d cells, %d cell types\n", ncol(query),
            ncol(reducedDim(query, "PhiSpace"))))

# Use unsmoothed PhiSpace scores
PhiScores <- reducedDim(query, "PhiSpace")
cell_types <- colnames(PhiScores)
n_types <- length(cell_types)

# Get pathology
pathology <- query$pathology_combined

# Color scales
cor_col <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

# Lineage annotation
lineages <- sapply(strsplit(cell_types, "\\("), function(x) {
  gsub("\\)", "", x[length(x)])
})
lineage_cols <- c(
  "immune" = "#66C2A5", "epithelial" = "#FC8D62",
  "endothelial" = "#8DA0CB", "mesenchymal" = "#E78AC3"
)

# Pathology colors (consistent with clusterEnrichment.R)
patho_cols <- c(Tumour = "black", Lymphoid = "#E41A1C",
                `BloodVessels` = "#984EA3", Normal = "#377EB8",
                Unannotated = "gray")

# ---- 2. Run analysis for k = 4, 5, 6, 7, 8, 9, 10 ----
k_values <- c(4, 5, 6, 7, 8, 9, 10)
all_results <- list()

for (k in k_values) {

  cat(sprintf("\n========== k = %d ==========\n", k))

  fig_dir <- file.path(fig_dir_base, sprintf("k%d", k))
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  # Get cluster labels
  clust_col <- paste0("PhiClust_k", k)
  clusters <- query[[clust_col]]
  cluster_ids <- sort(unique(clusters))
  cat(sprintf("  Clusters: %s\n", paste(cluster_ids, collapse = ", ")))

  # ---- 2a. Compute co-presence matrices ----
  cat("  Computing Spearman co-presence matrices...\n")

  cor_mat_list <- list()

  # Compute full pathology proportions per cluster
  path_comp <- table(clusters, pathology)
  path_prop <- prop.table(path_comp, margin = 1)

  cluster_info <- data.frame(
    cluster = character(),
    n_cells = integer(),
    stringsAsFactors = FALSE
  )

  for (cl in cluster_ids) {
    cells_in_cluster <- which(clusters == cl)
    n_cells <- length(cells_in_cluster)

    # Build row with pathology proportions
    info_row <- data.frame(cluster = cl, n_cells = n_cells, stringsAsFactors = FALSE)
    for (cat_name in names(patho_cols)) {
      if (cat_name %in% colnames(path_prop) && cl %in% rownames(path_prop)) {
        info_row[[cat_name]] <- path_prop[cl, cat_name]
      } else {
        info_row[[cat_name]] <- 0
      }
    }

    prop_str <- paste(sprintf("%s=%.0f%%", names(patho_cols),
                              sapply(names(patho_cols), function(x) info_row[[x]] * 100)),
                      collapse = ", ")
    cat(sprintf("    Cluster %s: %d cells (%s)\n", cl, n_cells, prop_str))

    scores <- PhiScores[cells_in_cluster, ]
    cor_mat <- cor(scores, method = "spearman")
    cor_mat_list[[cl]] <- cor_mat

    cluster_info <- rbind(cluster_info, info_row)
  }

  # ---- 2b. Individual heatmaps ----
  cat("  Generating heatmaps...\n")

  row_anno <- rowAnnotation(
    Lineage = lineages,
    col = list(Lineage = lineage_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

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
                             cl, info$n_cells, info$Tumor * 100),
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

  # ---- 2c. Combined panel ----
  n_clusters <- length(cluster_ids)
  ncol_grid <- min(4, n_clusters)
  nrow_grid <- ceiling(n_clusters / ncol_grid)

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
      column_title = sprintf("%s (%.0f%%T)", cl, info$Tumor * 100),
      column_title_gp = gpar(fontsize = 8),
      width = unit(3, "cm"),
      height = unit(3, "cm")
    )
  })

  png(file.path(fig_dir, "copresence_all_clusters.png"),
      width = ncol_grid * 3.5, height = nrow_grid * 4, units = "in", res = 300)
  pushViewport(viewport(layout = grid.layout(nrow_grid, ncol_grid)))
  for (i in seq_along(cluster_ids)) {
    row <- (i - 1) %/% ncol_grid + 1
    col <- (i - 1) %% ncol_grid + 1
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    draw(ht_list[[i]], newpage = FALSE)
    popViewport()
  }
  popViewport()
  dev.off()

  # ---- 2d. Vectorize and PCA ----
  cat("  PCA on co-presence patterns...\n")

  ut_idx <- which(upper.tri(cor_mat_list[[1]]), arr.ind = TRUE)
  n_pairs <- nrow(ut_idx)
  pair_names <- paste(cell_types[ut_idx[, 1]], cell_types[ut_idx[, 2]], sep = " <-> ")

  copresence_mat <- matrix(NA, nrow = length(cluster_ids), ncol = n_pairs)
  rownames(copresence_mat) <- cluster_ids
  colnames(copresence_mat) <- pair_names

  for (i in seq_along(cluster_ids)) {
    cl <- cluster_ids[i]
    cor_mat <- cor_mat_list[[cl]]
    copresence_mat[i, ] <- cor_mat[upper.tri(cor_mat)]
  }

  # PCA
  if (nrow(copresence_mat) > 2) {
    pca_res <- prcomp(copresence_mat, center = TRUE, scale. = FALSE)

    pca_df <- data.frame(
      cluster = rownames(copresence_mat),
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      n_cells = cluster_info$n_cells[match(rownames(copresence_mat), cluster_info$cluster)],
      stringsAsFactors = FALSE
    )
    # Add pathology proportion columns
    for (cat_name in names(patho_cols)) {
      pca_df[[cat_name]] <- cluster_info[[cat_name]][match(pca_df$cluster, cluster_info$cluster)]
    }

    var_explained <- summary(pca_res)$importance[2, ]
    var_pc1 <- round(var_explained[1] * 100, 1)
    var_pc2 <- round(var_explained[2] * 100, 1)

    # Pie radius proportional to sqrt(n_cells) for area-proportionality
    pie_scale <- diff(range(pca_df$PC1)) * 0.08
    pca_df$radius <- sqrt(pca_df$n_cells) / max(sqrt(pca_df$n_cells)) * pie_scale

    p_pca <- ggplot() +
      geom_scatterpie(
        aes(x = PC1, y = PC2, r = radius),
        data = pca_df,
        cols = names(patho_cols),
        color = NA
      ) +
      scale_fill_manual(values = patho_cols, name = "Pathology") +
      geom_text(data = pca_df, aes(x = PC1, y = PC2, label = cluster),
                vjust = -1, size = 3, fontface = "bold") +
      coord_equal() +
      labs(
        x = sprintf("PC1 (%.1f%%)", var_pc1),
        y = sprintf("PC2 (%.1f%%)", var_pc2),
        title = sprintf("PCA of co-presence patterns (k=%d)", k)
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(fig_dir, "copresence_PCA.png"), p_pca, width = 7, height = 5)

    # ---- PCA Loadings plots ----
    cat("  Generating PCA loading plots...\n")

    # Get loadings (rotation matrix)
    loadings <- pca_res$rotation

    # For PC1 and PC2, get top 15 positive and negative loadings
    for (pc_idx in 1:min(2, ncol(loadings))) {
      pc_name <- paste0("PC", pc_idx)
      pc_loadings <- loadings[, pc_idx]
      names(pc_loadings) <- pair_names

      # Sort by absolute value and get top 30
      top_idx <- order(abs(pc_loadings), decreasing = TRUE)[1:30]
      top_loadings <- pc_loadings[top_idx]

      # Create data frame for plotting
      load_df <- data.frame(
        pair = names(top_loadings),
        loading = as.numeric(top_loadings),
        stringsAsFactors = FALSE
      )
      load_df$pair <- factor(load_df$pair, levels = load_df$pair[order(load_df$loading)])

      p_load <- ggplot(load_df, aes(x = loading, y = pair, fill = loading)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        labs(
          x = sprintf("%s Loading", pc_name),
          y = "Cell type pair",
          title = sprintf("%s Loadings (k=%d, %.1f%% var)", pc_name, k, var_explained[pc_idx] * 100)
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )

      ggsave(file.path(fig_dir, sprintf("copresence_%s_loadings.png", pc_name)),
             p_load, width = 8, height = 8)
    }

    # Save loadings table
    loadings_df <- data.frame(
      pair = pair_names,
      PC1 = loadings[, 1],
      PC2 = loadings[, 2]
    )
    loadings_df <- loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ]
    write.csv(loadings_df, file.path(fig_dir, "copresence_PCA_loadings.csv"), row.names = FALSE)

  } else {
    pca_res <- NULL
  }

  # ---- 2e. Discriminant pairs ----
  cat("  Identifying discriminant co-presence pairs...\n")

  high_tumor_clusters <- cluster_info$cluster[cluster_info$Tumor > 0.90]
  low_tumor_clusters <- cluster_info$cluster[cluster_info$Tumor < 0.50]

  cat(sprintf("    High tumor (>90%%): %s\n",
              ifelse(length(high_tumor_clusters) > 0, paste(high_tumor_clusters, collapse = ", "), "none")))
  cat(sprintf("    Low tumor (<50%%): %s\n",
              ifelse(length(low_tumor_clusters) > 0, paste(low_tumor_clusters, collapse = ", "), "none")))

  discrim_df <- NULL
  if (length(high_tumor_clusters) > 0 && length(low_tumor_clusters) > 0) {
    high_mean <- colMeans(copresence_mat[high_tumor_clusters, , drop = FALSE])
    low_mean <- colMeans(copresence_mat[low_tumor_clusters, , drop = FALSE])
    diff_copresence <- high_mean - low_mean

    top_high <- sort(diff_copresence, decreasing = TRUE)[1:20]
    top_low <- sort(diff_copresence, decreasing = FALSE)[1:20]

    discrim_df <- data.frame(
      pair = c(names(top_high), names(top_low)),
      diff = c(top_high, top_low),
      direction = c(rep("High tumor", 20), rep("Low tumor", 20)),
      high_tumor_mean = c(high_mean[names(top_high)], high_mean[names(top_low)]),
      low_tumor_mean = c(low_mean[names(top_high)], low_mean[names(top_low)])
    )

    write.csv(discrim_df, file.path(fig_dir, "discriminant_pairs.csv"), row.names = FALSE)

    cat("\n    Top 5 pairs in HIGH tumor clusters:\n")
    for (i in 1:5) {
      cat(sprintf("      %s: diff=%.3f\n", names(top_high)[i], top_high[i]))
    }
    cat("    Top 5 pairs in LOW tumor clusters:\n")
    for (i in 1:5) {
      cat(sprintf("      %s: diff=%.3f\n", names(top_low)[i], top_low[i]))
    }
  }

  # Store results
  all_results[[as.character(k)]] <- list(
    cor_mat_list = cor_mat_list,
    copresence_mat = copresence_mat,
    cluster_info = cluster_info,
    pca_res = pca_res,
    discrim_df = discrim_df,
    pair_names = pair_names
  )
}

# ---- 3. Save all results ----
cat("\n\nSaving all results...\n")
qsave(all_results, file.path(output_dir, "coPresenceAnalysis_multiK.qs"))
cat(sprintf("  Saved to %s/coPresenceAnalysis_multiK.qs\n", output_dir))

# ---- 4. Summary comparison across k values ----
cat("\n=== Summary: Cluster Composition Across k ===\n")
for (k in k_values) {
  cat(sprintf("\nk = %d:\n", k))
  info <- all_results[[as.character(k)]]$cluster_info
  info <- info[order(-info$Tumor), ]
  for (i in 1:nrow(info)) {
    prop_str <- paste(sprintf("%s=%.0f%%", names(patho_cols),
                              sapply(names(patho_cols), function(x) info[[x]][i] * 100)),
                      collapse = ", ")
    cat(sprintf("  Cluster %s: %5d cells (%s)\n",
                info$cluster[i], info$n_cells[i], prop_str))
  }
}

cat("\n=== Co-Presence Analysis Complete ===\n")
cat("Outputs per k (k=4..10):\n")
cat("  figs/coPresence/k{4..10}/copresence_cluster*.png\n")
cat("  figs/coPresence/k{4..10}/copresence_all_clusters.png\n")
cat("  figs/coPresence/k{4..10}/copresence_PCA.png\n")
cat("  figs/coPresence/k{4..10}/copresence_PC1_loadings.png\n")
cat("  figs/coPresence/k{4..10}/copresence_PC2_loadings.png\n")
cat("  figs/coPresence/k{4..10}/copresence_PCA_loadings.csv\n")
cat("  figs/coPresence/k{4..10}/discriminant_pairs.csv\n")
cat("  output/coPresenceAnalysis_multiK.qs\n")
cat("\nDone!\n")
