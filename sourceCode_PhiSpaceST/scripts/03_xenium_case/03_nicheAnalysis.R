# nicheAnalysis.R
# Niche analysis for Xenium lung cancer data
# Step 5: Spatial smoothing, PCA, multi-k clustering, gene-expression comparison

suppressPackageStartupMessages({
  library(PhiSpace)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
})

# ---- Paths ----
xenium_dir <- "data/Xenium_NSCLC"
output_dir <- "output"
fig_dir <- "figs"
for (d in c(output_dir, fig_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- 1. Load data ----
cat("Loading Xenium SCE with PhiSpace scores...\n")
query <- qread(file.path(xenium_dir, "xenium_lung_sce.qs"))
cat(sprintf("  %d genes x %d cells\n", nrow(query), ncol(query)))
cat(sprintf("  PhiSpace scores: %d cell types\n",
            ncol(reducedDim(query, "PhiSpace"))))

query$y_centroid_flipped <- - query$y_centroid

x_coord <- "x_centroid"
y_coord <- "y_centroid_flipped"

# ---- 2. Spatial smoothing ----
cat("\nSpatial smoothing of PhiSpace scores...\n")
query <- spatialSmoother(
  query,
  reducedDim2smooth = "PhiSpace",
  smoothReducedDim = TRUE,
  x_coord = x_coord,
  y_coord = y_coord,
  k = 15,
  kernel = "gaussian"
)
cat(sprintf("  Created 'PhiSpace_smoothed': %d x %d\n",
            nrow(reducedDim(query, "PhiSpace_smoothed")),
            ncol(reducedDim(query, "PhiSpace_smoothed"))))

# ---- 3. PCA on smoothed scores ----
cat("\nPCA on smoothed PhiSpace scores...\n")
nPhiTypes <- ncol(reducedDim(query, "PhiSpace_smoothed"))
PhiPCRes <- getPC(
  reducedDim(query, "PhiSpace_smoothed"),
  ncomp = nPhiTypes - 1
)

# Elbow plot
png(file.path(fig_dir, "elbow_plot_smoothed.png"),
    width = 5, height = 4, units = "in", res = 300)
plot(1 - PhiPCRes$accuProps,
     xlab = "Number of PCs",
     ylab = "Residual variance proportion",
     main = "Elbow plot (smoothed PhiSpace scores)",
     pch = 16, cex = 0.6)
dev.off()
cat("  Saved elbow plot to figs/elbow_plot_smoothed.png\n")

# Select number of PCs (use 25 or nPhiTypes-1 if fewer)
nPC <- min(30, nPhiTypes - 1)
cat(sprintf("  Using %d PCs for clustering\n", nPC))
mat2clust <- PhiPCRes$scores[, 1:nPC]

# ---- 4. Multi-k clustering (k=2..10) with alignment ----
cat("\nK-means clustering with k = 2..10...\n")
clustResPath <- file.path(output_dir, "nicheClusterResults.qs")
if (!file.exists(clustResPath)) {

  set.seed(94863)
  outClusts <- lapply(2:10, function(kclust) {
    cat(sprintf("  k = %d...\n", kclust))
    kmeans(mat2clust, centers = kclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")
  })

  qsave(list(
    outClusts = outClusts,
    PhiPCRes = PhiPCRes,
    nPC = nPC
  ), clustResPath)
  cat(sprintf("  Saved clustering results to %s\n", clustResPath))
} else {
  cat(sprintf("  Loading existing results from %s\n", clustResPath))
  clustRes <- qread(clustResPath)
  outClusts <- clustRes$outClusts
}

# Extract cluster factors
clust_list <- lapply(outClusts, function(x) {
  factor(x$cluster, levels = sort(unique(x$cluster)))
})

# Align clusters backwards (k=10 down to k=2) using Hungarian algorithm
cat("  Aligning clusters across k values...\n")
for (x in length(clust_list):1) {
  clust <- clust_list[[x]]
  if (x < length(clust_list)) {
    clust_old <- clust_list[[x + 1]]
    clust <- vizOmics:::align_clusters(clust, clust_old) %>% as.factor()
    clust_list[[x]] <- factor(clust,
                              levels = as.character(sort(as.numeric(levels(clust)))))
  }
}

# Store all clustering results in colData
for (ii in seq_along(clust_list)) {
  k_val <- ii + 1  # k=2..10
  col_name <- paste0("PhiClust_k", k_val)
  query[[col_name]] <- as.character(clust_list[[ii]])
  cat(sprintf("  Added colData$%s\n", col_name))
}

# ---- 5. Visualise niches for all k values ----
cat("\nGenerating spatial niche plots for k = 2..10...\n")
clust_cols_10 <- c(
  "1" = "#E5D8BD", "2" = "#CCEBC5", "3" = "#FED9A6", "4" = "#B3CDE3",
  "5" = "gray",    "6" = "#FBB4AE", "7" = "#DECBE4", "8" = "#FFFFCC",
  "9" = "gray50",  "10" = "#A6CEE3"
)

outPlots <- vector("list", length(clust_list))
for (ii in seq_along(clust_list)) {
  k_val <- ii + 1
  col_name <- paste0("PhiClust_k", k_val)
  outPlots[[ii]] <- VizSpatial(
    query, x_coord = x_coord, y_coord = y_coord,
    colBy = col_name, ptSize = 0.3
  ) +
    scale_colour_manual(values = clust_cols_10) +
    guides(colour = guide_legend(
      override.aes = list(size = 2), nrow = 1
    )) +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.spacing = unit(0, "pt")
    ) +
    ggtitle(paste0("k = ", k_val))
}

p_all <- ggarrange(plotlist = outPlots, nrow = 3, ncol = 3, legend = "none")
ggsave(
  file.path(fig_dir, "niche_clusters_k2_to_k10.png"),
  p_all, width = 12, height = 8
)
cat("  Saved figs/niche_clusters_k2_to_k10.png\n")

# Also save individual legend for reference (k=8 as example)
p_legend <- VizSpatial(
  query, x_coord = x_coord, y_coord = y_coord,
  colBy = "PhiClust_k8", ptSize = 0.3
) +
  scale_colour_manual(values = clust_cols_10) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
  theme(legend.position = "top", legend.key.spacing = unit(0, "pt"))

ggsave(file.path(fig_dir, "niche_clusters_k=8.png"), 
       plot = p_legend + theme(legend.position = "none"), 
       width = 3.5, height = 2)

png(file.path(fig_dir, "niche_clusters_legend.png"),
    width = 6, height = 0.5, units = "in", res = 300)
op <- par(mar = rep(0, 4))
grid::grid.draw(ggpubr::get_legend(p_legend))
par(op)
dev.off()
cat("  Saved figs/niche_clusters_legend.png\n")

# pathology_combined is defined in runPhiSpace4Ref.R
p_patho <- VizSpatial(
  query, x_coord = x_coord, y_coord = y_coord,
  colBy = "pathology_combined", ptSize = 0.3
) +
  scale_colour_manual(values = c(Tumour = "black", Lymphoid = "#E41A1C",
                                   `BloodVessels` = "#984EA3", Normal = "#377EB8",
                                   Unannotated = "gray")) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
  theme(legend.position = "top", legend.key.spacing = unit(0, "pt"))
ggsave(file.path(fig_dir, "patho.png"), 
       plot = p_patho + theme(legend.position = "none"), 
       width = 3.5, height = 2)

png(file.path(fig_dir, "patho_legend.png"),
    width = 8, height = 0.5, units = "in", res = 300)
op <- par(mar = rep(0, 4))
grid::grid.draw(ggpubr::get_legend(p_patho))
par(op)
dev.off()

# ---- 6. Gene-expression-based clustering (comparison) ----
cat("\nGene-expression-based clustering for comparison...\n")
genClustResPath <- file.path(output_dir, "genExprClusterResults.qs")
if (!file.exists(genClustResPath)) {

  cat("  PCA on log1p gene expression...\n")
  GenPCRes <- getPC(t(assay(query, "log1p")), ncomp = 20)

  # Use k=8 as default (can be adjusted after inspecting multi-k plots)
  k_chosen <- 8
  cat(sprintf("  K-means with k = %d on gene expression PCs...\n", k_chosen))
  set.seed(94858)
  genClust_res <- kmeans(
    GenPCRes$scores, centers = k_chosen, iter.max = 500L, nstart = 50, algorithm = "Lloyd"
  )

  qsave(list(
    genClust_res = genClust_res,
    GenPCRes = GenPCRes,
    k_chosen = k_chosen
  ), genClustResPath)
  cat(sprintf("  Saved gene-expression clustering to %s\n", genClustResPath))
} else {
  cat(sprintf("  Loading existing results from %s\n", genClustResPath))
  genClustRes <- qread(genClustResPath)
  genClust_res <- genClustRes$genClust_res
  k_chosen <- genClustRes$k_chosen
}

# Align GenClust labels to PhiClust labels
phi_clust_col <- paste0("PhiClust_k", k_chosen)
query$GenClust <- vizOmics::align_clusters(
  as.character(genClust_res$cluster),
  query[[phi_clust_col]]
)
cat(sprintf("  Aligned GenClust to %s\n", phi_clust_col))

# ---- 7. Side-by-side comparison plot ----
cat("\nGenerating PhiClust vs GenClust comparison plot...\n")
p_phi <- VizSpatial(
  query, x_coord = x_coord, y_coord = y_coord,
  colBy = phi_clust_col, ptSize = 0.3
) +
  scale_colour_manual(values = clust_cols_10) +
  guides(colour = guide_legend(
    override.aes = list(size = 2), nrow = 1
  )) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.spacing = unit(0, "pt")
  ) +
  ggtitle("PhiSpace niches")

p_gen <- VizSpatial(
  query, x_coord = x_coord, y_coord = y_coord,
  colBy = "GenClust", ptSize = 0.3
) +
  scale_colour_manual(values = clust_cols_10) +
  guides(colour = guide_legend(
    override.aes = list(size = 2), nrow = 1
  )) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.spacing = unit(0, "pt")
  ) +
  ggtitle("Gene expression clusters")

p_compare <- ggarrange(p_phi, p_gen, ncol = 2, common.legend = TRUE, legend = "top")
ggsave(
  file.path(fig_dir, "PhiClust_vs_GenClust.png"),
  p_compare, width = 10, height = 5
)
cat("  Saved figs/PhiClust_vs_GenClust.png\n")

# ---- 8. Save updated SCE ----
cat("\nSaving updated SCE with cluster columns...\n")
qsave(query, file.path(xenium_dir, "xenium_lung_sce.qs"))
cat("  Saved to data/Xenium_NSCLC/xenium_lung_sce.qs\n")

# ---- Summary ----
cat("\n=== Niche Analysis Summary ===\n")
cat(sprintf("  Spatial smoothing: k=15, gaussian kernel\n"))
cat(sprintf("  PCA: %d PCs (of %d total)\n", nPC, nPhiTypes - 1))
cat(sprintf("  K-means: k=2..10 (PhiSpace), k=%d (gene expression)\n", k_chosen))
cat(sprintf("  Cluster columns added: PhiClust_k2..PhiClust_k10, GenClust\n"))
cat(sprintf("  Pathology breakdown:\n"))
print(table(query$pathology))
cat("\nFigures saved:\n")
cat("  figs/elbow_plot_smoothed.png\n")
cat("  figs/niche_clusters_k2_to_k10.png\n")
cat("  figs/niche_clusters_legend.png\n")
cat("  figs/PhiClust_vs_GenClust.png\n")
cat("\nDone!\n")
