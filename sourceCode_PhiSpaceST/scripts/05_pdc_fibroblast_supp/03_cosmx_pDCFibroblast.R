## ============================================================================
## pDC-Fibroblast co-presence validation in CosMx NSCLC
## Independent cohort validation for Reviewer 2, Comment 3
## ============================================================================

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(PhiSpace)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(ggpubr)
  library(qs)
  library(ComplexHeatmap)
})

source("../CaseCosMx/utils.R")
dat_dir <- paths$phispace_data_root
PhiAssay <- "log1p"

# =============================================================================
# 1. Load data
# =============================================================================

tissueNames <- tissueNames_CosMx

# Query SCE list
query_list <- qread(paste0(dat_dir, "data/CosMx_lung/SCE_obj/allTissue_SCE.qs"))

# PhiSpace results (4 lineage-specific references)
PhiResPath <- paste0(dat_dir, "output/Case3/CosMxAllLungsPhiRes4Refs_", PhiAssay, ".qs")
sc_list <- qread(PhiResPath)

# Combine scores across 4 lineages and simplify names (same as CaseCosMx_4Refs.Rmd)
for (ii in seq_along(tissueNames)) {
  tissueName <- tissueNames[ii]
  query <- query_list[[tissueName]]
  sc4Refs <- lapply(sc_list, function(x) x[[tissueName]])
  reducedDim(query, "PhiSpace") <- Reduce(cbind, sc4Refs)
  query_list[[ii]] <- query
}

# Simplify cell type names (same as CaseCosMx_4Refs.Rmd lines 74-86)
simpleNames <- colnames(reducedDim(query_list[[1]], "PhiSpace"))
simpleNames[c(2, 4, 9, 10, 12, 18, 21, 23, 24, 25, 26, 43)] <- c(
  "MoMacroph(immune)", "InflamMono(immune)", "AlveolarMacroph(immune)",
  "ProlifImm(immune)", "InterstMacroph(immune)", "ProlifEpi(epithelial)",
  "SecretSCGB3A2+(epithelial)", "TransAT2(epithelial)", "SecretSCGB1A1+/MUC5B+(epithelial)",
  "SecretSCGB1A1+/SCGB3A2+(epithelial)", "DiffCiliated(epithelial)", "MyoFB Act(mesenchymal)"
)
for (ii in seq_along(query_list)) {
  colnames(reducedDim(query_list[[ii]], "PhiSpace")) <- simpleNames
  query_list[[ii]] <- query_list[[ii]]
}

# Identify pDC and fibroblast columns
cTypeNames <- simpleNames
pDC_col <- grep("pDC", cTypeNames, value = TRUE)
FB_cols <- grep("FB|MyoFB", cTypeNames, value = TRUE)
cat("pDC column:", pDC_col, "\n")
cat("Fibroblast columns:", paste(FB_cols, collapse = ", "), "\n")

# Short labels for fibroblast subtypes
FB_shortNames <- gsub("\\(mesenchymal\\)", "", FB_cols)
names(FB_shortNames) <- FB_cols

# Minimum cells per niche to compute correlation
minNicheCells <- 50

# =============================================================================
# 1b. Spatial smoothing of PhiSpace scores
# =============================================================================
# CosMx is single-cell resolution; Visium spots are multicellular (~10-20 cells).
# To make co-presence comparable, smooth PhiSpace scores over spatial neighbours.
# Smoothed scores approximate the average score in a local neighbourhood,
# analogous to a Visium spot.

smoothed_path <- file.path(
  file.path(paths$output_root, "pDC"),
  "output", "cosmx_query_list_smoothed.qs"
)

if (file.exists(smoothed_path)) {
  cat("Loading pre-computed smoothed scores...\n")
  query_list <- qread(smoothed_path)
} else {
  cat("Normalising then spatially smoothing PhiSpace scores (k=10, uniform kernel)...\n")
  for (ii in seq_along(tissueNames)) {
    tissueName <- tissueNames[ii]
    cat(sprintf("  [%s] %d cells\n", tissueName, ncol(query_list[[tissueName]])))
    # Normalise first, then smooth
    reducedDim(query_list[[tissueName]], "PhiSpace_norm") <-
      normPhiScores(reducedDim(query_list[[tissueName]], "PhiSpace"))
    query_list[[tissueName]] <- spatialSmoother(
      query_list[[tissueName]],
      reducedDim2smooth = "PhiSpace_norm",
      smoothReducedDim = TRUE,
      x_coord = "sdimx",
      y_coord = "sdimy",
      k = 10,
      kernel = "uniform"
    )
  }
  qsave(query_list, smoothed_path)
  cat("Saved smoothed query list.\n")
}

# Use normalised-then-smoothed scores; no further normPhiScores needed
scoreDimName <- "PhiSpace_norm_smoothed"


# =============================================================================
# 2. Spatial heatmaps: representative sample
# =============================================================================

repSample <- "Lung5_Rep1"
query <- query_list[[repSample]]
scores <- reducedDim(query, scoreDimName)

spatialHeatmap <- function(sce, scores, col_name, short_title = col_name,
                           coordNames = c("sdimx", "sdimy"), ptSize = 0.1) {
  df <- as.data.frame(colData(sce)) %>%
    mutate(score = scores[, col_name]) %>%
    arrange(score)

  ggplot(df, aes(x = !!sym(coordNames[1]), y = !!sym(coordNames[2]))) +
    geom_point(aes(colour = score), size = ptSize, stroke = 0) +
    scale_colour_gradientn(colours = MATLAB_cols) +
    ggtitle(short_title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      legend.position = "none"
    )
}

panel_list <- list()

# Niche map (full niche labels)
niche_df <- as.data.frame(colData(query))
panel_list[["niche"]] <- ggplot(niche_df, aes(x = sdimx, y = sdimy)) +
  geom_point(aes(colour = niche), size = 0.1, stroke = 0) +
  scale_colour_manual(values = nicheCol) +
  ggtitle(paste0("Niches (", repSample, ")")) +
  theme_void() +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    legend.position = "none"
  )

# pDC
panel_list[["pDC"]] <- spatialHeatmap(query, scores, pDC_col, short_title = "pDC")

# Fibroblast subtypes
for (fb in FB_cols) {
  panel_list[[fb]] <- spatialHeatmap(query, scores, fb, short_title = FB_shortNames[fb])
}

# Co-localization
pDC_scores <- scores[, pDC_col]
FB_scores_mat <- scores[, FB_cols, drop = FALSE]
max_FB <- apply(FB_scores_mat, 1, max)
coloc <- pDC_scores * max_FB
coloc_mat <- matrix(coloc, ncol = 1, dimnames = list(names(coloc), "colocalization"))
panel_list[["coloc"]] <- spatialHeatmap(query, coloc_mat, "colocalization", short_title = "pDC-FB co-loc")

ncols <- length(panel_list)
p_spatial <- ggarrange(plotlist = panel_list, nrow = 1, ncol = ncols)
ggsave(
  "figs/cosmx_spatial_heatmaps.png", p_spatial,
  width = 2 * ncols, height = 2, dpi = 300
)
cat("Saved: figs/cosmx_spatial_heatmaps.png\n")


# =============================================================================
# 3. Per-sample, per-niche pDC-fibroblast Spearman correlations
# =============================================================================

cor_results <- data.frame(
  sample = character(),
  niche = character(),
  fibroblast = character(),
  correlation = numeric(),
  n_cells = integer(),
  stringsAsFactors = FALSE
)

for (ii in seq_along(tissueNames)) {
  tissueName <- tissueNames[ii]
  query <- query_list[[tissueName]]
  scores <- reducedDim(query, scoreDimName)
  niches <- query$niche
  nicheTypes <- unique(niches)

  for (nch in nicheTypes) {
    idx <- niches == nch
    n_cells <- sum(idx)
    if (n_cells < minNicheCells) next

    pDC_sc <- scores[idx, pDC_col]

    for (fb in FB_cols) {
      fb_sc <- scores[idx, fb]
      rho <- cor(pDC_sc, fb_sc, method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        sample = tissueName,
        niche = nch,
        fibroblast = FB_shortNames[fb],
        correlation = rho,
        n_cells = n_cells,
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(cor_results, "output/cosmx_pDC_fb_correlations.csv", row.names = FALSE)
cat("Saved: output/cosmx_pDC_fb_correlations.csv\n")

# Niche ordering for plots
niche_order <- c(
  "tumor interior", "tumor-stroma boundary",
  "stroma", "immune", "lymphoid structure",
  "myeloid-enriched stroma", "neutrophils",
  "plasmablast-enriched stroma", "macrophages"
)
niche_order <- intersect(niche_order, unique(cor_results$niche))
cor_results$niche <- factor(cor_results$niche, levels = niche_order)


# =============================================================================
# 4. Boxplot: per-niche pDC-fibroblast correlations
# =============================================================================

p_box <- ggplot(cor_results, aes(x = niche, y = correlation, fill = niche)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  facet_wrap(~ fibroblast, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = nicheCol) +
  theme_pubr(base_size = 7) +
  ylab("Spearman correlation with pDC") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, size = 5),
    strip.text = element_text(size = 6)
  )
ggsave(
  "figs/cosmx_correlation_boxplot.png", p_box,
  width = 2 * length(FB_cols), height = 3.5, dpi = 300
)
cat("Saved: figs/cosmx_correlation_boxplot.png\n")


# =============================================================================
# 5. Statistical tests: Kruskal-Wallis across niches
# =============================================================================

cat("\n=== Statistical tests (across niches, CosMx) ===\n")
stat_results <- list()
for (fb in unique(cor_results$fibroblast)) {
  sub <- cor_results %>% filter(fibroblast == fb)
  kw <- kruskal.test(correlation ~ niche, data = sub)
  stat_results[[fb]] <- list(kruskal_wallis_p = kw$p.value)
  cat(sprintf("  %s: KW p = %.4f\n", fb, kw$p.value))
}


# =============================================================================
# 6. Correlation heatmap: (sample x niche) rows, fibroblast subtypes columns
# =============================================================================

cor_wide <- cor_results %>%
  pivot_wider(names_from = fibroblast, values_from = correlation) %>%
  as.data.frame()
rownames(cor_wide) <- paste0(cor_wide$sample, " — ", cor_wide$niche)
cor_mat <- as.matrix(cor_wide[, -(1:3)])

ha <- rowAnnotation(
  Niche = as.character(cor_wide$niche),
  col = list(Niche = nicheCol[levels(cor_wide$niche)]),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Niche = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 5))
  )
)

png("figs/cosmx_correlation_heatmap.png", width = 6, height = 8, units = "in", res = 300)
ht <- Heatmap(
  cor_mat,
  name = "Spearman\ncorrelation",
  left_annotation = ha,
  row_split = cor_wide$niche,
  cluster_row_slices = FALSE,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)
  )
)
draw(ht)
dev.off()
cat("Saved: figs/cosmx_correlation_heatmap.png\n")


# =============================================================================
# 7. Summary statistics
# =============================================================================

cat("\n=== Mean pDC-fibroblast correlation by niche (CosMx) ===\n")
cor_results %>%
  group_by(niche, fibroblast) %>%
  summarise(
    mean_cor = mean(correlation),
    sd_cor = sd(correlation),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(niche) %>%
  print(n = 60)

cat("\nDone. All CosMx outputs saved to pDC/figs/ and pDC/output/\n")
