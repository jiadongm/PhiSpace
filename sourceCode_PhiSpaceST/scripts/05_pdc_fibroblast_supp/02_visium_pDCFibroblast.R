## ============================================================================
## pDC-Fibroblast co-presence analysis in Visium NSCLC (all 40 samples)
## Addressing Reviewer 2, Comment 3
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

source("../CaseVisium/utils.R")
pdc_dir <- file.path(paths$output_root, "pDC")

# =============================================================================
# 1. Load data (all 40 samples)
# =============================================================================

# Sample metadata
meta <- read.csv(file.path(pdc_dir, "data", "sample_metadata.csv"))
sampleNames <- meta$sample
cancerTypes <- meta$condition
names(cancerTypes) <- sampleNames

# PhiSpace results
PhiRes <- qread(file.path(pdc_dir, "output", "PhiRes_allSamples.qs"))

# Query SCE list (for spatial coordinates)
query_list <- qread(file.path(pdc_dir, "data", "allSamples.qs"))

# Identify pDC and fibroblast columns
cTypeNames <- colnames(normPhiScores(PhiRes$PhiSpaceScore[[sampleNames[1]]]))
pDC_col <- grep("Plasmacytoid", cTypeNames, value = TRUE)
FB_cols <- grep("[Ff]ibroblast", cTypeNames, value = TRUE)
cat("pDC column:", pDC_col, "\n")
cat("Fibroblast columns:", paste(FB_cols, collapse = ", "), "\n")
cat("Total samples:", length(sampleNames), "\n")

# Short labels for fibroblast subtypes (for plotting)
FB_shortNames <- gsub(" fibroblasts", "", FB_cols)
names(FB_shortNames) <- FB_cols


# =============================================================================
# 2. Spatial heatmaps: representative samples
# =============================================================================

repSamples <- c(
  Healthy = "D2_1",
  LUAD    = "P16_T1",
  LUSC    = "P11_T3"
)

# Helper: spatial heatmap for a single cell type score in one sample
spatialHeatmap <- function(sce, scores, col_name, short_title = col_name,
                           coordNames = c("x", "y"), ptSize = 0.8) {
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

# Build multi-panel figure: rows = conditions, cols = pDC + fibroblasts + co-localization
panel_list <- list()
for (cond in names(repSamples)) {
  samp <- repSamples[cond]
  sce <- query_list[[samp]]
  scores <- normPhiScores(PhiRes$PhiSpaceScore[[samp]])

  # pDC
  panel_list[[paste0(cond, "_pDC")]] <- spatialHeatmap(
    sce, scores, pDC_col, short_title = paste0("pDC (", cond, ")")
  )

  # Each fibroblast subtype
  for (fb in FB_cols) {
    panel_list[[paste0(cond, "_", fb)]] <- spatialHeatmap(
      sce, scores, fb, short_title = FB_shortNames[fb]
    )
  }

  # Co-localization: pDC score * max fibroblast score per spot
  pDC_scores <- scores[, pDC_col]
  FB_scores <- scores[, FB_cols, drop = FALSE]
  max_FB <- apply(FB_scores, 1, max)
  coloc <- pDC_scores * max_FB
  coloc_mat <- matrix(coloc, ncol = 1, dimnames = list(names(coloc), "colocalization"))
  panel_list[[paste0(cond, "_coloc")]] <- spatialHeatmap(
    sce, coloc_mat, "colocalization", short_title = "pDC-FB co-loc"
  )
}

ncols <- 1 + length(FB_cols) + 1
p_spatial <- ggarrange(
  plotlist = panel_list,
  nrow = length(repSamples), ncol = ncols
)
ggsave(
  "figs/visium_spatial_heatmaps.png", p_spatial,
  width = 2 * ncols, height = 2 * length(repSamples), dpi = 300
)
cat("Saved: figs/visium_spatial_heatmaps.png\n")


# =============================================================================
# 3. Per-sample pDC-fibroblast Spearman correlations (all 40 samples)
# =============================================================================

cor_results <- data.frame(
  sample = character(),
  condition = character(),
  patient = character(),
  fibroblast = character(),
  correlation = numeric(),
  stringsAsFactors = FALSE
)

for (ii in seq_along(sampleNames)) {
  samp <- sampleNames[ii]
  cond <- cancerTypes[samp]
  patient <- sub("_.*", "", samp)
  scores <- normPhiScores(PhiRes$PhiSpaceScore[[samp]])
  pDC_scores <- scores[, pDC_col]

  for (fb in FB_cols) {
    fb_scores <- scores[, fb]
    rho <- cor(pDC_scores, fb_scores, method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      sample = samp,
      condition = cond,
      patient = patient,
      fibroblast = FB_shortNames[fb],
      correlation = rho,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(cor_results, "output/visium_pDC_fb_correlations.csv", row.names = FALSE)
cat("Saved: output/visium_pDC_fb_correlations.csv\n")


# =============================================================================
# 4. Boxplot: pDC-fibroblast correlation by condition (Healthy vs LUAD vs LUSC)
# =============================================================================

cor_tumor <- cor_results %>%
  filter(condition %in% c("Healthy", "LUAD", "LUSC")) %>%
  mutate(condition = factor(condition, levels = c("Healthy", "LUAD", "LUSC")))

p_box <- ggplot(cor_tumor, aes(x = condition, y = correlation, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  facet_wrap(~ fibroblast, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = VisiumCancerCols) +
  theme_pubr(base_size = 7) +
  ylab("Spearman correlation with pDC") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 6)
  )
ggsave(
  "figs/visium_correlation_boxplot.png", p_box,
  width = 1.5 * length(FB_cols), height = 3, dpi = 300
)
cat("Saved: figs/visium_correlation_boxplot.png\n")

# Including Background samples
cor_all <- cor_results %>%
  mutate(condition = factor(condition, levels = c("Healthy", "Background", "LUAD", "LUSC")))

p_box_all <- ggplot(cor_all, aes(x = condition, y = correlation, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  facet_wrap(~ fibroblast, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = VisiumCancerCols) +
  theme_pubr(base_size = 7) +
  ylab("Spearman correlation with pDC") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 6)
  )
ggsave(
  "figs/visium_correlation_boxplot_all.png", p_box_all,
  width = 1.5 * length(FB_cols), height = 3, dpi = 300
)
cat("Saved: figs/visium_correlation_boxplot_all.png\n")


# =============================================================================
# 5. Statistical tests
# =============================================================================

# Non-tumor (Healthy + Background) vs Tumor (LUAD + LUSC)
cor_binary <- cor_results %>%
  mutate(
    group = ifelse(condition %in% c("Healthy", "Background"), "Non-tumor", "Tumor"),
    group = factor(group, levels = c("Non-tumor", "Tumor"))
  )

cat("\n=== Statistical tests (Non-tumor vs Tumor) ===\n")
stat_results <- list()
for (fb in unique(cor_binary$fibroblast)) {
  sub <- cor_binary %>% filter(fibroblast == fb)
  wt <- wilcox.test(correlation ~ group, data = sub)

  stat_results[[fb]] <- list(wilcox_p = wt$p.value)
  cat(sprintf("  %s: Wilcoxon p = %.6f\n", fb, wt$p.value))
}


# =============================================================================
# 6. Correlation heatmap: samples x fibroblast subtypes
# =============================================================================

cor_wide <- cor_results %>%
  pivot_wider(names_from = fibroblast, values_from = correlation) %>%
  as.data.frame()
rownames(cor_wide) <- cor_wide$sample
cor_mat <- as.matrix(cor_wide[, -(1:3)])

# Row annotation: condition
ha <- rowAnnotation(
  Condition = cor_wide$condition,
  col = list(Condition = VisiumCancerCols),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Condition = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))
  )
)

png("figs/visium_correlation_heatmap.png", width = 5, height = 6, units = "in", res = 300)
ht <- Heatmap(
  cor_mat,
  name = "Spearman\ncorrelation",
  left_annotation = ha,
  row_split = factor(
    cor_wide$condition,
    levels = c("Healthy", "Background", "LUAD", "LUSC")
  ),
  cluster_row_slices = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)
  )
)
draw(ht)
dev.off()
cat("Saved: figs/visium_correlation_heatmap.png\n")


# =============================================================================
# 7. Summary statistics
# =============================================================================

cat("\n=== Mean pDC-fibroblast correlation by condition ===\n")
cor_results %>%
  group_by(condition, fibroblast) %>%
  summarise(
    mean_cor = mean(correlation),
    sd_cor = sd(correlation),
    n = n(),
    .groups = "drop"
  ) %>%
  print(n = 50)

cat("\nDone. All outputs saved to pDC/figs/ and pDC/output/\n")
