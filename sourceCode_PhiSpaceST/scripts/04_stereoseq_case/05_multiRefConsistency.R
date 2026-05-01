## Multi-Reference Consistency Analysis
## Addressing Reviewer 2 Comment 4:
## How are common/unique cell types across references handled?
## How consistent are the results for the same cell types?
##
## Uses the Stereo-seq AML mouse spleen case study with 4 references:
## Spleen (18 types), BM (13 types), Neutro (11 types), CITE (28 types)

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(PhiSpace)
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(ggpubr)
  library(tidyr)
  library(qs)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("utils.R")
dat_dir <- paths$dawson_root

# ============================================================================
# Step 1: Load Data
# ============================================================================

cat("=== Step 1: Loading data ===\n")

# Per-reference PhiSpace scores for intermediate scRNA-seq
scPhiSc_list <- readRDS(
  paste0(dat_dir, "output/Mouse4_scRNA-seq_PhiSc_list.rds")
)

# Cell type name mapping (same as in Case_bridge_Mouse4.Rmd)
newNames <- c(
  "RBC", "CD8 T", "Plasma", "Mature B", "cDC",
  "NK", "Neutro", "CD4 T", "Cr2 B", "HSPC",
  "Trans B", "Mono", "T", "Pre-B cycl", "Endo",
  "Cd7+ NK", "Macro", "pDC", "Baso", "Pre-B",
  "Mono", "Naive B", "Granulo", "Macro", "Imm B",
  "HPC", "Late pro-B", "Imm NK", "ProEryThBla", "T",
  "ErythBla", "T3", "MAT 3", "MAT 2", "IMM 1",
  "MAT 4", "MAT 1", "IMM 2", "T1", "T2",
  "PreNeu", "MAT 5", "MZ B", "Trans B", "Mat B",
  "CD122+ B", "Ifit3+CD4 T", "CD4 T", "pDC", "CD8 T",
  "B1 B", "NKT", "Ifit3+ B", "NK", "GD T",
  "ICOS+ Tregs", "Tregs", "Ly6+ mono", "Neutro", "Cycl B/T",
  "cDC2", "Ly6- mono", "RedPulp macro", "RBC", "Mig DC",
  "Ifit3+CD8 T", "cDC1", "Act CD4 T", "Plasma", "MZ/Marco+ macro"
)
refSuffixes <- rep(
  c("(Spleen)", "(BM)", "(Neutro)", "(CITE)"),
  sapply(scPhiSc_list, ncol)
)
fullNames <- paste0(newNames, refSuffixes)

# Assign readable names to each per-reference score matrix
refNames <- c("Spleen", "BM", "Neutro", "CITE")
offset <- 0
for (i in seq_along(scPhiSc_list)) {
  nc <- ncol(scPhiSc_list[[i]])
  colnames(scPhiSc_list[[i]]) <- fullNames[(offset + 1):(offset + nc)]
  offset <- offset + nc
}

cat("Reference dimensions:\n")
for (nm in names(scPhiSc_list)) {
  cat(sprintf("  %s: %d cells x %d cell types\n",
              nm, nrow(scPhiSc_list[[nm]]), ncol(scPhiSc_list[[nm]])))
}

# ============================================================================
# Step 2: Cross-Reference Correlation of Overlapping Cell Types
# ============================================================================

cat("\n=== Step 2: Cross-reference correlation analysis ===\n")

# Combine all per-reference scores into one matrix
allScores <- do.call("cbind", scPhiSc_list)
cat(sprintf("Combined score matrix: %d cells x %d cell types\n",
            nrow(allScores), ncol(allScores)))

# Define overlapping cell type pairs
# Format: list of groups, each with (broad_name, list of (colname, ref) pairs)
overlap_groups <- list(
  list(name = "CD4 T",
       members = c("CD4 T(Spleen)", "CD4 T(CITE)")),
  list(name = "CD8 T",
       members = c("CD8 T(Spleen)", "CD8 T(CITE)")),
  list(name = "T (general)",
       members = c("T(Spleen)", "T(BM)")),
  list(name = "Mature B",
       members = c("Mature B(Spleen)", "Mat B(CITE)")),
  list(name = "Transitional B",
       members = c("Trans B(Spleen)", "Trans B(CITE)")),
  list(name = "NK",
       members = c("NK(Spleen)", "Imm NK(BM)", "NK(CITE)")),
  list(name = "Macrophage",
       members = c("Macro(Spleen)", "Macro(BM)",
                    "RedPulp macro(CITE)", "MZ/Marco+ macro(CITE)")),
  list(name = "Monocyte",
       members = c("Mono(Spleen)", "Mono(BM)",
                    "Ly6+ mono(CITE)", "Ly6- mono(CITE)")),
  list(name = "pDC",
       members = c("pDC(Spleen)", "pDC(CITE)")),
  list(name = "cDC",
       members = c("cDC(Spleen)", "cDC1(CITE)", "cDC2(CITE)")),
  list(name = "Neutrophil",
       members = c("Neutro(Spleen)", "Neutro(CITE)")),
  list(name = "Erythrocyte",
       members = c("RBC(Spleen)", "ErythBla(BM)", "RBC(CITE)")),
  list(name = "Plasma",
       members = c("Plasma(Spleen)", "Plasma(CITE)"))
)

# Compute pairwise Spearman correlations for all overlapping pairs
pair_results <- data.frame(
  group = character(),
  celltype1 = character(),
  celltype2 = character(),
  rho = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (grp in overlap_groups) {
  members <- grp$members
  # Verify all members exist in the score matrix
  members <- members[members %in% colnames(allScores)]
  if (length(members) < 2) next

  for (i in 1:(length(members) - 1)) {
    for (j in (i + 1):length(members)) {
      ct <- cor.test(
        allScores[, members[i]],
        allScores[, members[j]],
        method = "spearman"
      )
      pair_results <- rbind(pair_results, data.frame(
        group = grp$name,
        celltype1 = members[i],
        celltype2 = members[j],
        rho = ct$estimate,
        pvalue = ct$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}
rownames(pair_results) <- NULL

cat("\nOverlapping cell type correlations:\n")
print(pair_results[order(-pair_results$rho), ], row.names = FALSE)

# Save the table
write.csv(pair_results, "figs/multiRef/crossRef_overlap_correlations.csv",
          row.names = FALSE)

# --- Figure 2a: Full cross-reference correlation heatmap ---
cat("\nGenerating full 70x70 correlation heatmap...\n")

corMat <- cor(allScores, method = "spearman")

# Reference grouping annotation
refGroup <- sub(".*\\((.*)\\)", "\\1", colnames(allScores))
refCols <- c(
  Spleen = "#E41A1C", BM = "#377EB8",
  Neutro = "#4DAF4A", CITE = "#984EA3"
)
ha <- HeatmapAnnotation(
  Reference = refGroup,
  col = list(Reference = refCols),
  show_legend = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)
ra <- rowAnnotation(
  Reference = refGroup,
  col = list(Reference = refCols),
  show_legend = FALSE,
  simple_anno_size = unit(3, "mm")
)

# Simplify labels: remove reference suffix for readability
shortNames <- sub("\\(.*\\)", "", colnames(allScores))

png("figs/multiRef/crossRef_correlation_heatmap.png",
    width = 12, height = 10, units = "in", res = 300)
ht <- Heatmap(
  corMat,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.5, 0, 0.5, 1), c("blue", "white", "orange", "red")),
  top_annotation = ha,
  left_annotation = ra,
  row_labels = shortNames,
  column_labels = shortNames,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_names_max_width = unit(4, "cm"),
  column_names_max_height = unit(4, "cm"),
  column_split = factor(refGroup, levels = c("Spleen", "BM", "Neutro", "CITE")),
  row_split = factor(refGroup, levels = c("Spleen", "BM", "Neutro", "CITE")),
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  border = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)
draw(ht, merge_legend = TRUE)
dev.off()
cat("  Saved: figs/multiRef/crossRef_correlation_heatmap.png\n")

# --- Figure 2b: Bar plot of overlap pair correlations ---
cat("Generating overlap pair bar plot...\n")

pair_results$pair_label <- paste0(
  pair_results$celltype1, " vs\n", pair_results$celltype2
)
pair_results <- pair_results[order(pair_results$rho, decreasing = TRUE), ]
pair_results$pair_label <- factor(
  pair_results$pair_label,
  levels = rev(pair_results$pair_label)
)

p_bar <- ggplot(pair_results, aes(x = rho, y = pair_label)) +
  geom_bar(
    aes(fill = rho),
    stat = "identity"
  ) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  theme_pubr(base_size = 8) +
  xlab("Spearman rho") +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("figs/multiRef/crossRef_overlap_barplot.png",
       p_bar, width = 5, height = 6, dpi = 300)
cat("  Saved: figs/multiRef/crossRef_overlap_barplot.png\n")


# ============================================================================
# Step 3: Bridge PhiSpace scores to Stereo-seq
# ============================================================================

cat("\n=== Step 3: Bridging PhiSpace scores to Stereo-seq ===\n")

outPath <- "output/stereo_bin50_PhiSpace.qs"

if (!file.exists(outPath)) {
  # Load intermediate scRNA-seq with concatenated PhiSpace scores
  querySC <- qread(paste0(dat_dir, "mouse4_scRNAseq_sce.qs"))
  reducedDim(querySC, "PhiSpace") <- allScores
  querySC <- logTransf(querySC, targetAssay = "log1p", use_log1p = TRUE)

  # Load Stereo-seq query
  query <- qread(
    paste0(dat_dir,
           "stereo-seq_m4_additional_bin_sizes/mouse4_bin50_bc_sce.qs")
  )
  # Filter low-count bins
  cell2keep <- colnames(query)[
    query$total_counts > quantile(query$total_counts, 0.01)
  ]
  query <- query[, cell2keep]

  # Bridge: use scRNA-seq PhiSpace scores as response
  PhiRes <- PhiSpaceR_1ref(
    querySC,
    query,
    response = reducedDim(querySC, "PhiSpace"),
    refAssay = "log1p",
    nfeat = 500,
    regMethod = "PLS",
    scale = FALSE
  )
  reducedDim(query, "PhiSpace") <- normPhiScores(PhiRes$PhiSpaceScore)

  qsave(query, outPath)
  cat("  Saved annotated Stereo-seq:", outPath, "\n")
} else {
  query <- qread(outPath)
  cat("  Loaded existing annotated Stereo-seq:", outPath, "\n")
}

cat(sprintf("  Stereo-seq: %d bins x %d cell types\n",
            ncol(query), ncol(reducedDim(query, "PhiSpace"))))


# ============================================================================
# Step 4: Spatial Maps for Overlapping Cell Types
# ============================================================================

cat("\n=== Step 4: Spatial maps for overlapping cell types ===\n")

scores_st <- reducedDim(query, "PhiSpace")
coords <- colData(query)[, c("x", "y")] %>% as.data.frame()
coords$cell <- rownames(coords)

# Helper: spatial heatmap for a single cell type
spatialHeatmap <- function(scores, coords, ctype, colour_limits = NULL,
                           ptSize = 0.3, titleSize = 8) {
  vals <- scores[, ctype]
  vals <- PhiSpace:::censor(vals, quant = 0.95)
  pdat <- coords %>%
    mutate(score = vals) %>%
    arrange(score)

  p <- ggplot(pdat, aes(x, y)) +
    geom_point(aes(colour = score), size = ptSize, stroke = 0) +
    scale_colour_gradientn(
      colours = MATLAB_cols,
      limits = colour_limits
    ) +
    theme_void(base_size = titleSize) +
    ggtitle(ctype) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = titleSize, face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  return(p)
}

# Selected overlap groups for spatial comparison
spatial_overlap_groups <- list(
  list(name = "CD4 T",
       members = c("CD4 T(Spleen)", "CD4 T(CITE)")),
  list(name = "Mature B",
       members = c("Mature B(Spleen)", "Mat B(CITE)")),
  list(name = "NK",
       members = c("NK(Spleen)", "NK(CITE)")),
  list(name = "Macrophage",
       members = c("Macro(Spleen)", "Macro(BM)", "RedPulp macro(CITE)")),
  list(name = "Erythrocyte",
       members = c("RBC(Spleen)", "ErythBla(BM)", "RBC(CITE)"))
)

# Determine max columns needed
maxCols <- max(sapply(spatial_overlap_groups, function(g) length(g$members)))

# Generate plots
all_row_plots <- list()
for (grp in spatial_overlap_groups) {
  members <- grp$members[grp$members %in% colnames(scores_st)]

  # Use shared colour limits within each group
  combined_vals <- unlist(lapply(members, function(m) {
    PhiSpace:::censor(scores_st[, m], quant = 0.95)
  }))
  clims <- range(combined_vals)

  row_plots <- lapply(members, function(m) {
    spatialHeatmap(scores_st, coords, m, colour_limits = clims)
  })

  # Pad with blank plots if fewer than maxCols
  while (length(row_plots) < maxCols) {
    row_plots <- c(row_plots, list(
      ggplot() + theme_void() +
        theme(plot.background = element_rect(fill = "white", colour = NA))
    ))
  }
  all_row_plots <- c(all_row_plots, row_plots)
}

# Add a colour bar legend from one example plot
legend_plot <- spatialHeatmap(
  scores_st, coords, "CD4 T(Spleen)", ptSize = 0.3
) + theme(legend.position = "right")

p_spatial_overlap <- ggarrange(
  plotlist = all_row_plots,
  ncol = maxCols,
  nrow = length(spatial_overlap_groups)
)

ggsave("figs/multiRef/spatial_overlap_comparison.png",
       p_spatial_overlap, width = 2.5 * maxCols, height = 2.2 * length(spatial_overlap_groups),
       dpi = 300)
cat("  Saved: figs/multiRef/spatial_overlap_comparison.png\n")

# Also compute Stereo-seq level correlations for the overlap pairs
cat("\nStereo-seq spatial score correlations for overlap pairs:\n")
st_pair_results <- data.frame(
  group = character(),
  celltype1 = character(),
  celltype2 = character(),
  rho = numeric(),
  stringsAsFactors = FALSE
)
for (grp in spatial_overlap_groups) {
  members <- grp$members[grp$members %in% colnames(scores_st)]
  if (length(members) < 2) next
  for (i in 1:(length(members) - 1)) {
    for (j in (i + 1):length(members)) {
      rho <- cor(scores_st[, members[i]], scores_st[, members[j]],
                 method = "spearman")
      st_pair_results <- rbind(st_pair_results, data.frame(
        group = grp$name,
        celltype1 = members[i],
        celltype2 = members[j],
        rho = rho,
        stringsAsFactors = FALSE
      ))
    }
  }
}
print(st_pair_results[order(-st_pair_results$rho), ], row.names = FALSE)
write.csv(st_pair_results,
          "figs/multiRef/stereo_overlap_correlations.csv",
          row.names = FALSE)

# ============================================================================
# Step 5: Unique Cell Type Spatial Maps
# ============================================================================

cat("\n=== Step 5: Spatial maps for reference-unique cell types ===\n")

unique_types <- list(
  "Neutro-specific" = c("MAT 3(Neutro)", "MAT 1(Neutro)", "T1(Neutro)", "PreNeu(Neutro)"),
  "CITE-specific"   = c("ICOS+ Tregs(CITE)", "GD T(CITE)", "NKT(CITE)", "Mig DC(CITE)"),
  "BM-specific"     = c("HPC(BM)", "ProEryThBla(BM)", "Late pro-B(BM)", "Naive B(BM)")
)

unique_plots <- list()
for (refLabel in names(unique_types)) {
  ctypes <- unique_types[[refLabel]]
  ctypes <- ctypes[ctypes %in% colnames(scores_st)]
  for (ct in ctypes) {
    unique_plots <- c(unique_plots, list(
      spatialHeatmap(scores_st, coords, ct)
    ))
  }
}

nUniqueCols <- 4
nUniqueRows <- ceiling(length(unique_plots) / nUniqueCols)

p_unique <- ggarrange(
  plotlist = unique_plots,
  ncol = nUniqueCols,
  nrow = nUniqueRows
)

ggsave("figs/multiRef/spatial_unique_celltypes.png",
       p_unique,
       width = 2.5 * nUniqueCols, height = 2.2 * nUniqueRows,
       dpi = 300)
cat("  Saved: figs/multiRef/spatial_unique_celltypes.png\n")


cat("\n=== All done! ===\n")
