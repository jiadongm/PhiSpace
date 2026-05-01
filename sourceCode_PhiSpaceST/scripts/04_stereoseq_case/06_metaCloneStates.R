## Meta-Clone Cell Type Proportions & Cellular State Signatures
## Addressing Reviewer 2 Comment 5:
## 1. Quantitative cell type proportions per meta-clone (extending Fig 4D3)
## 2. Cellular state scoring via gene expression signatures (UCell + MSigDB Hallmark)

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
  library(UCell)
  library(msigdbr)
  library(grid)
})

source("utils.R")
dat_dir <- paths$dawson_root

# align_clusters: remap cluster labels to best match a reference clustering
align_clusters <- function(clust, clust_ref) {
  cont_table <- table(clust, clust_ref)
  mapping <- apply(cont_table, 1, which.max)
  new_labels <- colnames(cont_table)[mapping[as.character(clust)]]
  return(factor(new_labels))
}

# ============================================================================
# Step 0: Load Data
# ============================================================================

cat("=== Step 0: Loading data ===\n")

# Load annotated Stereo-seq query (from multiRefConsistency.R)
query <- qread("output/stereo_bin50_PhiSpace.qs")
cat(sprintf("  Stereo-seq query: %d bins x %d genes\n", ncol(query), nrow(query)))
cat(sprintf("  PhiSpace scores: %d bins x %d cell types\n",
            nrow(reducedDim(query, "PhiSpace")),
            ncol(reducedDim(query, "PhiSpace"))))

# Load barcode clustering results
outClusts <- readRDS(paste0(dat_dir, "output/barcode-seq_clustering.rds"))

# Reconstruct clust_list with alignment (same as Case_bridge_Mouse4.Rmd lines 260-284)
clust_list <- lapply(
  outClusts,
  function(x) {
    factor(x$cluster, levels = sort(unique(x$cluster)))
  }
)
for (x in length(clust_list):1) {
  if (x < length(clust_list)) {
    clust_old <- clust_list[[x + 1]]
    clust <- align_clusters(clust_list[[x]], clust_old) %>% as.factor()
    clust_list[[x]] <- factor(clust, levels = as.character(sort(as.numeric(levels(clust)))))
  }
}

# Load barcode-seq assay to get bin IDs (same as Case_bridge_Mouse4.Rmd lines 192-213)
bcRaw <- read.table(
  paste0(dat_dir, "stereo-seq_m4_additional_bin_sizes/mouse4_bin50_bc_counts.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
bcRaw$bin_y <- -bcRaw$bin_y
bcRaw$y_center <- -bcRaw$y_center
bcRaw$barcode <- gsub("GFP_", "", bcRaw$barcode)

sortedBC <- table(bcRaw$barcode) %>% sort(decreasing = TRUE)
sortedBC_filt <- sortedBC[sortedBC >= 50]
bcFilt <- bcRaw %>%
  filter(barcode %in% names(sortedBC_filt)) %>%
  mutate(barcode = as.factor(barcode))

# Spatially enhanced barcode-seq assay (xyRange = 500)
binIDs <- unique(bcFilt$cell_id)
barAssayPath <- paste0(dat_dir, "output/barAssay_xyRange500.rds")
if (file.exists(barAssayPath)) {
  barAssay <- readRDS(barAssayPath)
  cat("  Loaded cached barAssay\n")
} else {
  barAssay <- sapply(
    seq_along(binIDs),
    function(x) {
      binID <- binIDs[x]
      coords <- bcFilt %>% filter(cell_id == binID)
      coords <- coords[1, c("x_center", "y_center")]
      xcond <- (abs(bcFilt$x_center - coords[1, 1]) < 500)
      ycond <- (abs(bcFilt$y_center - coords[1, 2]) < 500)
      idx <- xcond & ycond
      table(bcFilt[idx, "barcode"])
    }
  ) %>% t() %>% `rownames<-`(binIDs)
  saveRDS(barAssay, barAssayPath)
  cat("  Computed and cached barAssay\n")
}

# ============================================================================
# Step 1: Define Meta-Clone Assignments
# ============================================================================

cat("\n=== Step 1: Defining meta-clone assignments ===\n")

# clustIdx = 7 corresponds to k=9 clustering (index into list of k=2..10)
clustIdx <- 7
clustBar <- clust_list[[clustIdx]] %>%
  `names<-`(rownames(barAssay))

# PhiSpace scores
phiScores <- reducedDim(query, "PhiSpace") %>% as.data.frame()

# Assign meta-clone labels
clustBar <- clustBar[intersect(rownames(phiScores), names(clustBar))]
clust <- rep("background", nrow(phiScores))
names(clust) <- rownames(phiScores)
clust[names(clustBar)] <- as.character(clustBar)

# Select meta-clones (excluding 1, 4, 9 as in original analysis)
selectedClust <- c("2", "3", "5", "6", "7", "8", "10", "background")
clust <- clust[clust %in% selectedClust]
clust <- factor(clust, levels = selectedClust)

# Build plot data
plot_dat <- phiScores[names(clust), ] %>%
  mutate(
    cluster = clust,
    x = query[, names(clust)]$x,
    y = query[, names(clust)]$y
  )

cTypes <- colnames(reducedDim(query, "PhiSpace"))
cat(sprintf("  %d bins assigned to meta-clones/background\n", length(clust)))
cat("  Meta-clone sizes:\n")
print(table(clust))

# Reference grouping for cell types
refSuffixes <- sub(".*\\((.*)\\)", "\\1", cTypes)

# ============================================================================
# Step 2: Cell Type Proportion Analysis
# ============================================================================

cat("\n=== Step 2: Cell type proportion analysis ===\n")

# --- 2A: Mean PhiSpace score heatmap per meta-clone ---
cat("  Computing mean PhiSpace scores per meta-clone...\n")

meanScores <- plot_dat %>%
  group_by(cluster) %>%
  summarise(across(all_of(cTypes), mean)) %>%
  as.data.frame()
rownames(meanScores) <- meanScores$cluster
meanScores$cluster <- NULL
meanScores <- t(as.matrix(meanScores))  # cell types x meta-clones

# Save quantitative table
write.csv(meanScores, "figs/metaClone/mean_scores_table.csv")
cat("  Saved: figs/metaClone/mean_scores_table.csv\n")

# Heatmap annotations
refGroup <- sub(".*\\((.*)\\)", "\\1", rownames(meanScores))
refCols <- c(
  Spleen = "#E41A1C", BM = "#377EB8",
  Neutro = "#4DAF4A", CITE = "#984EA3"
)
nBins <- table(clust)[colnames(meanScores)]

# Row annotation: reference
ra <- rowAnnotation(
  Reference = refGroup,
  col = list(Reference = refCols),
  show_legend = TRUE,
  simple_anno_size = unit(3, "mm")
)

# Column annotation: number of bins
ha <- HeatmapAnnotation(
  `N bins` = anno_text(
    as.character(nBins),
    gp = gpar(fontsize = 7),
    just = "center"
  ),
  simple_anno_size = unit(4, "mm")
)

# Simplify row labels
shortNames <- sub("\\(.*\\)", "", rownames(meanScores))

# Scale scores by row (z-score per cell type) for better visualisation
meanScores_scaled <- t(scale(t(meanScores)))

png("figs/metaClone/celltype_heatmap.png",
    width = 8, height = 14, units = "in", res = 300)
ht <- Heatmap(
  meanScores_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  left_annotation = ra,
  top_annotation = ha,
  row_labels = shortNames,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 0,
  row_names_max_width = unit(5, "cm"),
  row_split = factor(refGroup, levels = c("Spleen", "BM", "Neutro", "CITE")),
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  column_title = "Mean PhiSpace score per meta-clone (z-scaled)"
)
draw(ht, merge_legend = TRUE)
dev.off()
cat("  Saved: figs/metaClone/celltype_heatmap.png\n")

# --- 2B: Stacked bar plot of dominant cell type proportions ---
cat("  Computing dominant cell type proportions...\n")

# For each bin, assign dominant cell type = highest PhiSpace score
scoresMat <- reducedDim(query, "PhiSpace")[names(clust), ]
dominantType <- cTypes[apply(scoresMat, 1, which.max)]
names(dominantType) <- names(clust)

# Compute proportions per meta-clone
propDat <- data.frame(
  cluster = clust,
  dominant = dominantType,
  stringsAsFactors = FALSE
) %>%
  group_by(cluster, dominant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Keep top cell types for readability (aggregate rare ones as "Other")
topTypes <- propDat %>%
  group_by(dominant) %>%
  summarise(total = sum(n)) %>%
  arrange(-total) %>%
  slice_head(n = 15) %>%
  pull(dominant)

propDat <- propDat %>%
  mutate(dominant_label = ifelse(dominant %in% topTypes, dominant, "Other")) %>%
  group_by(cluster, dominant_label) %>%
  summarise(prop = sum(prop), .groups = "drop")

# Colour palette for dominant types
nCols <- length(unique(propDat$dominant_label))
domCols <- c(
  scales::hue_pal()(min(nCols - 1, 15)),
  "grey50"  # for "Other"
)
names(domCols) <- c(
  sort(setdiff(unique(propDat$dominant_label), "Other")),
  "Other"
)

p_bar <- ggplot(propDat, aes(x = cluster, y = prop, fill = dominant_label)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = domCols, name = "Dominant\ncell type") +
  theme_pubr(base_size = 10) +
  xlab("Meta-clone") +
  ylab("Proportion of bins") +
  ggtitle("Dominant cell type composition per meta-clone") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm")
  )

ggsave("figs/metaClone/dominant_celltype_barplot.png",
       p_bar, width = 10, height = 5, dpi = 300)
cat("  Saved: figs/metaClone/dominant_celltype_barplot.png\n")


# ============================================================================
# Step 3: UCell Scoring of Hallmark Gene Sets
# ============================================================================

cat("\n=== Step 3: UCell scoring of Hallmark gene sets ===\n")

# Get mouse Hallmark gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  select(gs_name, gene_symbol)

# Selected biologically relevant Hallmark sets
selected_hallmarks <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_COMPLEMENT"
)

# Filter to selected sets and convert to list format for UCell
sig_list <- hallmark_sets %>%
  filter(gs_name %in% selected_hallmarks) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# Check gene overlap
query_genes <- rownames(query)
for (nm in names(sig_list)) {
  overlap <- sum(sig_list[[nm]] %in% query_genes)
  total <- length(sig_list[[nm]])
  cat(sprintf("  %s: %d/%d genes (%.0f%%)\n", nm, overlap, total, 100 * overlap / total))
}

# Run UCell scoring
cat("  Running UCell scoring...\n")
ucellPath <- "output/ucell_hallmark_scores.rds"
if (file.exists(ucellPath)) {
  ucell_scores <- readRDS(ucellPath)
  cat("  Loaded cached UCell scores\n")
} else {
  # UCell expects a matrix with genes as rows and cells as columns
  count_mat <- assay(query, "counts")
  ucell_scores <- ScoreSignatures_UCell(
    count_mat,
    features = sig_list,
    name = ""  # no suffix
  )
  saveRDS(ucell_scores, ucellPath)
  cat("  Computed and cached UCell scores\n")
}

cat(sprintf("  UCell scores: %d bins x %d signatures\n",
            nrow(ucell_scores), ncol(ucell_scores)))

# Clean up column names (remove trailing underscore if present)
colnames(ucell_scores) <- gsub("_$", "", colnames(ucell_scores))

# Short display names for signatures
sig_short_names <- c(
  "HALLMARK_G2M_CHECKPOINT" = "G2M checkpoint",
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "HALLMARK_MYC_TARGETS_V1" = "MYC targets",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFa/NF-kB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "IFN-alpha",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFN-gamma",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6/JAK/STAT3",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "OxPhos",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_COMPLEMENT" = "Complement"
)

# ============================================================================
# Step 4: Spatial Heatmaps of Signature Scores
# ============================================================================

cat("\n=== Step 4: Spatial heatmaps of signature scores ===\n")

coords <- colData(query)[, c("x", "y")] %>% as.data.frame()

# Spatial heatmap helper (from multiRefConsistency.R)
spatialHeatmap <- function(scores_vec, coords, title, ptSize = 0.3, titleSize = 8) {
  vals <- PhiSpace:::censor(scores_vec, quant = 0.95)
  pdat <- coords %>%
    mutate(score = vals) %>%
    arrange(score)

  p <- ggplot(pdat, aes(x, y)) +
    geom_point(aes(colour = score), size = ptSize, stroke = 0) +
    scale_colour_gradientn(colours = MATLAB_cols) +
    theme_void(base_size = titleSize) +
    ggtitle(title) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = titleSize, face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  return(p)
}

# Generate spatial heatmaps for all 11 signatures
sig_names <- intersect(selected_hallmarks, colnames(ucell_scores))
sig_plots <- lapply(sig_names, function(sig) {
  short_name <- sig_short_names[sig]
  spatialHeatmap(ucell_scores[colnames(query), sig], coords, short_name, ptSize = 0.5, titleSize = 7)
})

nSigCols <- 4
nSigRows <- ceiling(length(sig_plots) / nSigCols)

p_sig_spatial <- ggarrange(
  plotlist = sig_plots,
  ncol = nSigCols, nrow = nSigRows
)

ggsave("figs/metaClone/signature_spatial_maps.png",
       p_sig_spatial, width = 2.5 * nSigCols, height = 2.2 * nSigRows, dpi = 300)
cat("  Saved: figs/metaClone/signature_spatial_maps.png\n")

# ============================================================================
# Step 5: Meta-Clone Signature Enrichment Heatmap
# ============================================================================

cat("\n=== Step 5: Meta-clone signature enrichment ===\n")

# Build signature data for bins with meta-clone assignments
ucell_clust <- ucell_scores[names(clust), sig_names, drop = FALSE]

# Mean signature score per meta-clone
meanSig <- data.frame(
  cluster = clust,
  ucell_clust,
  check.names = FALSE
) %>%
  group_by(cluster) %>%
  summarise(across(all_of(sig_names), mean)) %>%
  as.data.frame()
rownames(meanSig) <- meanSig$cluster
meanSig$cluster <- NULL
meanSig <- t(as.matrix(meanSig))  # signatures x meta-clones

# Short names for rows
rownames(meanSig) <- sig_short_names[rownames(meanSig)]

# Z-scale per signature for visualisation
meanSig_scaled <- t(scale(t(meanSig)))

# Wilcoxon test + fold change enrichment analysis
sigEnrich_pval <- sapply(sig_names, function(sig) {
  sc_split <- split(ucell_clust[, sig], clust)
  bkgrd <- sc_split$background
  sc_split[["background"]] <- NULL
  sapply(sc_split, function(x) wilcox.test(x, bkgrd)$p.value)
}) %>% t()

sigEnrich_fc <- sapply(sig_names, function(sig) {
  sc_split <- split(ucell_clust[, sig], clust)
  bkgrd <- sc_split$background
  sc_split[["background"]] <- NULL
  sapply(sc_split, function(x) (mean(x) - mean(bkgrd)) / sd(bkgrd))
}) %>% t()

sigScores_sig <- sigEnrich_fc * (-log10(sigEnrich_pval))
rownames(sigScores_sig) <- sig_short_names[sig_names]

# Heatmap of significance scores (fc * -log10(pval))
png("figs/metaClone/signature_metaClone_heatmap.png",
    width = 7, height = 5, units = "in", res = 300)
ht_sig <- Heatmap(
  sigScores_sig,
  name = "Enrichment\nscore",
  col = colorRamp2(
    c(min(sigScores_sig), 0, max(sigScores_sig)),
    c("blue", "white", "red")
  ),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  cluster_columns = FALSE,
  border = TRUE,
  show_column_dend = FALSE,
  column_title = "Signature enrichment per meta-clone\n(fold change x -log10 p-value vs background)"
)
draw(ht_sig)
dev.off()
cat("  Saved: figs/metaClone/signature_metaClone_heatmap.png\n")

# Also save the z-scaled mean signature heatmap
png("figs/metaClone/signature_metaClone_meanZ.png",
    width = 7, height = 5, units = "in", res = 300)
ht_meanZ <- Heatmap(
  meanSig_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  cluster_columns = FALSE,
  border = TRUE,
  show_column_dend = FALSE,
  column_title = "Mean signature score per meta-clone (z-scaled)"
)
draw(ht_meanZ)
dev.off()
cat("  Saved: figs/metaClone/signature_metaClone_meanZ.png\n")

# ============================================================================
# Step 6: Cell Type vs Signature Correlation
# ============================================================================

cat("\n=== Step 6: Cell type vs signature correlation ===\n")

# Select biologically interesting cell types for correlation
selected_ctypes <- c(
  "T1(Neutro)", "MAT 1(Neutro)", "IMM 1(Neutro)", "PreNeu(Neutro)",
  "ICOS+ Tregs(CITE)", "CD4 T(CITE)", "CD8 T(CITE)", "NK(CITE)",
  "HPC(BM)", "ErythBla(BM)", "ProEryThBla(BM)",
  "RedPulp macro(CITE)", "Neutro(Spleen)", "Macro(Spleen)"
)
# Keep only those that exist
selected_ctypes <- selected_ctypes[selected_ctypes %in% cTypes]

# Compute Spearman correlations: selected cell types x signatures
ct_scores <- reducedDim(query, "PhiSpace")[, selected_ctypes]
sig_scores_mat <- ucell_scores[colnames(query), sig_names]

cor_mat <- cor(ct_scores, sig_scores_mat, method = "spearman")
colnames(cor_mat) <- sig_short_names[colnames(cor_mat)]
# Simplify row names
rownames(cor_mat) <- sub("\\(.*\\)", "", rownames(cor_mat))

png("figs/metaClone/signature_celltype_correlation.png",
    width = 8, height = 6, units = "in", res = 300)
ht_cor <- Heatmap(
  cor_mat,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  border = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title = "Spearman correlation: cell type scores vs signature scores"
)
draw(ht_cor)
dev.off()
cat("  Saved: figs/metaClone/signature_celltype_correlation.png\n")


cat("\n=== All done! ===\n")
cat("Output files:\n")
cat("  figs/metaClone/celltype_heatmap.png\n")
cat("  figs/metaClone/dominant_celltype_barplot.png\n")
cat("  figs/metaClone/mean_scores_table.csv\n")
cat("  figs/metaClone/signature_spatial_maps.png\n")
cat("  figs/metaClone/signature_metaClone_heatmap.png\n")
cat("  figs/metaClone/signature_metaClone_meanZ.png\n")
cat("  figs/metaClone/signature_celltype_correlation.png\n")
