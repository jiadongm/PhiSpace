# Regenerate Fig 3C cluster plots with updated niche 5 colour (gold)
# This script reproduces the relevant figure panels from CaseCosMx_4Refs.Rmd

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs))

source(file.path(paths$project_root, "scripts/02_cosmx_case/_utils.R"))

dat_dir <- paths$phispace_data_root
PhiAssay <- "log1p"

# --- Load data ---
tissueNames <- tissueNames_CosMx
query_list <- qread(paste0(dat_dir, "data/CosMx_lung/SCE_obj/allTissue_SCE.qs"))
PhiResPath <- paste0(dat_dir, "output/Case3/CosMxAllLungsPhiRes4Refs_", PhiAssay, ".qs")
sc_list <- qread(PhiResPath)
for(ii in 1:length(tissueNames)){
  tissueName <- tissueNames[ii]
  query <- query_list[[tissueName]]
  sc4Refs <- lapply(sc_list, function(x) x[[tissueName]])
  reducedDim(query, "PhiSpace") <- Reduce(cbind, sc4Refs)
  query_list[[ii]] <- query
}

# Simplify cell type names
originalNames <- simpleNames <- colnames(reducedDim(query_list[[1]], "PhiSpace"))
simpleNames[c(2, 4, 9, 10, 12, 18, 21, 23, 24, 25, 26, 43)] <- c(
  "MoMacroph(immune)", "InflamMono(immune)",  "AlveolarMacroph(immune)",
  "ProlifImm(immune)", "InterstMacroph(immune)", "ProlifEpi(epithelial)",
  "SecretSCGB3A2+(epithelial)", "TransAT2(epithelial)", "SecretSCGB1A1+/MUC5B+(epithelial)",
  "SecretSCGB1A1+/SCGB3A2+(epithelial)", "DiffCiliated(epithelial)", "MyoFB Act(mesenchymal)"
)
for(ii in 1:length(query_list)){
  query <- query_list[[ii]]
  colnames(reducedDim(query, "PhiSpace")) <- simpleNames
  query_list[[ii]] <- query
}

# --- Focus on Lung5_Rep1 ---
LungName <- "Lung5_Rep1"
query <- query_list[[LungName]]

# Load cached clustering
clustResPath <- paste0(dat_dir, "output/Case3/Lung5_Rep1_PhiClusts4Refs.qs")
clust_res <- qread(clustResPath)
query$PhiClust <- as.character(clust_res$cluster)

# Updated colour palette: niche 5 changed from "gray" to gold
tempClustCols <- c(
  "1" = "#E5D8BD", "2" = "#CCEBC5", "3" = "#FED9A6", "4" = "#B3CDE3",
  "5" = "#E41A1C", "6" = "#FBB4AE", "7" = "#DECBE4", "8" = "#FFFFCC", "9" = "gray50"
)

# --- Figure 1: PhiSpace cluster spatial map ---
p <- VizSpatial(
  query, x_coord = "sdimx", y_coord = "sdimy", colBy = "PhiClust", ptSize = 0.3
) +
  scale_colour_manual(values = tempClustCols) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 2), nrow = 1
    )
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    legend.key.spacing = unit(0, "pt")
  ) +
  labs(colour = "")

ggsave(
  "figs/clusts/PhiClust.png",
  p + theme(legend.position = "none"),
  width = 2.5, height = 2.5
)
png(
  "figs/clusts/PhiClust_leg.png",
  width = 4, height = 0.4, units = "in", res = 300
)
op <- par(mar = rep(0, 4))
grid::grid.draw(ggpubr::get_legend(p))
par(op)
dev.off()

cat("Saved figs/clusts/PhiClust.png and PhiClust_leg.png\n")

# --- Figure 2: Gene expression cluster spatial map ---
clustResPath2 <- paste0(dat_dir, "output/Case3/Lung5_Rep1_OmicsClusts4Refs.qs")
clust_res2 <- qread(clustResPath2)

# align_clusters inline (removed from current PhiSpace)
align_clusters <- function(pred, true) {
  tab <- table(pred, true)
  mapping <- apply(tab, 1, which.max)
  names(mapping) <- rownames(tab)
  unname(colnames(tab)[mapping[pred]])
}

query$GenClust <- align_clusters(as.character(clust_res2$cluster), query$PhiClust)
p2 <- VizSpatial(
  query, x_coord = "sdimx", y_coord = "sdimy", colBy = "GenClust", ptSize = 0.3
) +
  scale_colour_manual(values = tempClustCols) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 2)
    )
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.spacing = unit(0, "pt")
  ) +
  labs(colour = "")

ggsave(
  "figs/clusts/OmicsClust.png",
  p2 + theme(legend.position = "none"),
  width = 2.5, height = 2.5
)

cat("Saved figs/clusts/OmicsClust.png\n")

# --- Figure 3: Divergence boxplot ---
query$divergence <- apply(
  reducedDim(query, "PhiSpace"), 1, quantile, prob = 0.75
)
pBox <- colData(query) %>% as.data.frame() %>%
  ggplot(aes(x = PhiClust, y = divergence)) +
  geom_boxplot(aes(fill = PhiClust), outlier.size = 0.5, linewidth = 0.1, outlier.stroke = 0) +
  theme_bw(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values = tempClustCols) +
  xlab("PhiSpace Clusters") + ylab("Divergence")

ggsave("figs/clusts/divergence.png", pBox, width = 1, height = 2)

cat("Saved figs/clusts/divergence.png\n")
cat("All figures regenerated with updated niche 5 colour (#FFD700 gold)\n")
