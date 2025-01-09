# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))


dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/PhiSpace-ST_submit/" # replace with your own directory

### scRNA-seq from healthy and fibrotic lungs (by lineages)
# Available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227136
seu <- readRDS(
  paste0(dat_dir, "data/LungFibrosis/scRNA-seq/GSE227136_ILD_immune_Seurat.rds")
)
seu <- UpdateSeuratObject(seu)
sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = seu[["RNA"]]@counts),
  colData = seu@meta.data
)
sce <- PhiSpace::subsample(sce, "manual_annotation_1", proportion = 0.1, minCellNum = 50)
sce <- PhiSpace::zeroFeatQC(sce)
sce <- PhiSpace::scranTransf(sce)
qsave(sce, paste0(dat_dir, "data/LungFibrosis/scRNA-seq/immune_sce.qs"))

seu <- readRDS(
  paste0(dat_dir, "data/LungFibrosis/scRNA-seq/GSE227136_ILD_epithelial_Seurat.rds")
)
seu <- UpdateSeuratObject(seu)
sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = seu[["RNA"]]@counts),
  colData = seu@meta.data
)
sce <- PhiSpace::subsample(sce, "manual_annotation_1", proportion = 0.1, minCellNum = 50)
sce <- PhiSpace::zeroFeatQC(sce)
sce <- PhiSpace::scranTransf(sce)
qsave(sce, paste0(dat_dir, "data/LungFibrosis/scRNA-seq/epithelial_sce.qs"))

seu <- readRDS(
  paste0(dat_dir, "data/LungFibrosis/scRNA-seq/GSE227136_ILD_endothelial_Seurat.rds")
)
seu <- UpdateSeuratObject(seu)
sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = seu[["RNA"]]@counts),
  colData = seu@meta.data
)
sce <- PhiSpace::subsample(sce, "manual_annotation_1", proportion = 0.1, minCellNum = 50)
sce <- PhiSpace::zeroFeatQC(sce)
sce <- PhiSpace::scranTransf(sce)
qsave(sce, paste0(dat_dir, "data/LungFibrosis/scRNA-seq/endothelial_sce.qs"))

seu <- SeuratDisk::LoadH5Seurat(
  paste0(dat_dir, "data/LungFibrosis/scRNA-seq/GSE227136_ILD_mesenchymal_Seurat.rds")
)
sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts = seu[["RNA"]]@counts),
  colData = seu@meta.data
)
sce <- PhiSpace::subsample(sce, "manual_annotation_1", proportion = 0.1, minCellNum = 50)
sce <- PhiSpace::zeroFeatQC(sce)
sce <- PhiSpace::scranTransf(sce)
qsave(sce, paste0(dat_dir, "data/LungFibrosis/scRNA-seq/mesenchymal_sce.qs"))




## CosMx: Convert original Giotto objects to SCE objects

# Tried: "Lung5_Rep2"
# "Lung13" is special as no cell type map correspond perfectly to cancer niche (low proportion of cancer cells)
# Lung12
# Lung6 Squamous cell

library(Giotto)
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/PhiSpace-ST_submit/" # replace with your own directory
load(paste0(dat_dir, "data/CosMx_lung/SMI_Giotto_Object.RData")) 


## Separate lung tissues
table(gem@cell_metadata$rna$Run_Tissue_name)
# ! note that Lung5_Rep1 is also called lung5-3
# Lung12     Lung13 Lung5_Rep1 Lung5_Rep2 Lung5_Rep3      Lung6 Lung9_Rep1 Lung9_Rep2 
# 71304      81236      98002     105800      97809      89975      87606     139504 

tissueName <- "Lung9_Rep2" # replace this by different tissue names listed above
col_dat <- gem@cell_metadata$rna %>%
  filter(Run_Tissue_name == tissueName) %>%
  as.data.frame()
rownames(col_dat) <- col_dat$cell_ID
locs <- as.data.frame(gem@spatial_locs$raw)
rownames(locs) <- locs$cell_ID
locs <- locs[rownames(col_dat),-3]
col_dat <- cbind(col_dat, locs)
## SCE object
query <- SingleCellExperiment(
  list(counts = gem@expression$rna$raw[, rownames(col_dat)]),
  colData = col_dat 
)
clusts <- quickCluster(query, assay.type = "counts")
query <- computeSumFactors(query, assay.type = "counts", cluster = clusts)
# If having non-positive size factors
if(min(query$sizeFactor) <= 0){
  smallestPositive <- query$sizeFactor[which(sort(query$sizeFactor) > 0)[1]]
  query$sizeFactor[query$sizeFactor <= 0] <- smallestPositive
}
query <- logNormCounts(query, assay.type = "counts")
saveRDS(query, paste0(dat_dir, "data/CosMx_lung/", tissueName, "_SCE.rds"))






















