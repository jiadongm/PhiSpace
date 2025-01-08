# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects


### Process Azimuth Lung Reference
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
# Azimuth Lung v2
# https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
reference <- readRDS(paste0(dat_dir, "data/LungRef/AzimuthLung2.0.rds"))
rowdat <- reference@assays$RNA@meta.features

# Look at annotations
reference$cell_type %>% table() %>% length()
reference$ann_finest_level %>% table() %>% names()
reference$disease %>% table()
reference$ann_coarse_for_GWAS_and_modeling %>% table() %>% names()
plotSankey3(colData(reference)[,YtrainName] %>% as.character(),
            reference$ann_level_4 %>% as.character(),
            reference$ann_finest_level %>% as.character())


reference <- SingleCellExperiment(
  list(counts = reference@assays$RNA$counts),
  colData = reference@meta.data
)

all(rownames(reference) == rownames(rowdat))
rownames(reference) <- rowdat$feature_name

qsave(reference, paste0(dat_dir, "data/LungRef/AzimuthLung2.0_sce.qs"))




## Subsetting reference
# Reference
YtrainName <- refLabName <- "ann_finest_level"
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
reference <- qread(
  paste0(
    dat_dir,
    "data/LungRef/AzimuthLung2.0_sce.qs"
  )
)
ref_sub <- subsample(
  reference,
  key = "ann_finest_level",
  proportion = 0.1,
  minCellNum = 50
)
ref_sub <- zeroFeatQC(ref_sub)
ref_sub <- logTransf(
  ref_sub,
  use_log1p = TRUE,
  targetAssay = "log1p"
)
qsave(ref_sub, refPath)








## CosMx: Convert original Giotto objects to SCE objects

# Tried: "Lung5_Rep2"
# "Lung13" is special as no cell type map correspond perfectly to cancer niche (low proportion of cancer cells)
# Lung12
# Lung6 Squamous cell

library(Giotto)
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
load(paste0(dat_dir, "data/CosMx_lung/SMI_Giotto_Object.RData")) 


## Separate lung tissues
table(gem@cell_metadata$rna$Run_Tissue_name)
# ! note that Lung5_Rep1 is also called lung5-3
# Lung12     Lung13 Lung5_Rep1 Lung5_Rep2 Lung5_Rep3      Lung6 Lung9_Rep1 Lung9_Rep2 
# 71304      81236      98002     105800      97809      89975      87606     139504 

tissueName <- "Lung9_Rep2"
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






















