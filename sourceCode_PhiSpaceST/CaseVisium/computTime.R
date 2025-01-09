# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
source("/data/projects/punim0613/JiaDong/PhiSpace/CaseVisium/utils.R")
PhiAssay <- "log1p"


# Reference
YtrainName <- refLabName <- "ann_finest_level"
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
reference <- qread(refPath)
reference <- zeroFeatQC(reference)
# Query
tissueName <- tissueNames_Visium[1]
visListPath <- paste0(
  dat_dir, "data/Visium_NSCLC/allTissues.qs"
)
qu_list <- qread(visListPath)
qu_list <- lapply(
  qu_list, zeroFeatQC
)
query <- qu_list[[tissueName]]
rm(qu_list); gc()


#################### PhiSpace ###############################################
tik <- Sys.time()
# Feature selection
impScPath <- paste0(
  dat_dir,
  "output/Case3/refImpScores.rds"
)
impScores <- readRDS(impScPath)
selectedFeat <- selectFeat(impScores, 500)$selectedFeat

# PhiSpace
PhiRes <- PhiSpaceR_1ref(
  reference, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS", 
  center = T,
  scale = F,
  assay2rank = "counts"
)
toc <- Sys.time()
toc - tik

#################### RCTD ###############################################
# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(spacexr))
suppressPackageStartupMessages(library(Matrix))


# Query
VisTissueName <- tissueName
visPath <- paste0(
  paste0(
    dat_dir,
    "data/Visium_NSCLC/",
    VisTissueName,
    "/sce.rds"
  )
)
qu_vis <- readRDS(visPath)
qu_vis <- zeroFeatQC(qu_vis)


# Feature selection
impScPath <- paste0(
  dat_dir,
  "output/Case3/refImpScores.rds"
)
impScores <- readRDS(impScPath)
selectedFeat <- selectFeat(impScores, 500)$selectedFeat

cGenes <- intersect(
  selectedFeat, rownames(qu_vis)
)

# run RCTD
tik <- Sys.time()
refRCTD <- Reference(
  counts = assay(reference[cGenes,], "counts"),
  cell_types = colData(reference)[,YtrainName] %>% `names<-`(colnames(reference))
)
quRCTD <- SpatialRNA(
  coords = colData(qu_vis)[,c("x","y")],
  counts = assay(qu_vis[cGenes,], "counts")
)
myRCTD <- create.RCTD(
  quRCTD,
  reference = refRCTD,
  max_cores = 4
)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
toc <- Sys.time()
toc - tik
























