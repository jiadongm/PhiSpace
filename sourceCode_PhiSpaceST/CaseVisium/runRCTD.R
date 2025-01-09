
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
suppressPackageStartupMessages(library(spacexr))
suppressPackageStartupMessages(library(Matrix))

source("~/PhiSpaceR/Case_CosMx/utils.R")
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

PhiAssay <- "log1p"



# Reference
YtrainName <- refLabName <- "ann_finest_level"
# reference <- readRDS(paste0(dat_dir, "output/Case3/Azimuth2.0_PB_ref.rds"))
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
reference <- qread(refPath)


tissueNames <- c(
  "D2_1", "D2_2", 
  "P10_T1", "P10_B1", 
  "P11_T3", "P11_B1", 
  "P15_B1", "P15_T2", 
  "P16_B1", "P16_T1",
  "P17_B1", "P17_T1", 
  "P19_B2", "P19_T1", 
  "P24_B1", "P24_T1", 
  "P25_B2", "P25_T1"
)


for(ii in 1:length(tissueNames)){
  
  # Query
  VisTissueName <- tissueNames[[ii]]
  
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
  
  # RCTD
  RCTDpath <- paste0(
    dat_dir,
    "output/Case3/RCTD/NSCLC_Visium",
    VisTissueName,
    "_RCTDres.rds"
  )
  cGenes <- intersect(
    selectedFeat,
    rownames(qu_vis)
  )
  
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
  
  saveRDS(myRCTD, RCTDpath)
}
