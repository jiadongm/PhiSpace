

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

source("~/PhiSpaceR/CaseVisium/utils.R")
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

PhiAssay <- "log1p"


# Reference
YtrainName <- refLabName <- "ann_finest_level"
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
if(!file.exists(refPath)){
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
  
} else {
  
  reference <- qread(refPath)
}

impScPath <- paste0(
  dat_dir,
  "output/Case3/refImpScores.rds"
)
if(!file.exists(impScPath)){
  
  tuneRes <- tunePhiSpace(
    reference = reference,
    assayName = "log1p",
    phenotypes = YtrainName,
    tune_ncomp = F,
    tune_nfeat = F
  )
  impScores <- tuneRes$impScores
  saveRDS(impScores, impScPath)
} else {
  
  impScores <- readRDS(impScPath)
}

# Feature selection 
selectedFeat <- selectFeat(impScores, 500)$selectedFeat




tissueNames <- tissueNames_Visium
visListPath <- paste0(
  dat_dir, "data/Visium_NSCLC/allTissues.qs"
)
qu_list <- qread(visListPath)
qu_list <- lapply(
  qu_list, zeroFeatQC
)


PhiRes <- PhiSpaceR_1ref(
  reference, 
  qu_list, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS", 
  center = T,
  scale = F,
  assay2rank = "counts"
)


PhiResPath <- paste0(
  dat_dir,
  "output/Case3/PhiSpace/combo_PhiRes.qs"
)
qsave(PhiRes, PhiResPath)













