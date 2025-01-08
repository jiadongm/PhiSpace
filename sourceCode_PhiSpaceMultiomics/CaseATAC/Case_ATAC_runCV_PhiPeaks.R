### PhiSpace --------------------------------------------------------------------
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"


## Load data (see PrepareData.R for data processing and normalisation)
bridgeATAC_list <- readRDS(paste0(dat_dir,"data/NeurIPS2021/bridgeATAC_TFIDF_list_cor.rds"))
bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))
reference <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.rna_for_refMap_sce.rds"))


## RNA modality: No pseudo-bulking
PhiSpaceAssay <- "logcounts"
YtrainName <- "celltype.l2"

# Tuning
if(F){
  tune_res <- tunePhiSpace(
    reference = reference,
    assayName = PhiSpaceAssay,
    phenotypes = YtrainName,
    tune_nfeat = FALSE,
    regMethod = "PLS",
    Kfolds = 10
  )
}

PhiRes <- PhiSpaceR_1ref(
  reference, 
  query = bridgeRNA_list, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  regMethod = "PLS",
  scale = FALSE
)
rm(bridgeRNA_list); gc()
if(F){
  YrefHat <- PhiRes$YrefHat
  YrefHat_norm <- normPhiScores(YrefHat)
  YrefHats <- list(
    YrefHat = YrefHat,
    YrefHat_norm = YrefHat_norm
  )
  if(F) saveRDS(YrefHats, paste0(dat_dir, "output/Case2/YrefHats.rds"))
}

## ATAC modality
batchNames <- names(bridgeATAC_list)
PhiSpaceAnn_cv <- vector("list", length(batchNames))
names(PhiSpaceAnn_cv) <- batchNames
for(jj in 1:length(batchNames)){
  
  
  cat("Sim= ", jj, "/", length(batchNames), "\n")
  
  # Batch names
  bridgeBatchName <- batchNames[jj]
  queryBatchNames <- setdiff(batchNames, bridgeBatchName)
  
  
  
  PhiSpaceAssay <- "data"
  if(F){
    tune_res <- tunePhiSpace(
      reference = bridgeATAC_list[[bridgeBatchName]],
      assayName = PhiSpaceAssay,
      YY = PhiRes$PhiSpaceScore[[bridgeBatchName]],
      regMethod = "PLS",
      Kfolds = 10
    )
    # ncomp =18, nfeat = 2000
  }
  
  if(T){
    response <- normPhiScores(PhiRes$PhiSpaceScore[[bridgeBatchName]])
  } else {
    resposne <- PhiRes$PhiSpaceScore[[bridgeBatchName]]
  }
  
  PhiRes_ATAC <- PhiSpaceR_1ref(
    reference = bridgeATAC_list[[bridgeBatchName]], 
    query = bridgeATAC_list[queryBatchNames], 
    response = response, 
    PhiSpaceAssay = PhiSpaceAssay, 
    regMethod = "PLS",
    center = T, # Convention in TF-IDF processing is to not center (preserve sparsity)
    scale = FALSE
  )
  
  
  PhiSpaceAnn_cv[[bridgeBatchName]] <- lapply(
    PhiRes_ATAC$PhiSpaceScore,
    function(x){
      list(
        queryScore = x,
        queryScore_norm = normPhiScores(x)
      )
    }
  )
  
}
saveRDS(PhiSpaceAnn_cv, paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_TF-IDF_cor.rds"))