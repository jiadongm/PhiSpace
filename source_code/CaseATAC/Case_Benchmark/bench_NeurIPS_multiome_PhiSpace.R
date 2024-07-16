.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"


################ PhiSpace ##################
suppressPackageStartupMessages(library(PhiSpace))

# Load data: 13 batches of scRNA-seq
bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))

YtrainName <- "cellType"
allBatchNames <- names(bridgeRNA_list)
SR_ann_res <- vector(
  "list", length(allBatchNames)
) %>%
  `names<-`(allBatchNames)
for(refIdx in 1:length(allBatchNames)){
  
  # Define reference dataset
  refBatchName <- allBatchNames[refIdx]
  reference <- bridgeRNA_list[[refBatchName]]
  
  # Remaining batches are queries
  queryBatchNames <- setdiff(
    allBatchNames, refBatchName
  )
  
  temp_res <- vector(
    "list", length(queryBatchNames)
  ) %>%
    `names<-`(queryBatchNames)
  for(queryIdx in 1:length(queryBatchNames)){
    
    # Define query dataset
    queryBatchName <- queryBatchNames[queryIdx]
    query <- bridgeRNA_list[[queryBatchName]]
    
    classOriginal <- colData(query)[, YtrainName] %>% as.character()
    
    # PhiSpace prediction
    PhiAssay <- "logcounts"
    regMethod <- "PLS"
    time_ <-
      system.time(
        PhiRes <- PhiSpaceR_1ref(
          reference = reference, 
          query = query,
          phenotypes = YtrainName,
          PhiSpaceAssay = PhiAssay,
          regMethod = regMethod
        )
      )
    cat("PhiSpace took", round(time_[3]), "seconds to compute.", "\n")
    
    # Class err
    classSR <- PhiRes$PhiSpaceScore %>% 
      normPhiScores() %>%
      getClass()
    errSR <- PhiSpace:::classErr(classSR, classOriginal)$err
    cat("PhiSpace errors =", errSR, '\n')
    
    # Store results
    temp_res[[queryBatchName]] <- list(
      PhiRes = PhiRes,
      errSR = errSR,
      runTime = time_
    )
    
  }
  
  # Store results
  SR_ann_res[[refBatchName]] <- temp_res
}

saveRDS(
  SR_ann_res,
  paste0(
    dat_dir,
    "output/scRNAseq/NeurIPS2021_multiome_Phi_ann_res_cor.rds"
  )
)