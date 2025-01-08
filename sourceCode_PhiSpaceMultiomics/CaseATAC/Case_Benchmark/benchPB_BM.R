.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")

############################## Pseudo-bulk: Seurat bone marrow ###############################
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(SingleR))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

bmmc_path <- paste0(
  dat_dir,
  "output/benchPB/BM_query020000.rds"
)

if(!file.exists(bmmc_path)){
  # SeuratData::InstallData("bmcite")
  library(SeuratData)
  query0 <- bmcite
  query0 <- SingleCellExperiment(
    assays = list(counts = query0@assays$RNA$counts),
    colData = query0@meta.data 
  )
  table(query0$celltype.l2)
  dim(query0)
  query0 <- zeroFeatQC(query0, "counts")
  dim(query0)
  saveRDS(query0, bmmc_path)
} else {
  
  query0 <- readRDS(bmmc_path)
}
YtrainName <- 'celltype.l2'


## CV evaluation of error rates
Kfolds <- 5
clTable <- table(colData(query0)[,YtrainName])
sort(clTable <- table(colData(query0)[,YtrainName]))
seed_vec <- 5202056:5202075
re_list <- Errs_list <- vector("list", length(seed_vec))
# jj <- 1
for(jj in 1:length(seed_vec)){
  seed <- seed_vec[jj]
  
  
  ## CV computation of classification errors
  # Create stratified per-class data partition
  set.seed(seed)
  splits <- vector("list", length(clTable))
  for(clIdx in 1:length(clTable)){
    cl <- names(clTable)[clIdx]
    allIdx <- which(colData(query0)[,YtrainName] == cl)
    permuIdx <- sample(allIdx) # random permutation of allIdx
    foldSize <- floor(length(allIdx)/Kfolds) # input of split f below
    splits[[clIdx]] <- PhiSpace:::split2(permuIdx, Kfolds)
  }
  splitIdx <- vector("list", Kfolds)
  for(kk in 1:Kfolds){
    splitIdx[[kk]] <- sapply(1:length(splits),
                             function(x){
                               splits[[x]][[kk]]
                             }) %>% unlist()
  }
  
  # Main loop for CV
  re_jj <- vector("list", Kfolds)
  Errs <- matrix(NA, 4, Kfolds)
  rownames(Errs) <- c("SRover", "SRbal", 
                      "NPover", "NPbal")
  for(kk in 1:Kfolds){
    # kk <- 1
    cat("sim=", jj, "; Fold idex= ", kk, "/", Kfolds, "\n")
    
    queryIdx <- splitIdx[[kk]]
    refIdx <- setdiff(1:ncol(query0), queryIdx)
    reference <- query0[,refIdx]
    query <- query0[,queryIdx]
    # dim(reference); dim(query)
    
    
    ## Pseudo-bulk
    time_ <-
      system.time(
        reference <- pseudoBulk(reference,
                                phenotypes = YtrainName,
                                clusterid = YtrainName,
                                assay = "counts",
                                nPool = 15,
                                resampSizes = 100,
                                seed = seed+kk)
      )
    cat("PB aggregation took", round(time_[3]), "seconds to compute.", "\n")
    
    colData(reference)[,YtrainName] <- getClass(
      reducedDim(reference)
    ) 
    
    
    ## Scran transform
    reference <- computeSumFactors(
      reference,
      cluster = quickCluster(reference, assay.type = 'data'),
      assay.type = 'data'
    )
    reference <- logNormCounts(reference, assay.type = 'data')
    query <- computeSumFactors(query, cluster = quickCluster(query))
    query <- logNormCounts(query)
    
    # reference <- RankTransf(reference, "logcounts")
    # query <- RankTransf(query, "logcounts")
    
    
    #################################### SingleR ##############################
    singleRassay <- 'logcounts'
    # singleRassay <- 'rank'
    sr_train <- trainSingleR(assay(reference, singleRassay),
                             labels = colData(reference)[, YtrainName])
    time_ <-
      system.time(
        sr_re <- classifySingleR(assay(query, singleRassay),
                                 sr_train,
                                 fine.tune = T)
      )
    cat("SingleR took", round(time_[3]), "seconds to compute.", "\n")
    classSR <- sr_re$labels
    
    
    ################################## PhiSpace ####################################
    PhiAssay <- 'logcounts'
    # PhiAssay <- 'rank'
    regMethod <- "PLS"
    
    if(F){
      tuneRes <- tunePhiSpace(
        reference = reference,
        phenotypes = YtrainName,
        regMethod = "PLS"
      )
    }
    
    PhiRes <- PhiSpaceR_1ref(
      reference = reference, 
      query = query,
      phenotypes = YtrainName,
      PhiSpaceAssay = PhiAssay,
      regMethod = regMethod
    )
    
    PhiSc_norm <- normPhiScores(
      PhiRes$PhiSpaceScore
    )
    
    classOriginal <- colData(query)[, YtrainName]
    errSR <- PhiSpace:::classErr(classSR, classOriginal)$err
    errPhi <- PhiSpace:::classErr(getClass(PhiSc_norm), classOriginal)$err
    cat("SingleR errors =", errSR, '\n',
        "PhiSpace errors =", errPhi, '\n'
    )
    
    
    ## Store pred scores (eg to compute other errors)
    Errs[,kk] <- c(errSR, errPhi)
    
    re_jj[[kk]] <- PhiRes
    
  }
  
  
  re_list[[jj]] <- re_jj
  Errs_list[[jj]] <- Errs 
  
}

re_list_path <- paste0(
  dat_dir,
  "output/scRNAseq/BMMC_re_list.rds"
)
if(!file.exists(re_list_path)){
  saveRDS(
    re_list, re_list_path
  )
}

errs_list_path <- paste0(
  dat_dir,
  "output/scRNAseq/BMMC_Errs_list.rds"
)
if(!file.exists(errs_list_path)){
  saveRDS(
    Errs_list, errs_list_path
  )
}



