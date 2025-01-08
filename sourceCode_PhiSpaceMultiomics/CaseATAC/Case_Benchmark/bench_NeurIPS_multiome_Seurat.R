.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"


################ Seurat ##################
suppressPackageStartupMessages(library(Seurat))

# Load data: 13 batches of scRNA-seq
dataPath <- paste0(
  dat_dir,
  "data/NeurIPS2021/bridgeRNA_list_Seu_cor.rds"
)
if(!file.exists(dataPath)){
  
  suppressPackageStartupMessages(library(SingleCellExperiment))
  suppressPackageStartupMessages(library(SeuratObject))
  
  bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))
  
  for(idx in 1:length(bridgeRNA_list)){
    
    
    sce <- bridgeRNA_list[[idx]]
    
    seu <- CreateSeuratObject(
      counts = assay(sce, "counts"),
      meta.data = colData(sce) %>% as.data.frame()
    )
    
    seu[["RNA"]]$data <- assay(sce, "logcounts")
    
    bridgeRNA_list[[idx]] <- seu
  }
  
  saveRDS(bridgeRNA_list, dataPath)
} else {
  
  obj.multi.batches <- readRDS(dataPath)
}



YtrainName <- "cellType"
allBatchNames <- names(obj.multi.batches)
SR_ann_res <- vector(
  "list", length(allBatchNames)
) %>%
  `names<-`(allBatchNames)
for(refIdx in 1:length(allBatchNames)){
  
  # Define reference dataset
  refBatchName <- allBatchNames[refIdx]
  reference <- obj.multi.batches[[refBatchName]]
  
  # Remaining batches are queries
  queryBatchNames <- setdiff(
    allBatchNames, refBatchName
  )
  
  
  ## Seurat pre-processing
  DefaultAssay(reference) <- "RNA"
  reference <- FindVariableFeatures(reference)
  reference <- ScaleData(reference)
  reference <- RunPCA(reference)
  reference <- FindNeighbors(reference, dims = 1:30)
  reference <- FindClusters(reference)
  
  
  
  temp_res <- vector(
    "list", length(queryBatchNames)
  ) %>%
    `names<-`(queryBatchNames)
  for(queryIdx in 1:length(queryBatchNames)){
    
    # Define query dataset
    queryBatchName <- queryBatchNames[queryIdx]
    query <- obj.multi.batches[[queryBatchName]]
    
    classOriginal <- query@meta.data[, YtrainName] %>% as.character()
    
    # Seurat prediction
    time_ <-
      system.time(
        ref.anchors <- FindTransferAnchors(
          reference = reference,
          query = query,
          dims = 1:30,
          reference.reduction = "pca"
        )
      )
    time__ <- system.time(
      predictions <- TransferData(
        anchorset = ref.anchors,
        refdata = reference$cellType,
        dims = 1:30
      )
    )
    computTime <- time_ + time__
    
    # Class err
    classSR <- predictions$predicted.id
    errSR <- PhiSpace:::classErr(classSR, classOriginal)$err
    cat("Seurat errors =", errSR, '\n')
    
    # Store results
    temp_res[[queryBatchName]] <- list(
      predictions = predictions,
      errSR = errSR,
      runTime = computTime
    )
    
  }
  
  # Store results
  SR_ann_res[[refBatchName]] <- temp_res
}

saveRDS(
  SR_ann_res,
  paste0(
    dat_dir,
    "output/scRNAseq/NeurIPS2021_multiome_Seu_ann_res_cor.rds"
  )
)








































