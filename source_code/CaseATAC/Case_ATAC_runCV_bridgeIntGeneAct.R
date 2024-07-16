# Using gene activity scores as ATAC features to run Seurat bridge integration

.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")


suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

options(future.globals.maxSize = 8000 * 1024^2)


## Seurat Bridge integration ----------------------------------------------------------------
## Loading data (see benchBridge_BMMC_BridgeInt.R for data preparation)
# BMMC scRNAseq reference & bimoadal datasets
obj.rna <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.rna_for_refMap.rds"))

# Store normalised gene activity scores
multiPath <- paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches.geneAct_cor.rds")
if(!file.exists(multiPath)){
  
  obj.multi.batches <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches_cor.rds"))
  bridgeATAC_list <- readRDS(
    paste0(
      dat_dir, 
      "data/NeurIPS2021/bridgeATAC_list_cor.rds"
    )
  )
  
  for(ii in 1:length(obj.multi.batches)){
    
    obj.multi <- obj.multi.batches[[ii]]
    bridgeATAC <- bridgeATAC_list[[ii]]
    
    obj.multi[["ATAC"]] <- CreateAssayObject(
      data = assay(bridgeATAC, "logcounts")
    )
    
    obj.multi.batches[[ii]] <- obj.multi
  }
  
  saveRDS(obj.multi.batches, multiPath)
} else {
  
  obj.multi.batches <- readRDS(multiPath)
}



## CV iteration
# In each outer iteration (indexed by jj), use 1 batch as bridge
# In each inner iteration (indexed by kk), use 1 batch as query
batchNames <- names(obj.multi.batches)
bridgeIntAnn_cv <- vector("list", length(batchNames))
names(bridgeIntAnn_cv) <- batchNames
for(jj in 1:length(batchNames)){
  
  # Batch names
  bridgeBatchName <- batchNames[jj]
  queryBatchNames <- setdiff(batchNames, bridgeBatchName)
  bridgeIntAnn <- vector("list", length(queryBatchNames))
  names(bridgeIntAnn) <- queryBatchNames
  
  # Prepare Bridge
  obj.multi <- obj.multi.batches[[bridgeBatchName]]
  # normalize multiome RNA
  DefaultAssay(obj.multi) <- "RNA"
  obj.multi <- SCTransform(obj.multi, verbose = FALSE)
  # normalize multiome ATAC
  DefaultAssay(obj.multi) <- "ATAC"
  obj.multi[["ATAC"]]$counts <- obj.multi[["ATAC"]]$data
  obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
  
  # Prepare bridge 
  dims.atac <- 1:50
  dims.rna <- 1:50
  DefaultAssay(obj.multi) <-  "RNA"
  DefaultAssay(obj.rna) <- "SCT"
  obj.rna.ext <- PrepareBridgeReference(
    reference = obj.rna,
    bridge = obj.multi,
    reference.reduction = "spca",
    reference.dims = dims.rna,
    normalization.method = "SCT"
  )
  
  for(kk in 1:length(queryBatchNames)){
    
    cat("Sim= ", jj, "/", length(batchNames), "; Fold ", kk, "/", length(queryBatchNames), "\n")
    
    # Query
    queryBatchName <- queryBatchNames[kk]
    obj.atac <- obj.multi.batches[[queryBatchName]]
    DefaultAssay(obj.atac) <- "ATAC"
    
    
    ## Bridge integration
    bridge.anchor <- FindBridgeTransferAnchors(
      extended.reference = obj.rna.ext,
      query = obj.atac,
      reduction = "lsiproject",
      dims = dims.atac
    )
    obj.atac <- MapQuery(
      anchorset = bridge.anchor, 
      reference = obj.rna.ext, 
      query = obj.atac, 
      refdata = list(
        l1 = "celltype.l1",
        l2 = "celltype.l2"
      ),
      reduction.model = "wnn.umap"
    )
    
    ## Store bridge classif results
    bridgeIntAnn[[kk]] <- list(
      metadata = obj.atac@meta.data,
      scores.l1 = obj.atac@assays$prediction.score.l1$data,
      scores.l2 = obj.atac@assays$prediction.score.l2$data
    )
  }
  
  bridgeIntAnn_cv[[bridgeBatchName]] <- bridgeIntAnn
}
saveRDS(bridgeIntAnn_cv, paste0(dat_dir, "output/Case2/bridgeIntAnnGeneAct_cv_cor.rds"))
