.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")

### Run this on RPC server 'PhiSpace'
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"


## Load data (see PrepareData.R for data processing and normalisation)
PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc_cor.rds"))
PhiSpaceAnnPeaks_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_TF-IDF_cor.rds"))
reference <- colData(readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.rna_for_refMap_sce.rds")))
# Seurat BridgeInt annotations
bridgeIntAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnn_cv_cor.rds"))
bridgeIntAnnGeneAct_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnnGeneAct_cv_cor.rds"))
cellTypeTable <- readRDS(paste0(dat_dir, "output/Case2/cellTypeTable.rds"))
# Cell type dictionaries
{
  refTypes_l2 <- unique(reference$celltype.l2) %>% sort
  refTypes_l1 <- sapply(
    refTypes_l2,
    function(x){
      unique(reference$celltype.l1[reference$celltype.l2 == x])
    }
  )
  ref_lookup <- data.frame(
    l2 = refTypes_l2,
    l1 = refTypes_l1
  )
  ## Query original annotations (ILC is hard to classify)
  queryTypes_l2 <- rownames(cellTypeTable) %>% sort
  queryTypes_l1 <- c(
    "B cell", "Mono/DC", "Mono/DC", "T cell", "T cell",
    "T cell", "T cell", "Mono/DC", "Progenitor cells", "Progenitor cells",
    "Progenitor cells", "Progenitor cells", "ILC", "Progenitor cells", "Progenitor cells",
    "B cell", "NK", "Progenitor cells", "Mono/DC", "B cell",
    "Progenitor cells", "B cell"
  )
  query_lookup <- data.frame(
    l2 = queryTypes_l2,
    l1 = queryTypes_l1
  )
  
  purity <- function(clust, origin){
    
    clustLabs <- unique(clust)
    purity_per_clust <- sapply(
      clustLabs,
      function(x){
        idx <- (clust==x)
        purity_per_clust <- max(table(origin[idx]))/sum(idx)
        clustSize <- max(table(origin[idx]))
        
        return(
          c(
            purity_per_clust,
            clustSize
          )
        )
      }
    )
    
    avePurity <- mean(purity_per_clust[1,])
    weiPurity <- sum(purity_per_clust[1,] * purity_per_clust[2,]/sum(purity_per_clust[2,]))
    
    return(
      list(
        purity_per_clust = purity_per_clust,
        avePurity = avePurity,
        weiPurity = weiPurity
      )
    )
  }
}

## Prepare data for TFIDF SVD clustering
atacOnly_path <- paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches.ATAC.only_cor.rds")
if(!file.exists(atacOnly_path)){
  obj.multi.batches <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches_cor.rds"))
  for(ii in 1:length(obj.multi.batches)){
    
    DefaultAssay(obj.multi.batches[[ii]]) <- "ATAC"
    obj.multi.batches[[ii]] <- DietSeurat(
      obj.multi.batches[[ii]],
      assays = "ATAC"
    )
    
  }
  
  saveRDS(obj.multi.batches, atacOnly_path)
} else {
  
  obj.multi.batches <- readRDS(atacOnly_path)
}


NclusVec <- c(6, 22)
batchNames <- colnames(cellTypeTable)
ncomp <- length(reference$celltype.l2 %>% unique)
Nsim <- length(batchNames)
out <- vector("list", Nsim)
for(sim in 1:Nsim){
  
  cat("sim =", sim, "/", Nsim, "\n")
  
  bridgeBatchName <- batchNames[sim]
  queryBatchNames <- setdiff(batchNames, bridgeBatchName)
  
  
  
  ## Clustering using TFIDF SVD
  comboATAC <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches.ATAC.only_cor.rds"))[queryBatchNames]
  comboATAC <- SeuratObject:::merge.Seurat(
    x = comboATAC[[1]],
    y = comboATAC[2:length(comboATAC)]
  )
  gc()
  comboATAC <- RunTFIDF(comboATAC)
  comboATAC <- FindTopFeatures(comboATAC, min.cutoff = 'q0')
  comboATAC <- RunSVD(comboATAC, n = ncomp + 1)
  dat2clust_TFIDF <- comboATAC@reductions$lsi@cell.embeddings[,2:(ncomp+1)]
  rm(comboATAC); gc()
  
  ## Clustering using PCs
  comboATAC <- do.call(
    cbind,
    readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_list_cor.rds"))[queryBatchNames]
  )
  gc()
  # comboATAC <- scranTransf(comboATAC)
  comboATAC <- logTransf(comboATAC)
  gc()
  dat2clust_PC <- getPC(
    t(assay(comboATAC, "data")), 
    ncomp = ncomp, 
    center = T, 
    scale = F
  )$scores %>%
    as.matrix()
  comboATAC <- colData(comboATAC)
  gc()
  
  
  
  ## Clustering using Seurat scores
  dat2clust_BI <- bridgeIntAnn_cv[[bridgeBatchName]]
  dat2clust_BI <- lapply(
    1:length(dat2clust_BI),
    function(x){
      t(dat2clust_BI[[x]]$scores.l2)
    }
  )
  dat2clust_BI <- do.call(rbind, dat2clust_BI)
  # Gene activity scores version
  dat2clust_BIgeneAct <- bridgeIntAnnGeneAct_cv[[bridgeBatchName]]
  dat2clust_BIgeneAct <- lapply(
    1:length(dat2clust_BIgeneAct),
    function(x){
      t(dat2clust_BIgeneAct[[x]]$scores.l2)
    }
  )
  dat2clust_BIgeneAct <- do.call(rbind, dat2clust_BIgeneAct)
  
  
  ## Clustering using PhiSpace scores
  dat2clust_Phi <- PhiSpaceAnn_cv[[bridgeBatchName]]
  dat2clust_Phi <- lapply(
    1:length(dat2clust_Phi),
    function(x){
      dat2clust_Phi[[x]]$queryScore_norm
    }
  )
  dat2clust_Phi <- do.call(rbind, dat2clust_Phi)
  # Peaks version
  dat2clust_PhiPeaks <- PhiSpaceAnnPeaks_cv[[bridgeBatchName]]
  dat2clust_PhiPeaks <- lapply(
    1:length(dat2clust_PhiPeaks),
    function(x){
      dat2clust_PhiPeaks[[x]]$queryScore_norm
    }
  )
  dat2clust_PhiPeaks <- do.call(rbind, dat2clust_PhiPeaks)
  
  
  ## Quality of clusters
  origin <- comboATAC$cellType %>% as.character()
  origin_l1 <- query_lookup$l1[match(comboATAC$cellType, query_lookup$l2)]
  originClust <- comboATAC$cellType %>% as.numeric()
  originClust_l1 <- origin_l1 %>% as.factor() %>% as.numeric()
 
  
  clustRes_sim <- lapply(
    1:length(NclusVec),
    function(jj){
      
      Nclust <- NclusVec[jj]
      
      TFIDFkmeans <- kmeans(dat2clust_TFIDF, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      PCkmeans <- kmeans(dat2clust_PC, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      Phikmeans <- kmeans(dat2clust_Phi, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      PhiPeakskmeans <- kmeans(dat2clust_PhiPeaks, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      BIkmeans <- kmeans(dat2clust_BI, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      BIgeneActkmeans <- kmeans(dat2clust_BIgeneAct, centers = Nclust, iter.max = 500, nstart = 50, algorithm = "Lloyd")$cluster
      
      
      # purTFIDFl1_res <- purity(TFIDFkmeans, origin_l1)
      # purPCl1_res <- purity(PCkmeans, origin_l1)
      # purPhil1_res <-  purity(Phikmeans, origin_l1)
      # purBIl1_res <- purity(BIkmeans, origin_l1)
      # purTFIDF_res <- purity(TFIDFkmeans, origin)
      # purPC_res <- purity(PCkmeans, origin)
      # purPhi_res <-  purity(Phikmeans, origin) 
      # purBI_res <-  purity(BIkmeans, origin) 
      
      list(
        TFIDFkmeans = TFIDFkmeans,
        PCkmeans = PCkmeans,
        Phikmeans = Phikmeans,
        PhiPeakskmeans = PhiPeakskmeans,
        BIkmeans = BIkmeans,
        BIgeneActkmeans = BIgeneActkmeans
      )
      
      # c(
      #   # ARI
      #   ARI_TFIDF = fossil::adj.rand.index(TFIDFkmeans, originClust),
      #   ARI_PC = fossil::adj.rand.index(PCkmeans, originClust),
      #   ARI_Phi = fossil::adj.rand.index(Phikmeans, originClust),
      #   ARI_BI = fossil::adj.rand.index(BIkmeans, originClust),
      #   # ARI l1
      #   ARI_TFIDF_l1 = fossil::adj.rand.index(TFIDFkmeans, originClust_l1),
      #   ARI_PC_l1 = fossil::adj.rand.index(PCkmeans, originClust_l1),
      #   ARI_Phi_l1 = fossil::adj.rand.index(Phikmeans, originClust_l1),
      #   ARI_BI_l1 = fossil::adj.rand.index(BIkmeans, originClust_l1),
      #   # purity
      #   purTFIDF_ave = purTFIDF_res$avePurity,
      #   purTFIDF_wei = purTFIDF_res$weiPurity,
      #   purPC_ave = purPC_res$avePurity,
      #   purPC_wei = purPC_res$weiPurity,
      #   purPhi_ave = purPhi_res$avePurity,
      #   purPhi_wei = purPhi_res$weiPurity,
      #   purBI_ave = purBI_res$avePurity,
      #   purBI_wei = purBI_res$weiPurity,
      #   # purity l1
      #   purTFIDFl1_ave = purTFIDFl1_res$avePurity,
      #   purTFIDFl1_wei = purTFIDFl1_res$weiPurity,
      #   purPCl1_ave = purPCl1_res$avePurity,
      #   purPCl1_wei = purPCl1_res$weiPurity,
      #   purPhil1_ave = purPhil1_res$avePurity,
      #   purPhil1_wei = purPhil1_res$weiPurity,
      #   purBIl1_ave = purBIl1_res$avePurity,
      #   purBIl1_wei = purBIl1_res$weiPurity
      # )
    }
  )
  
  
  out[[sim]] <- list(
    origin = origin,
    origin_l1 = origin_l1,
    originClust = originClust,
    originClust_l1,
    clustRes_sim = clustRes_sim
  )
  
  
}
saveRDS(out, paste0(dat_dir, "output/Case2/clus_res_cor.rds"))























#################### Old: do not run ####################
# if(F){
#   devtools::install_github("jiadongm/PhiSpace/pkg")
#   
#   library(PhiSpace)
#   library(ggplot2)
#   library(dplyr)
#   library(magrittr)
#   library(ggpubr)
#   library(tidyr)
#   
#   rm(list=ls()); gc()
#   setwd("~/PhiSpaceR")
#   dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
#   
#   
#   ## Load data (see PrepareData.R for data processing and normalisation)
#   bridgeATAC_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_list.rds"))
#   PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc.rds"))
#   
#   reference <- colData(readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.rna_for_refMap_sce.rds")))
#   # Bridge (Choose a bridge that yields comparable performances of Seurat and PhiSpace)
#   batchNames <- names(bridgeATAC_list)
#   bridgeBatchName <- "s1d1"
#   queryBatchNames <- setdiff(batchNames, bridgeBatchName)
#   # Seurat BridgeInt annotations
#   bridgeIntAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnn_cv.rds"))
#   originAnn_l <- readRDS(paste0(dat_dir, "output/Case2/originAnn_l.rds"))
#   cellTypeTable <- readRDS(paste0(dat_dir, "output/Case2/cellTypeTable.rds"))
#   
#   refTypes_l2 <- unique(reference$celltype.l2) %>% sort
#   refTypes_l1 <- sapply(
#     refTypes_l2,
#     function(x){
#       unique(reference$celltype.l1[reference$celltype.l2 == x])
#     }
#   )
#   ref_lookup <- data.frame(
#     l2 = refTypes_l2,
#     l1 = refTypes_l1
#   )
#   ## Query original annotations (ILC is hard to classify)
#   queryTypes_l2 <- rownames(cellTypeTable) %>% sort
#   queryTypes_l1 <- c(
#     "B cell", "Mono/DC", "Mono/DC", "T cell", "T cell",
#     "T cell", "T cell", "Mono/DC", "Progenitor cells", "Progenitor cells",
#     "Progenitor cells", "Progenitor cells", "ILC", "Progenitor cells", "Progenitor cells",
#     "B cell", "NK", "Progenitor cells", "Mono/DC", "B cell",
#     "Progenitor cells", "B cell"
#   )
#   query_lookup <- data.frame(
#     l2 = queryTypes_l2,
#     l1 = queryTypes_l1
#   )
#   
#   
#   
#   # Subsample for clustering
#   Nsim <- 50
#   purity <- function(clust, origin){
#     
#     clustLabs <- unique(clust)
#     purity_per_clust <- sapply(
#       clustLabs,
#       function(x){
#         idx <- (clust==x)
#         purity_per_clust <- max(table(origin[idx]))/sum(idx)
#         clustSize <- max(table(origin[idx]))
#         
#         return(
#           c(
#             purity_per_clust,
#             clustSize
#           )
#         )
#       }
#     )
#     
#     avePurity <- mean(purity_per_clust[1,])
#     weiPurity <- sum(purity_per_clust[1,] * purity_per_clust[2,]/sum(purity_per_clust[2,]))
#     
#     return(
#       list(
#         purity_per_clust = purity_per_clust,
#         avePurity = avePurity,
#         weiPurity = weiPurity
#       )
#     )
#   }
#   ncomp <- length(reference$celltype.l2 %>% unique)
#   out <- vector("list", Nsim)
#   NclusVec <- c(5, 6, 10, 15, 22, 30)
#   for(sim in 1:Nsim){
#     
#     cat("sim =", sim, "/", Nsim, "\n")
#     
#     set.seed(53522 + sim)
#     
#     comboATAC <- lapply(
#       1:length(queryBatchNames),
#       function(i){
#         
#         queryBatchName <- queryBatchNames[i]
#         queryBatch <- bridgeATAC_list[[queryBatchName]]
#         
#         idx <- sample(1:ncol(queryBatch), 1000)
#         
#         return(
#           list(
#             queryBatch[, idx],
#             idx
#           )
#         )
#         
#       }
#     )
#     
#     idx_list <- lapply(
#       comboATAC,
#       function(x){
#         x[[2]]
#       }
#     )
#     
#     comboATAC <- lapply(
#       comboATAC,
#       function(x){
#         x[[1]]
#       }
#     )
#     comboATAC <- do.call(cbind, comboATAC)
#     
#     
#     
#     
#     ## Clustering using PCs
#     comboATAC <- scranTransf(comboATAC)
#     if(T){
#       dat2clust_PC <- getPC(t(assay(comboATAC, "logcounts")), ncomp = ncomp, center = T, scale = F)$scores
#     } else {
#       dat2clust_PC <- getPC(assay(comboATAC), ncomp = ncomp, center = F, scale = F)$loadings
#     }
#     
#     
#     ## Clustering using PhiSpace scores
#     dat2clust <- PhiSpaceAnn_cv[[bridgeBatchName]]
#     dat2clust <- lapply(
#       1:length(dat2clust),
#       function(x){
#         dat2clust[[x]]$queryScore_norm[idx_list[[x]],]
#       }
#     )
#     dat2clust_Phi <- do.call(rbind, dat2clust)
#     
#     
#     ## Quality of clusters
#     origin_l1 <- query_lookup$l1[match(comboATAC$cellType, query_lookup$l2)]
#     originClust <- comboATAC$cellType %>% as.numeric()
#     originClust_l1 <- origin_l1 %>% as.factor() %>% as.numeric()
#     origin <- comboATAC$cellType %>% as.character()
#     
#     out[[sim]] <- vector("list", length(NclusVec))
#     for(jj in 1:length(NclusVec)){
#       
#       Nclust <- NclusVec[jj]
#       
#       PCkmeans <- kmeans(dat2clust_PC, centers = Nclust, iter.max = 100, nstart = 50)$cluster
#       Phikmeans <- kmeans(dat2clust_Phi, centers = Nclust, iter.max = 100, nstart = 50)$cluster
#       
#       
#       purPCl1_res <- purity(PCkmeans, origin_l1)
#       purPhil1_re <-  purity(Phikmeans, origin_l1) 
#       purPC_res <- purity(PCkmeans, origin)
#       purPhi_re <-  purity(Phikmeans, origin) 
#       
#       out[[sim]][[jj]] <- c(
#         ARI_PC = fossil::adj.rand.index(PCkmeans, originClust),
#         ARI_Phi = fossil::adj.rand.index(Phikmeans, originClust),
#         ARI_PC_l1 = fossil::adj.rand.index(PCkmeans, originClust_l1),
#         ARI_Phi_l1 = fossil::adj.rand.index(Phikmeans, originClust_l1),
#         purPC_ave = purPC_res$avePurity,
#         purPC_wei = purPC_res$weiPurity,
#         purPhi_ave = purPhi_re$avePurity,
#         purPhi_wei = purPhi_re$weiPurity,
#         purPCl1_ave = purPCl1_res$avePurity,
#         purPCl1_wei = purPCl1_res$weiPurity,
#         purPhil1_ave = purPhil1_re$avePurity,
#         purPhil1_wei = purPhil1_re$weiPurity
#       )
#       
#       if(F){
#         plotSankey3(PCkmeans, origin_l1, Phikmeans)
#       }
#       
#     }
#     
#   }
#   # if(T){
#   #   saveRDS(out, paste0(dat_dir, "output/Case2/clus_res_genAct.rds"))
#   # } else {
#   #   saveRDS(out, paste0(dat_dir, "output/Case2/clus_res.rds"))
#   # }
# }


































