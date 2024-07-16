#### Cross-validation: use each of 13 batches as bridge
devtools::install_github("jiadongm/PhiSpace/pkg")

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"



if(T){
  # Load Seurat BridgeInt and PhiSpace annotations
  # bridgeIntAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnn_cv_cor.rds"))
  bridgeIntAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnnGeneAct_cv_cor.rds"))
  
  if(T){
    
    
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc_cor.rds"))
    
    
    
  } else {
    
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc.rds"))
    
    ## Using gene activity scores
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct.rds"))
    
    # Use normPhiScore as response (recommended)
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc.rds"))
    # Both PCA, not good
    PhiSpaceAnn_cv <- readRDS(
      paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_normPhiSc_method=PCA.rds")
    )
    
    # PB (not good)
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_genAct_PB_refOnly.rds"))
    
    ## Using peaks
    # Used normPhiSC
    PhiSpaceAnn_cv <- readRDS(
      paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_TF-IDF.rds")
    )
    
    # Default settings for RNA and ATAC
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv.rds"))
    # Alternative: SCT norm for RNA
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_ncomp_SCT.rds"))
    # Default for RNA, tunned ATAC
    PhiSpaceAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/PhiSpaceAnn_cv_tunedATAC.rds"))
  }
  
  
  ### Classification accuracy: convert cell types to major cell types
  ## Reference
  originAnn_l <- readRDS(paste0(dat_dir, "output/Case2/originAnn_l.rds"))
  cellTypeTable <- readRDS(paste0(dat_dir, "output/Case2/cellTypeTable.rds"))
  # YrefHats <- readRDS(paste0(dat_dir, "output/Case2/YrefHats.rds"))
  if(F){
    bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list.rds"))
    # Cell types in each file
    originAnn_l <-
      sapply(
        bridgeRNA_list,
        function(x){
          colData(x)
        }
      )
    saveRDS(originAnn_l, paste0(dat_dir, "output/Case2/originAnn_l.rds"))
  }
  reference <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.rna_for_refMap_sce.rds"))
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
  query_origin_lookup <- data.frame(
    l2 = queryTypes_l2,
    l1 = queryTypes_l1
  )
  
  
  nBatches <- length(bridgeIntAnn_cv)
  bridgeNames <- names(bridgeIntAnn_cv)
  errs_l <- vector("list", nBatches*(nBatches - 1))
  for(jj in 1:length(bridgeNames)){
    
    bridgeName <- bridgeNames[jj]
    bridgeIntAnn_l <- bridgeIntAnn_cv[[bridgeName]]
    PhiSpaceAnn_l <- PhiSpaceAnn_cv[[bridgeName]]
    queryNames <- names(bridgeIntAnn_l)
    
    for(kk in 1:length(queryNames)){
      
      queryName <- queryNames[kk]
      
      # Translate Seurat bridgeInt labels
      bridgeIntAnn <- bridgeIntAnn_l[[queryName]]$metadata$predicted.l2
      bridgeIntAnn_l1 <- ref_lookup$l1[match(bridgeIntAnn, ref_lookup$l2)]
      
      nNAs <- bridgeIntAnn_l1 %>% is.na() %>% sum()
      if(nNAs) stop("NAs appeared in converting l2 to l1 cell type lables.")
      
      # Translate PhiSpace labels
      PhiSpaceAnn <- PhiSpaceAnn_l[[queryName]]$queryScore_norm %>% getClass
      PhiSpaceAnn_l1 <- ref_lookup$l1[match(PhiSpaceAnn, ref_lookup$l2)]
      
      nNAs <- PhiSpaceAnn_l1 %>% is.na() %>% sum()
      if(nNAs) stop("NAs appeared in converting l2 to l1 cell type lables.")
      
      # Original annotations
      originAnn <- originAnn_l[[queryName]]$cellType
      originAnn_l1 <- query_origin_lookup$l1[
        match(originAnn, query_origin_lookup$l2)
      ]
      
      # Some cell types do not have good correspondence to l1 cell types
      types2move <- c("ILC")
      idx <- !(originAnn %in% types2move)
      
      
      ## PhiSpace errs by other methods
      
      
      errs_ind <- (jj-1)*nBatches + kk
      errs_l[[errs_ind]] <-
        data.frame(
          errs = c(
            PhiSpace:::classErr(bridgeIntAnn_l1[idx], originAnn_l1[idx])$err,
            PhiSpace:::classErr(PhiSpaceAnn_l1[idx], originAnn_l1[idx])$err
          ),
          typeOfErr = c("overall", "balanced", "overall", "balanced"),
          bridge = rep(bridgeName, 4),
          query = rep(queryName, 4),
          method = c("SeuratBridge", "SeuratBridge", "PhiSpace", "PhiSpace")
        )  
      
    }
    
    if(F){
      plotSankey3(bridgeIntAnn, originAnn, PhiSpaceAnn)
    }
  }
  
  ## Visualise errors
  errs_long <- do.call(rbind, errs_l)
  # BridgeInt vs PhiSpace (just shown balanced error here)
  errs_long %>%
    filter(typeOfErr == "balanced") %>%
    ggplot(aes(method, errs, fill = method)) +
    geom_boxplot() +
    scale_fill_manual(
      values = c(
        "#E41A1C", "#377EB8"
      )
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("Error rate") +
    guides(
      fill = guide_legend(
        title = element_blank()
      )
    )
}
ggsave(
  "classErr.pdf",
  width = 5,
  height = 5
)

# Line plots
outPlots <- vector("list", nBatches)
for(i in 1:nBatches){
  
  bridgeName <- bridgeNames[i]
  outPlots[[i]] <-
    errs_long %>%
    filter(bridge == bridgeName) %>%
    ggplot(aes(query, errs, colour = method, shape = typeOfErr)) +
    geom_point(size = 2) +
    theme(axis.text.x = element_text(angle = 45),
          axis.title = element_blank()) +
    ggtitle(paste0("Bridge = ", bridgeName))
}
ggarrange(plotlist = outPlots, common.legend = T)















######################## OLD ############################
### Analyse results ------------------------
## Dictionaries
if(F){
  reference0 <- readRDS("data/NeurIPS2021/obj.rna_for_refMap_sce.rds")
  refTypeDict <- 
    reference0@colData[,c("celltype.l1", "celltype.l2")] %>%
    as.data.frame() %>%
    distinct(celltype.l2, .keep_all = T)
  rownames(refTypeDict) <- NULL
  rm(reference0); gc()
  
  bridgeRNA_list <- readRDS("data/NeurIPS2021/bridgeRNA_list.rds")
  toDelete <- c("ILC")
  for(i in 1:length(bridgeRNA_list)){
    temp <- bridgeRNA_list[[i]]
    if(any(temp$cellType %in% toDelete)){
      bridgeRNA_list[[i]] <- temp[,!(temp$cellType %in% toDelete)]
    }
  }
  labs <- lapply(bridgeRNA_list,
                 function(x){
                   x$cellType
                 })
  labs <- do.call(c, labs)
  bridgeTypeDict <- data.frame(
    original = sort(as.character(unique(labs))),
    celltype.l1 = c("B cell", "Mono/DC", "Mono/DC", "T cell", 
                    "T cell", "T cell", "T cell", "Mono/DC", 
                    "Progenitor cells", "Progenitor cells", "Progenitor cells", "Progenitor cells",
                    "Progenitor cells", "Progenitor cells", "B cell", "NK", 
                    "Progenitor cells", "Mono/DC", "B cell", "Progenitor cells", "B cell")
  )
  bridgeTypeDict
  
  saveRDS(refTypeDict, "data/NeurIPS2021/refTypeDict.rds")
  saveRDS(bridgeTypeDict, "data/NeurIPS2021/bridgeTypeDict.rds")
}
refTypeDict <- readRDS(paste0(dat_dir, "data/NeurIPS2021/refTypeDict.rds"))
bridgeTypeDict <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeTypeDict.rds"))
AzimuthTypeDict <- readRDS(paste0(dat_dir, "data/NeurIPS2021/AzimuthTypeDict.rds"))

## Load data needed for annotation
YtrainName <- "celltype.l2"
reference <- readRDS(paste0(dat_dir, "data/NeurIPS2021/ref_pb.rds"))
refLabs <- colData(reference)[, YtrainName]
bridgeIntAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/bridgeIntAnn_cv.rds"))
NPintsAnn_cv <- readRDS(paste0(dat_dir, "output/Case2/NPintsAnn_cv.rds"))

## Some options
# Use normalised scores (NPintsScore_norm) or not 
useNormScore <- T
# Which is the ground truth
AzimuthAsTruth <- T
if(AzimuthAsTruth){
  originLabName <- "Azimuth.l2"
} else {
  originLabName <- "cellType"
}

## Loop for bridges
JJ <- length(bridgeIntAnn_cv)
errs_l <- vector("list", JJ)
names(errs_l) <- names(NPintsAnn_cv)
for(jj in 1:JJ){
  
  BridgeBatchName <- names(NPintsAnn_cv)[jj]
  NPintsAnn <- NPintsAnn_cv[[BridgeBatchName]]
  bridgeIntAnn <- bridgeIntAnn_cv[[BridgeBatchName]]
  queryBatchNames <- names(NPintsAnn)
  
  # Define output structure
  errs_names <- list(
    names(NPintsAnn), 
    c("BridgeIntOver", "BridgeIntBal", 
      "NPintsOver", "NPintsBal",
      "NP1nnOver", "NP1nnBal"
      # "NPldaOver", "NPldaBal"
    )
  )
  errs <- matrix(
    NA, 
    length(errs_names[[1]]), 
    length(errs_names[[2]])) %>%
    `dimnames<-`(errs_names)
  
  ## Loop for query datasets
  for(i in 1:length(NPintsAnn)){
    
    cat(jj, "th bridge;", i, "th query", "\n")
    
    # Get predicted classes
    refScoreName <- ifelse(useNormScore,
                           "YrefHat_norm",
                           "YrefHat")
    queryScoreName <- ifelse(useNormScore,
                             "queryScore_norm",
                             "queryScore")
    batchName <- queryBatchNames[i]
    bridgeIntLabs <- bridgeIntAnn[[batchName]]$predicted.l2
    originLabs <- bridgeATAC_list[[batchName]]@colData[,originLabName]
    NPintsLabs <- getClass(NPintsAnn[[batchName]][[queryScoreName]])
    class1nn <- as.character(
      FNN::knn(NPintsAnn[[batchName]][[refScoreName]], 
               NPintsAnn[[batchName]][[queryScoreName]], 
               k = 1, 
               cl = colData(reference)[, YtrainName]
      )
    )
    # mod_fit <- MASS::lda(NPintsAnn[[batchName]][[refScoreName]],
    #                      colData(reference)[,YtrainName])
    # classLDA <- as.character(predict(mod_fit, NPintsAnn[[batchName]][[queryScoreName]])$class)
    
    ## Translate fine labels to major cell type labels
    if(AzimuthAsTruth){
      originLabs_l1 <- translateLabel(originLabs, AzimuthTypeDict, "predicted.celltype.l2", "predicted.celltype.l1")
    } else {
      originLabs_l1 <- translateLabel(originLabs, refTypeDict, "original", "celltype.l1")
    }
    bridgeIntLabs_l1 <- translateLabel(bridgeIntLabs, refTypeDict, "celltype.l2", "celltype.l1")
    NPintsLabs_l1 <- translateLabel(NPintsLabs, refTypeDict, "celltype.l2", "celltype.l1")
    class1nn_l1 <- translateLabel(class1nn, refTypeDict, "celltype.l2", "celltype.l1")
    # classLDA_l1 <- translateLabel(classLDA, refTypeDict, "celltype.l2", "celltype.l1")
    
    errs[batchName,] <- c(
      classErr(bridgeIntLabs_l1, originLabs_l1)$err,
      classErr(NPintsLabs_l1, originLabs_l1)$err,
      classErr(class1nn_l1, originLabs_l1)$err
      # classErr(classLDA_l1, originLabs_l1)$err
    )
    
  }
  
  errs_l[[BridgeBatchName]] <- errs
}

if(F) saveRDS(errs_l, paste0(dat_dir, "output/Case2/bmmc_errs_list_", Sys.Date(), ".rds"))



### Viz -----------------------------------------------------------------------
errs_l_date <- "2023-10-12"
errs_l <- readRDS(paste0(dat_dir, "output/Case2/bmmc_errs_list_", errs_l_date, ".rds"))
## Overall: one box plot for each method
temp <- lapply(
  1:length(errs_l),
  function(x){
    BridgeBatchName <- names(errs_l)[x]
    out <- errs_l[[BridgeBatchName]]
    out <- out %>% 
      data.frame() %>%
      mutate(queryBatch = rownames(out),
             bridgeBatch = rep(BridgeBatchName, nrow(out))) 
    rownames(out) <- NULL
    return(out)
  }
) 
temp <- do.call(rbind, temp) %>%
  pivot_longer(-c(queryBatch, bridgeBatch), names_to = "method", values_to = "err") 
lvls <- c("BridgeIntOver", 
          "NPintsOver", 
          "NP1nnOver", 
          # "NP3nnOver", 
          # "NPldaOver",
          "BridgeIntBal",
          "NPintsBal",
          "NP1nnBal"
          # "NP3nnBal", 
          # "NPldaBal"
)
p_box <-
  temp %>%
  mutate(method = factor(method, levels = lvls)) %>%
  filter(method %in% c("NPintsBal", "BridgeIntBal")) %>%
  ggplot(aes(y = err, x = method, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#377EB8", "#E41A1C"), labels = c("Bridge", "PhiSpace")) +
  ylab("Balanced classification error") + xlab("Methods") +
  scale_x_discrete(labels= c("Bridge", "PhiSpace")) 
p_box

## By bridge batch
batchSizes <- readRDS(paste0(dat_dir, "data/NeurIPS2021/batchSizes_without_ILC.rds"))
outPlots <- vector("list", length(errs_l))
names(outPlots) <- names(errs_l)
for(BridgeBatchName in names(errs_l)){
  
  plot_dat <- errs_l[[BridgeBatchName]]
  
  # order batches by their sizes for plotting 
  orderBatch <- names(sort(batchSizes[rownames(plot_dat)])) 
  
  plot_dat <-
    plot_dat[orderBatch, ] %>%
    t() %>%
    as.data.frame() %>%
    mutate(method = rownames(.))
  p <-
    plot_dat %>%
    filter(method %in% c("NPintsBal", "BridgeIntBal")) %>%
    GGally::ggparcoord(columns = 1:(ncol(plot_dat)-1), 
                       groupColumn = "method", 
                       scale = "globalminmax",
                       showPoints = T) +
    scale_x_discrete(name = "", guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "", limits = c(0,0.5)) + 
    scale_colour_manual(values = c("#377EB8", "#E41A1C"), labels = c("Bridge", "PhiSpace")) +
    ggtitle(
      paste0(
        "Bridge = ", 
        BridgeBatchName, 
        " (#cells = ", 
        batchSizes[[BridgeBatchName]], 
        ")")
    ) 
  outPlots[[BridgeBatchName]] <- p
}
ggarrange(plotlist = outPlots,
          nrow = 4, ncol = 4, common.legend = T) 







## Old: analysis of reults for one bridge
if(F){
  bridgeIntAnn <- readRDS("output/Case2/bridgeIntAnn_s4d1_as_bridge.rds")
  NPintsAnn <- readRDS("output/Case2/NPintsAnn_s4d1_as_bridge.rds")
  errs <- matrix(NA, length(NPintsAnn), 10)
  dimnames(errs) <- list(
    names(NPintsAnn), 
    c("BridgeIntOver", "BridgeIntBal", "NPintsOver", "NPintsBal",
      "NP1nnOver", "NP1nnBal", "NP3nnOver", "NP3nnBal",
      "NPldaOver", "NPldaBal")
  )
  queryBatchNames <- names(NPintsAnn)
  for(i in 1:length(NPintsAnn)){
    
    batchName <- queryBatchNames[i]
    bridgeIntLabs <- bridgeIntAnn[[batchName]]$predicted.l2
    originLabs <- bridgeIntAnn[[batchName]]$cellType
    NPintsLabs <- getClass(NPintsAnn[[batchName]]$queryScore_norm)
    class1nn <- as.character(
      FNN::knn(NPintsAnn[[batchName]]$YrefHat_norm, 
               NPintsAnn[[batchName]]$queryScore_norm, 
               k = 1, 
               cl = colData(reference)[, YtrainName]
      )
    )
    class3nn <- as.character(
      FNN::knn(NPintsAnn[[batchName]]$YrefHat_norm, 
               NPintsAnn[[batchName]]$queryScore_norm, 
               k = 3, 
               cl = colData(reference)[, YtrainName]
      )
    )
    lda_fit <- MASS::lda(NPintsAnn[[batchName]]$YrefHat_norm, 
                         colData(reference)[, YtrainName])
    classLDA <- as.character(predict(lda_fit, NPintsAnn[[batchName]]$queryScore_norm)$class)
    
    originLabs_l1 <- translateLabel(originLabs, bridgeTypeDict, "original", "celltype.l1")
    bridgeIntLabs_l1 <- translateLabel(bridgeIntLabs, refTypeDict, "celltype.l2", "celltype.l1")
    NPintsLabs_l1 <- translateLabel(NPintsLabs, refTypeDict, "celltype.l2", "celltype.l1")
    class1nn_l1 <- translateLabel(class1nn, refTypeDict, "celltype.l2", "celltype.l1")
    class3nn_l1 <- translateLabel(class3nn, refTypeDict, "celltype.l2", "celltype.l1")
    classLDA_l1 <- translateLabel(classLDA, refTypeDict, "celltype.l2", "celltype.l1")
    
    errs[batchName,] <- c(
      classErr(bridgeIntLabs_l1, originLabs_l1)$err,
      classErr(NPintsLabs_l1, originLabs_l1)$err,
      classErr(class1nn_l1, originLabs_l1)$err,
      classErr(class3nn_l1, originLabs_l1)$err,
      classErr(classLDA_l1, originLabs_l1)$err
    )
    
    
    if(F){
      plotSankey3(bridgeIntLabs, originLabs, NPintsLabs)
      plotSankey3(bridgeIntLabs_l1, originLabs, NPintsLabs_l1)
      plotSankey3(bridgeIntLabs_l1, originLabs_l1, NPintsLabs_l1)
    }
  }
  plot_dat <-
    errs %>%
    t() %>%
    as.data.frame() %>%
    mutate(method = rownames(.))
  p1 <-
    plot_dat %>%
    filter(grepl("Bal", method)) %>%
    GGally::ggparcoord(columns = 1:(ncol(plot_dat)-1), 
                       groupColumn = "method", 
                       scale = "globalminmax",
                       showPoints = T) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "")
  
  errs <- errs_l[[11]]
  plot_dat <-
    errs %>%
    t() %>%
    as.data.frame() %>%
    mutate(method = rownames(.))
  p2 <-
    plot_dat %>%
    filter(grepl("Bal", method)) %>%
    GGally::ggparcoord(columns = 1:(ncol(plot_dat)-1), 
                       groupColumn = "method", 
                       scale = "globalminmax",
                       showPoints = T) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "")
  
  ggarrange(p1, p2, nrow = 2, ncol = 1, common.legend = T)
}