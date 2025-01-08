
# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

suppressPackageStartupMessages(library(spacexr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))

suppressPackageStartupMessages(library(parallel))

dat_dir <- "/home/unimelb.edu.au/jmao1/PhiSpace/"

# Reference -------------------------------------------------------
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
YtrainName <- refLabName <- "ann_finest_level"
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

  
  VisTissueName <- tissueNames[ii]
  
  # Query -----------------------------------------------------------
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
  
  # Results ------------------------------------------------------------
  
  # RCTD
  RCTDpath <- paste0(
    dat_dir,
    "output/Case3/RCTD/NSCLC_Visium", VisTissueName, "_RCTDres.rds"
  )
  myRCTD <- readRDS(RCTDpath)
  
  # PhiSpace
  PhiResPath <- paste0(
    dat_dir,
    "output/Case3/PhiSpace/combo_PhiRes.qs"
  )
  PhiRes <- qread(PhiResPath)
  PhiScores <- PhiRes$PhiSpaceScore[[VisTissueName]]
  PhiScores_norm <- normPhiScores(PhiScores)
  
  # TACCO
  tacco_res <- read.csv(
    paste0(
      dat_dir,
      "output/Case3/TACCO/Visium_NSCLC_", VisTissueName, "_TACCO.csv"
    ),
    row.names = 1
  )
  
  # cell2location
  q05path <- paste0(
    dat_dir,
    "output/Case3/c2l/cell2location_", VisTissueName, "/cell2location_map/q05_abundance.csv"
  )
  q05_res <- read.csv(q05path, row.names = 1)
  colnames(q05_res) <- colnames(PhiScores_norm)
  
  
  # RCTD did some filtering of the spots
  qu_vis <- qu_vis[,rownames(myRCTD@results$weights)]
  reducedDim(qu_vis, "RCTD") <- spacexr::normalize_weights(myRCTD@results$weights)
  reducedDim(qu_vis, "PhiSpace") <- PhiScores_norm[colnames(qu_vis),]
  reducedDim(qu_vis, "TACCO") <- tacco_res[colnames(qu_vis),] %>%
    `colnames<-`(colnames(reducedDim(qu_vis, "PhiSpace")))
  reducedDim(qu_vis,"c2l_q05") <- q05_res[colnames(qu_vis),]
  
  
  # Marker analysis -------------------------------------------------------------
  # cell2location
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_c2l.qs"
  )
  
  
  tempDEGfun <- function(cType){
    
    sc <- reducedDim(qu_vis, "c2l_q05")[,cType] %>%
      `names<-`(rownames(reducedDim(qu_vis, "c2l_q05")))
    
    # Choose a cutoff and select all spots with AT0 scores higher than cutoff.
    (cutoff <- quantile(reducedDim(qu_vis, "c2l_q05")[,cType], 0.8))
    idx <- names(sc)[sc > cutoff]
    groupIdx <- rep("other", ncol(qu_vis)) %>%
      `names<-`(colnames(qu_vis))
    groupIdx[idx] <- cType
    colData(qu_vis)[,"groupIdx"] <- groupIdx
    
    toKeep <- rep(T, nrow(qu_vis))
    
    seu <- CreateSeuratObject(
      counts = assay(qu_vis[toKeep,], "counts"),
      meta.data = as.data.frame(colData(qu_vis[toKeep,]))
    )
    seu[["RNA"]]$data <- assay(qu_vis[toKeep,], "log1p")
    Idents(seu) <- "groupIdx"
    
    deg.results <- FindAllMarkers(
      seu,
      only.pos = TRUE,
      logfc.threshold = 0.5,
      min.pct = 0.25,
      min.diff.pct = 0,
      test.use = "MAST",
      verbose = F
    )
    
    return(deg.results[deg.results$cluster==cType,])
  }
  
  tik <- Sys.time()
  outList <- mclapply(
    colnames(reducedDim(qu_vis, "c2l_q05")),
    tempDEGfun,
    mc.cores = 10
  )
  names(outList) <- colnames(reducedDim(qu_vis, "c2l_q05"))
  Sys.time() - tik
  
  qsave(outList, degResPath)
  
  
  # PhiSpace
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes.qs"
  )
  
  tempDEGfun <- function(cType){
    sc <- reducedDim(qu_vis, "PhiSpace")[,cType]
    
    # Choose a cutoff and select all spots with AT0 scores higher than cutoff.
    (cutoff <- quantile(reducedDim(qu_vis, "PhiSpace")[,cType], 0.8))
    idx <- names(sc)[sc > cutoff]
    groupIdx <- rep("other", ncol(qu_vis)) %>%
      `names<-`(colnames(qu_vis))
    groupIdx[idx] <- cType
    colData(qu_vis)[,"groupIdx"] <- groupIdx
    
    toKeep <- rep(T, nrow(qu_vis))
    
    seu <- CreateSeuratObject(
      counts = assay(qu_vis[toKeep,], "counts"),
      meta.data = as.data.frame(colData(qu_vis[toKeep,]))
    )
    seu[["RNA"]]$data <- assay(qu_vis[toKeep,], "log1p")
    Idents(seu) <- "groupIdx"
    
    deg.results <- FindAllMarkers(
      seu,
      only.pos = TRUE,
      logfc.threshold = 0.5,
      min.pct = 0.25,
      min.diff.pct = 0,
      test.use = "MAST",
      verbose = F
    )
    
    deg.results <- deg.results[deg.results$cluster==cType,]
    
    return(deg.results)
  }
  
  tik <- Sys.time()
  outList <- mclapply(
    colnames(reducedDim(qu_vis, "c2l_q05")),
    tempDEGfun,
    mc.cores = 10
  )
  names(outList) <- colnames(reducedDim(qu_vis, "c2l_q05"))
  Sys.time() - tik
  
  qsave(outList, degResPath)
  
  
  # TACCO
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_TACCO.qs"
  )
  tempDEGfun <- function(cType){
    
    sc <- reducedDim(qu_vis, "TACCO")[,cType] %>%
      `names<-`(rownames(reducedDim(qu_vis, "TACCO")))
    
    # Choose a cutoff and select all spots with AT0 scores higher than cutoff.
    (cutoff <- quantile(reducedDim(qu_vis, "TACCO")[,cType], 0.8))
    idx <- names(sc)[sc > cutoff]
    groupIdx <- rep("other", ncol(qu_vis)) %>%
      `names<-`(colnames(qu_vis))
    groupIdx[idx] <- cType
    colData(qu_vis)[,"groupIdx"] <- groupIdx
    
    toKeep <- rep(T, nrow(qu_vis))
    
    seu <- CreateSeuratObject(
      counts = assay(qu_vis[toKeep,], "counts"),
      meta.data = as.data.frame(colData(qu_vis[toKeep,]))
    )
    seu[["RNA"]]$data <- assay(qu_vis[toKeep,], "log1p")
    Idents(seu) <- "groupIdx"
    
    deg.results <- FindAllMarkers(
      seu,
      only.pos = TRUE,
      logfc.threshold = 0.5,
      min.pct = 0.25,
      min.diff.pct = 0,
      test.use = "MAST",
      verbose = F
    )
    
    deg.results <- deg.results[deg.results$cluster==cType,]
    
    return(deg.results)
  }
  
  tik <- Sys.time()
  outList <- mclapply(
    colnames(reducedDim(qu_vis, "c2l_q05")),
    tempDEGfun,
    mc.cores = 10
  )
  names(outList) <- colnames(reducedDim(qu_vis, "c2l_q05"))
  Sys.time() - tik
  
  qsave(outList, degResPath)
  
  
  # RCTD
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_RCTD.qs"
  )
  
  tempDEGfun <- function(cType){
    
    sc <- reducedDim(qu_vis, "RCTD")[,cType]
    
    # Choose a cutoff and select all spots with AT0 scores higher than cutoff.
    (cutoff <- quantile(reducedDim(qu_vis, "RCTD")[,cType], 0.8))
    idx <- names(sc)[sc > cutoff]
    groupIdx <- rep("other", ncol(qu_vis)) %>%
      `names<-`(colnames(qu_vis))
    groupIdx[idx] <- cType
    colData(qu_vis)[,"groupIdx"] <- groupIdx
    
    toKeep <- rep(T, nrow(qu_vis))
    
    seu <- CreateSeuratObject(
      counts = assay(qu_vis[toKeep,], "counts"),
      meta.data = as.data.frame(colData(qu_vis[toKeep,]))
    )
    seu[["RNA"]]$data <- assay(qu_vis[toKeep,], "log1p")
    Idents(seu) <- "groupIdx"
    
    deg.results <- FindAllMarkers(
      seu,
      only.pos = TRUE,
      logfc.threshold = 0.5,
      min.pct = 0.25,
      min.diff.pct = 0,
      test.use = "MAST",
      verbose = F
    )
    
    deg.results <- deg.results[deg.results$cluster==cType,]
    
    return(deg.results)
  }
  
  tik <- Sys.time()
  outList <- mclapply(
    colnames(reducedDim(qu_vis, "c2l_q05")),
    tempDEGfun,
    mc.cores = 10
  )
  names(outList) <- colnames(reducedDim(qu_vis, "c2l_q05"))
  Sys.time() - tik
  
  qsave(outList, degResPath)
}
  


