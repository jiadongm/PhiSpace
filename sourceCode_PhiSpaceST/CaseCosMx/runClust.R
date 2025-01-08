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

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
PhiAssay <- "log1p"
source("~/PhiSpaceR/Case_CosMx/utils.R")


# Query data
tissueNames <- c(
  "Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3",
  "Lung6",
  "Lung9_Rep1", "Lung9_Rep2",
  "Lung12",     "Lung13"
)

PhiResPath <- paste0(
  dat_dir,
  "output/Case3/CosMx/CosMxAllLungsPhiRes.qs"
)
PhiRes_list <- qread(PhiResPath)

tissueNames <- tissueNames_CosMx
query_list <- vector("list", length(tissueNames)) %>%
  `names<-`(tissueNames)
for(ii in 1:length(tissueNames)){
  tissueName <- tissueNames[ii]
  qs_dir <- paste0(dat_dir, "data/CosMx_lung/SCE_obj/", tissueName, "_SCE.qs")
  query <- qread(qs_dir)
  query <- zeroFeatQC(query)
  query <- logTransf(
    query,
    use_log1p = TRUE,
    targetAssay = "log1p"
  )
  reducedDim(query, "PhiSpace") <- PhiRes_list$PhiSpaceScore[[ii]]
  query_list[[ii]] <- query
}

# Clustering
for(ii in 1:length(tissueNames)){
  
  tissueName <- tissueNames[ii]
  
  query <- query_list[[tissueName]]
  toClust <- reducedDims(query)[["PhiSpace"]]
  nclust <- length(unique(query$cell_type))
  clust_Phi <- kmeans(toClust, centers = nclust, iter.max = 200, nstart = 20)
  toClust_cens <- apply(toClust, 2, censor, quant = 0.5)
  clust_PhiCens <- kmeans(toClust_cens, centers = nclust, iter.max = 200, nstart = 20)
  qsave(
    clust_Phi,
    paste0(
      dat_dir, "output/Case3/CosMx/CosMx_", tissueName, "_clust_Phi.qs"
    )
  )
  qsave(
    clust_PhiCens,
    paste0(
      dat_dir, "output/Case3/CosMx/CosMx_", tissueName, "_clust_PhiCens.qs"
    )
  )
  
  toClust <- t(
    assay(query, "log1p")
  )
  pca_res <- getPC(toClust, ncomp = ncol(toClust_cens))
  clust_log1p <- kmeans(pca_res$scores, centers = nclust, iter.max = 200, nstart = 20)
  qsave(
    clust_log1p,
    paste0(
      dat_dir, "output/Case3/CosMx/CosMx_", tissueName, "_clust_log1p.qs"
    )
  )
}






































