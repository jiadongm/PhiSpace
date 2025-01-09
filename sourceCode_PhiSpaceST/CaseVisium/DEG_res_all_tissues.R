# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

suppressPackageStartupMessages(library(Seurat))

source("~/PhiSpaceR/Case_CosMx/utils.R")
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

PhiAssay <- "log1p"
Mk2Use <- c("exp", "noexp")[1]

tissueNames <- tissueNames_Visium



# Azimuth markers (ground truth)
AzimuthMarkers <- readRDS(
  paste0(dat_dir, "data/LungRef/AzimuthLungMarkers.rds")
)
AzimuthMkExpanded <- qread(
  paste0(
    dat_dir, "data/LungRef/AzimuthLungMarkersExpanded.qs"
  )
)
if(!file.exists(AzExpandPath)){

  geneCor <- sparse.cor(
    t(assay(reference[selectedFeat,], "log1p"))
  )
  # rank each column
  geneCorRked <- selectFeat(geneCor, absVal = F)$orderedFeatMat
  AzimuthMkExpanded <- lapply(
    AzimuthMarkers,
    function(x){

      idx <- intersect(x, colnames(geneCorRked))
      unique(as.character(t(geneCorRked[, idx]) ))[1:100]
    }
  )
  qsave(AzimuthMkExpanded, AzExpandPath)
} else {

  AzimuthMkExpanded <- qread(AzExpandPath)
}

# AzimuthMkExpanded <- qread(
#   paste0(
#     dat_dir, "output/Case3/DEG_ref_sets.qs"
#   )
# )


methodNames <- c("PhiSpace",  "RCTD", "cell2location", "TACCO")
out <- matrix(NA, length(tissueNames), length(methodNames)) %>%
  `colnames<-`(methodNames) 
outList <- vector("list", length(tissueNames))
# TissueNames defined in utils.R
for(ii in 1:length(tissueNames)){
  
  
  VisTissueName <- tissueNames[ii]
  

  # Load DEG results
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes.qs"
  )
  outList_Phi <- qread(degResPath)
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_TACCO.qs"
  )
  outList_TACCO <- qread(degResPath)
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_RCTD.qs"
  )
  outList_RCTD <- qread(degResPath)
  degResPath <- paste0(
    dat_dir,
    "output/Case3/DEG/Visium", VisTissueName, "_degRes_c2l.qs"
  )
  outList_c2l <- qread(degResPath)
  
  
  # AzimuthMk2Use <- AzimuthMkExpanded
  # cTypes <- names(outList_Phi)
  AzimuthMk2Use <- AzimuthMarkers
  cTypes <- intersect(
    names(outList_Phi), names(AzimuthMk2Use)
  )
  
  mks <- matrix(NA, length(cTypes), 4) %>%
    `rownames<-`(cTypes) %>%
    `colnames<-`(methodNames)
  
  for(i in 1:length(cTypes)){
    
    cType <- cTypes[i]
    mks[cType,] <- c(
      PhiSpace = tempVenn(
        outList_Phi[[cType]],
        AzimuthMarkers = AzimuthMk2Use,
        plotVenn = F
      )$intersection,
      RCTD = tempVenn(
        outList_RCTD[[cType]],
        AzimuthMarkers = AzimuthMk2Use,
        plotVenn = F
      )$intersection,
      cell2location = tempVenn(
        outList_c2l[[cType]], 
        AzimuthMarkers = AzimuthMk2Use,
        plotVenn = F
      )$intersection,
      TACCO = tempVenn(
        outList_TACCO[[cType]], 
        AzimuthMarkers = AzimuthMk2Use,
        plotVenn = F
      )$intersection 
    )
  }
  
  # Fit ZINB
  suppressWarnings(fitRes <- (apply(mks, 2, fitZINB)))
  outList[[ii]] <- sapply(
    fitRes,
    function(x){
      unlist(x)
    }
  ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      tissueNames = rep(VisTissueName, nrow(.)),
      method = rownames(.)
    )
  
  out[ii,] <- mks %>% colMeans()
}

saveRDS(
  list(out, outList), 
  paste0(dat_dir, "output/Case3/DEG_summ", Mk2Use, ".rds")
)


  































