# See https://github.com/meiosis97/Sincast for a tutorial on Sincast

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

query <- readRDS(paste0(dat_dir, "output/Case1/query_Rosa_sub.rds"))
reference <- readRDS(paste0(dat_dir,"data/stemformatics/ref_dc.rds"))
query <- RankTransf(query, "counts")
reference <- RankTransf(reference, "data", sparse = F)

PhiResPath <- paste0(
  dat_dir,
  "output/Case1/PhiRes_RosaNoSubset.rds"
)
PhiRes <- readRDS(PhiResPath)
selectedFeat <- rownames(PhiRes$atlas_re$reg_re$loadings)
# inclusion <- as.logical(rowData(reference)$inclusion)
# selectedFeat <- rownames(reference)[inclusion]
# selectedFeat <- intersect(rownames(query), selectedFeat)

suppressPackageStartupMessages(library(Sincast))
SinResPath <- paste0(
  dat_dir,
  "output/Case1/SinImputRosa.rds"
)

querySin <- rcTransform(query)[selectedFeat, ]
refSin <- reference[selectedFeat, ]

rm(reference, query, PhiRes); gc()



SinRunTime <- Sys.time()

# refSin <- featureWeighting(refSin, clusterid = "Cell Type")
# c(refSin, querySin) %<-% filterData(refSin, querySin)
# refSin <- makeAtlas(reference = refSin, col.by = "Cell Type", vis.atlas = F)
querySin <- sincastImp(querySin, col.by = "mainTypes", saveGrpah = F)
querySin <- postScale(querySin)

(SinRunTime <- Sys.time() - SinRunTime)

idx <- colnames(readRDS(paste0(dat_dir, "output/Case1/query_Rosa_sub.rds")))
querySin <- querySin[,idx]
metadata(querySin) <- list()

saveRDS(querySin, SinResPath)
