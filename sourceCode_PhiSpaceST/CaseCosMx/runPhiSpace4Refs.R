# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/PhiSpace-ST_submit/" # replace with your own directory
source(paste0(dat_dir, "CaseCosMx/utils.R"))

PhiAssay <- "log1p"

# Load reference
YtrainName <- c("manual_annotation_1")
ref_list <- lapply(
  c("immune", "epithelial", "endothelial", "mesenchymal"),
  function(lineage){
    
    qread(paste0(dat_dir, "data/LungFibrosis/scRNA-seq/", lineage, "_sce.qs"))
  }
)
names(ref_list) <- c("immune", "epithelial", "endothelial", "mesenchymal")

# Process query datasets
tissueNames <- tissueNames_CosMx
query_list <- vector("list", length(tissueNames)) %>% `names<-`(tissueNames)
for(ii in 1:length(tissueNames)){
  tissueName <- tissueNames[ii]
  qs_dir <- paste0(dat_dir, "data/CosMx_lung/SCE_obj/", tissueName, "_SCE.qs")
  query <- qread(qs_dir)
  query <- zeroFeatQC(query)
  query <- logTransf(
    query,
    use_log1p = TRUE,
    targetAssay = PhiAssay
  )
  query_list[[ii]] <- query
}
qsave(query_list, paste0(dat_dir, "data/CosMx_lung/SCE_obj/allTissue_SCE.qs"))

# Annotation
sc_list <- vector("list", length(ref_list))
names(sc_list) <- names(ref_list)
for(ii in 1:length(ref_list)){
  
  lineage <- names(ref_list)[ii]
  reference <- ref_list[[lineage]]
  reference <- logTransf(reference, targetAssay = PhiAssay, use_log1p = T)
  
  PhiRes <- PhiSpaceR_1ref(
    reference, 
    query_list, 
    phenotypes = YtrainName, 
    PhiSpaceAssay = PhiAssay, 
    regMethod = "PLS", 
    center = T,
    scale = F
  )
  
  sc <- PhiRes$PhiSpaceScore
  sc <- lapply(
    sc,
    function(x){
      xnorm <- normPhiScores(x)
      colnames(xnorm) <- paste0(colnames(xnorm), "(", lineage, ")")
      return(xnorm)
    }
  )
  
  sc_list[[lineage]] <- sc
}

PhiResPath <- paste0(
  dat_dir, "output/Case3/CosMxAllLungsPhiRes4Refs_", PhiAssay, ".qs"
)
qsave(sc_list, PhiResPath)








