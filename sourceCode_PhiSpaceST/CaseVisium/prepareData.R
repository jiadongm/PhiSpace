# Load packages and predefined quantities
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs)) # quick read and write of R objects

### Process Azimuth Lung Reference
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
# Azimuth Lung v2
# https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
reference <- readRDS(paste0(dat_dir, "data/LungRef/AzimuthLung2.0.rds"))
rowdat <- reference@assays$RNA@meta.features

# Look at annotations
reference$cell_type %>% table() %>% length()
reference$ann_finest_level %>% table() %>% names()
reference$disease %>% table()
reference$ann_coarse_for_GWAS_and_modeling %>% table() %>% names()
plotSankey3(colData(reference)[,YtrainName] %>% as.character(),
            reference$ann_level_4 %>% as.character(),
            reference$ann_finest_level %>% as.character())


reference <- SingleCellExperiment(
  list(counts = reference@assays$RNA$counts),
  colData = reference@meta.data
)
all(rownames(reference) == rownames(rowdat))
rownames(reference) <- rowdat$feature_name
qsave(reference, paste0(dat_dir, "data/LungRef/AzimuthLung2.0_sce.qs"))

## Subsetting reference
# Reference
YtrainName <- refLabName <- "ann_finest_level"
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
reference <- qread(
  paste0(
    dat_dir,
    "data/LungRef/AzimuthLung2.0_sce.qs"
  )
)
ref_sub <- subsample(
  reference,
  key = "ann_finest_level",
  proportion = 0.1,
  minCellNum = 50
)
ref_sub <- zeroFeatQC(ref_sub)
ref_sub <- logTransf(
  ref_sub,
  use_log1p = TRUE,
  targetAssay = "log1p"
)
qsave(ref_sub, refPath)

# Compute feature importance scores for feature selection
refPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs"
)
reference <- qread(refPath)
impScPath <- paste0(
  dat_dir,
  "output/Case3/refImpScores.rds"
)
if(!file.exists(impScPath)){
  
  tuneRes <- tunePhiSpace(
    reference = reference,
    assayName = "log1p",
    phenotypes = YtrainName,
    tune_ncomp = F,
    tune_nfeat = F
  )
  impScores <- tuneRes$impScores
  saveRDS(impScores, impScPath)
} else {
  
  impScores <- readRDS(impScPath)
}


# Create h5ad objects for benchmarking methods in Python
refSelectPath <- paste0(
  dat_dir,
  "data/LungRef/AzimuthLung2.0_sce_0.1sub_selectFeat.h5ad"
)
if(!file.exists(refSelectPath)){
  
  ref_selected <- reference[selectedFeat, ]
  write.csv(selectedFeat, paste0(dat_dir, "output/Case3/LungRef0.1selectedFeat.csv"), row.names = F)
  
  library(reticulate)
  use_condaenv("TACCO_env")
  ad <- import("anndata")
  
  ref_ad <- ad$AnnData(
    X = t(
      assay(ref_selected, "counts")
    )
  )
  ref_ad$var_names <- rownames(ref_selected)
  ref_ad$obs_names <- colnames(ref_selected)
  ref_ad$obs <- colData(ref_selected) %>% as.data.frame()
  ref_ad$write_h5ad(
    refSelectPath
  )
}


###### List of query Visium objects
# The data can be downloaded from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530
query_list <- lapply(
  tissueNames_Visium,
  function(VisTissueName){
    
    visPath <- paste0(
      dat_dir,
      "data/Visium_NSCLC/",
      VisTissueName,
      "/sce.rds"
    )
    if(!file.exists(visPath)){
      
      # run this in command line: mkdir spatial && tar -xvf XXX-spatial.tar -C spatial 
      qu_vis <- Seurat::Load10X_Spatial(
        paste0(
          dat_dir,
          "data/Visium_NSCLC/",
          VisTissueName
        ),
        filename = paste0(
          VisTissueName,
          "-filtered_feature_bc_matrix.h5"
        )
      )
      spat_info <- data.frame(
        x = qu_vis@images$slice1$centroids@coords[,2],
        y = - qu_vis@images$slice1$centroids@coords[,1]
      )
      qu_vis <- Seurat::AddMetaData(
        qu_vis,
        spat_info
      )
      qu_vis <- SingleCellExperiment(
        list(
          counts = qu_vis[["Spatial"]]$counts
        ),
        colData = qu_vis@meta.data %>% as.data.frame()
      )
      # qu_vis <- scranTransf(
      #   qu_vis
      # )
      qu_vis <- logTransf(qu_vis)
      
      assay(qu_vis, "log1p") <- log1p(
        assay(qu_vis, "counts")
      )
      
      saveRDS(qu_vis, visPath)
      
      # source("~/PhiSpaceR/CaseCosMx/CaseCosMx_prepareH5AD.R")
    } else {
      
      qu_vis <- readRDS(visPath)
    }
    
    return(qu_vis)
  }
)
names(query_list) <- tissueNames_Visium