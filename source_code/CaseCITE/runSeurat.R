## Run Seurat to replicate PhiSpace results
.libPaths("/data/gpfs/projects/punim0613/JiaDong/Rlibs")

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

suppressPackageStartupMessages(library(Seurat))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"


## RNA
reference <- readRDS(
  paste0(
    dat_dir,
    "data/Covid/Haniffa/HaniffaRNA.rds"
  )
)
query <- readRDS(
  paste0(
    dat_dir,
    "data/Covid/Vento-Tormo/RNA.rds"
  )
)

commonFeat <- intersect(
  rownames(reference),
  rownames(query)
)
reference <- reference[commonFeat,]
query <- query[commonFeat,]


# Convert to Seurat objects
reference <- CreateSeuratObject(
  counts = assay(
    reference, "counts"
  ), 
  data = assay(
    reference, "logcounts"
  ),
  meta.data = colData(reference) %>% as.data.frame()
)
query <- CreateSeuratObject(
  counts = assay(
    query, "counts"
  ), 
  data = assay(
    query, "logcounts"
  ),
  meta.data = colData(query) %>% as.data.frame()
)



## Seurat V3 reference mapping 
SeuResPath <- paste0(
  dat_dir,
  "output/Case_CITE-seq/Covid/doublePhenoRNASeuRes.rds"
)
SeuResPathRef <- paste0(
  dat_dir,
  "output/Case_CITE-seq/Covid/doublePhenoRNASeuResRef.rds"
)

# Prepare ref
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, npcs = 60) # 60 cell types, matching PhiSpace
reference <- FindNeighbors(reference, dims = 1:60)
reference <- FindClusters(reference)
# Prepare query
queryAnchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = 1:60,
  reference.reduction = "pca"
)
refAnchors <- FindTransferAnchors(
  reference = reference,
  query = reference,
  dims = 1:60,
  reference.reduction = "pca"
)
# Prediction
ctype_pred <- TransferData(
  anchorset = queryAnchors, 
  refdata = reference$full_clustering,
  dims = 1:60
)
disease_pred <- TransferData(
  anchorset = queryAnchors, 
  refdata = reference$Status_on_day_collection_summary,
  dims = 1:60
)
ctype_pred_ref <- TransferData(
  anchorset = refAnchors, 
  refdata = reference$full_clustering,
  dims = 1:60
)
disease_pred_ref <- TransferData(
  anchorset = refAnchors, 
  refdata = reference$Status_on_day_collection_summary,
  dims = 1:60
)

saveRDS(
  list(
    ctype_pred,
    disease_pred
  ), 
  SeuResPath
)
saveRDS(
  list(
    ctype_pred_ref,
    disease_pred_ref
  ), 
  SeuResPathRef
)


## ADT
reference <- readRDS(
  paste0(
    dat_dir,
    "data/Covid/Haniffa/HaniffaADT.rds"
  )
)
query <- readRDS(
  paste0(
    dat_dir,
    "data/Covid/Vento-Tormo/ADT.rds"
  )
)

commonFeat <- intersect(
  rownames(reference),
  rownames(query)
)
reference <- reference[commonFeat,]
query <- query[commonFeat,]

# Convert to Seurat objects
reference <- CreateSeuratObject(
  counts = assay(
    reference, "counts"
  ), 
  data = assay(
    reference, "data"
  ),
  meta.data = colData(reference) %>% as.data.frame()
)
query <- CreateSeuratObject(
  counts = assay(
    query, "counts"
  ), 
  data = assay(
    query, "data"
  ),
  meta.data = colData(query) %>% as.data.frame()
)

## Seurat V3 reference mapping 
SeuResPath <- paste0(
  dat_dir,
  "output/Case_CITE-seq/Covid/doublePhenoADTSeuRes.rds"
)
SeuResPathRef <- paste0(
  dat_dir,
  "output/Case_CITE-seq/Covid/doublePhenoADTSeuResRef.rds"
)
# Prepare ref
reference <- FindVariableFeatures(reference, nfeatures = nrow(reference))
reference <- ScaleData(reference)
reference <- RunPCA(reference, npcs = 27) # 60 cell types, matching PhiSpace
reference <- FindNeighbors(reference, dims = 1:27)
reference <- FindClusters(reference)
# Prepare query
queryAnchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  dims = 1:27,
  reference.reduction = "pca"
)
refAnchors <- FindTransferAnchors(
  reference = reference,
  query = reference,
  dims = 1:27,
  reference.reduction = "pca"
)
# Prediction
ctype_pred <- TransferData(
  anchorset = queryAnchors, 
  refdata = reference$initial_clustering,
  dims = 1:27
)
disease_pred <- TransferData(
  anchorset = queryAnchors, 
  refdata = reference$Status_on_day_collection_summary,
  dims = 1:27
)
ctype_pred_ref <- TransferData(
  anchorset = refAnchors, 
  refdata = reference$initial_clustering,
  dims = 1:27
)
disease_pred_ref <- TransferData(
  anchorset = refAnchors, 
  refdata = reference$Status_on_day_collection_summary,
  dims = 1:27
)

saveRDS(
  list(
    ctype_pred,
    disease_pred
  ), 
  SeuResPath
)
saveRDS(
  list(
    ctype_pred_ref,
    disease_pred_ref
  ), 
  SeuResPathRef
)























