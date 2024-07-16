

### Update: 26 April 2024. Realised that previously, multiome object did not contain RNA counts
# So this means we probabaly used integrated and normalised RNA data

##1 Regenerate RNA counts object --------------------------------------------------
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad(
  paste0(
    dat_dir,
    "data/NeurIPS2021/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
  )
)

gc()

countMat <- as(t(adata$layers["counts"]), 'CsparseMatrix') 
dimnames(countMat) <- list(
  adata$var_names$values,
  adata$obs_names$values
)
ann <- adata$obs

rm(adata); gc()

multiome <- SeuratObject::CreateSeuratObject(
  SeuratObject::CreateAssayObject(
    countMat
  )
)

multiome@meta.data <- as.data.frame(ann)

saveRDS(
  multiome,
  paste0(
    dat_dir,
    "data/NeurIPS2021/multiome_cor.rds"
  )
)

##2 BridgeRNA and bridgeATAC separation
# Need feature_types to extract GEX, ATAC_gene_names
feature_types <- adata$var$feature_types
saveRDS(
  feature_types, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/feature_types_cor.rds"
  )
)
ATAC_gene_names <- adata$uns$ATAC_gene_activity_var_names

# ATAC assay of bridge
bridge_ATAC <- as(
  t(adata$obsm["ATAC_gene_activity"]),
  'CsparseMatrix'
)
rownames(bridge_ATAC) <- ATAC_gene_names

# RNA assay of bridge
bridge_RNA <- multiome[["RNA"]]$counts[feature_types == "GEX", ]



bridgeRNA <- SingleCellExperiment(
  assays = list(counts = bridge_RNA),
  colData = multiome@meta.data[,!grepl("ATAC", colnames(multiome@meta.data))]
)
bridgeATAC <- SingleCellExperiment(
  assays = list(counts = bridge_ATAC),
  colData = multiome@meta.data[,grepl("ATAC", colnames(multiome@meta.data))]
)

YtrainName <- "cell_type"
colData(bridgeATAC)[,YtrainName] <- colData(bridgeRNA)[,YtrainName]

saveRDS(
  bridgeRNA, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/bridgeRNA_cor.rds"
  )
)
saveRDS(
  bridgeATAC, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/bridgeATAC_cor.rds"
  )
)


rm(list = ls()); gc()



##3 Generate obj.multi.batches object -----------------------------------------
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
bridgeRNA <- readRDS(
  paste0(
    dat_dir,
    "data/NeurIPS2021/bridgeRNA_cor.rds"
  )
)
feature_types <- readRDS(
  paste0(
    dat_dir,
    "data/NeurIPS2021/feature_types_cor.rds"
  )
)
obj.multi.temp <- readRDS(
  paste0(
    dat_dir,
    "data/NeurIPS2021/multiome_cor.rds"
  )
)
orig.meta.data <- obj.multi.temp@meta.data
all(colnames(bridgeRNA) == colnames(obj.multi.temp@assays$RNA$counts))

## Add RNA assay
obj.multi0 <- CreateSeuratObject(counts = assay(bridgeRNA, "counts"))

## Add the ATAC-seq assay
atac_counts <- obj.multi.temp@assays$RNA$counts[feature_types == "ATAC",]
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# Get gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Change style to UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# Add in ATAC-seq data as ChromatinAssay object
obj.multi0[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  # sep = c(":", "-"),
  genome = 'hg38',
  # fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
rm(obj.multi.temp); gc()
saveRDS(
  orig.meta.data, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/orig.meta.data_cor.rds"
  )
)
saveRDS(
  obj.multi0, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/obj.multi_cor.rds"
  )
)

### Split batches

## Batch indices
all(rownames(obj.multi0@meta.data) == rownames(orig.meta.data))
obj.multi0$batch <- orig.meta.data$batch
obj.multi.batches <- SplitObject(obj.multi0, split.by = "batch")
rm(obj.multi0); gc()


# Normalise
obj.multi.batches <- lapply(X = obj.multi.batches, FUN = NormalizeData)

# Original annotations
for (i in 1:length(obj.multi.batches)) {
  
  batchName <- names(obj.multi.batches)[i]
  temp <- orig.meta.data %>%
    filter(batch == batchName)
  if(all(rownames(obj.multi.batches[[i]]@meta.data) == rownames(temp))){
    
    obj.multi.batches[[i]]@meta.data[,"cellType"] <- temp$cell_type
    
  }
  
}
saveRDS(
  obj.multi.batches, 
  paste0(
    dat_dir,
    "data/NeurIPS2021/obj.multi.batches_cor.rds"
  )
)


####################################################################################
## Two versions of ATAC: gene activity and TF-IDF normalised -------------------------
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))

dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"

# Gene activity
obj.multi.batches <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches_cor.rds"))
bridgeATAC <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_cor.rds"))

# Create sce objects
bridgeATAC_list <- vector("list", length(obj.multi.batches))
batchNames <- names(obj.multi.batches)
names(bridgeATAC_list) <- batchNames
for(i in 1:length(obj.multi.batches)){
  
  batchName <- batchNames[i]
  batchColData <- obj.multi.batches[[batchName]]@meta.data
  bridgeATAC_list[[batchName]] <- SingleCellExperiment(
    assays = list(counts = assay(bridgeATAC, "counts")[,rownames(batchColData)]),
    colData = batchColData
  )
}
saveRDS(bridgeATAC_list, paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_list_cor.rds"))
# scran normalise
library(scran)
bridgeATAC_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_list_cor.rds"))
for(i in 1:length(bridgeATAC_list)){
  
  cat("Scran normalise", i, "th\n")
  batchName <- names(bridgeATAC_list)[[i]]
  query <- bridgeATAC_list[[batchName]]
  query <- computeSumFactors(query, cluster = quickCluster(query))
  query <- logNormCounts(query)
  bridgeATAC_list[[batchName]] <- query
}
saveRDS(bridgeATAC_list, paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_list_cor.rds"))


# TF-IDF
all(
  startsWith(
    rownames(obj.multi.batches[[1]]@assays$ATAC$counts),
    "chr"
  )
)

# Create sce objects bridgeATAC_list
batchNames <- names(obj.multi.batches)
bridgeATAC_list <- vector("list", length(obj.multi.batches)) %>%
  setNames(batchNames)
for(i in 1:length(obj.multi.batches)){
  
  batchName <- batchNames[i]
  batchObj <- obj.multi.batches[[batchName]]
  batchObj <- Signac::RunTFIDF(batchObj)
  batchColData <- batchObj@meta.data
  bridgeATAC_list[[batchName]] <- SingleCellExperiment(
    assays = list(data = batchObj[["ATAC"]]$data),
    colData = batchColData
  )
}
saveRDS(bridgeATAC_list, paste0(dat_dir, "data/NeurIPS2021/bridgeATAC_TFIDF_list_cor.rds"))




### reference dataset -----------------------------------------------------------
obj.rna <- readRDS(
  paste0(
    dat_dir,
    "data/NeurIPS2021/obj.rna_for_refMap.rds"
  )
)
reference <- SingleCellExperiment(
  list(counts = obj.rna@assays$RNA$counts),
  colData = obj.rna@meta.data
)
refTypeDict <- 
  reference@colData[,c("celltype.l1", "celltype.l2")] %>%
  as.data.frame() %>%
  distinct(celltype.l2, .keep_all = T)
refTypeDict
saveRDS(reference, paste0("data/NeurIPS2021/obj.rna_for_refMap_sce.rds"))


# Create bridgeRNA_list
obj.multi.batches <- readRDS(paste0(dat_dir, "data/NeurIPS2021/obj.multi.batches_cor.rds"))
bridgeRNA <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_cor.rds"))
bridgeRNA_list <- vector("list", length(obj.multi.batches))
batchNames <- names(obj.multi.batches)
names(bridgeRNA_list) <- batchNames
for(i in 1:length(obj.multi.batches)){
  
  batchName <- batchNames[i]
  batchColData <- obj.multi.batches[[batchName]]@meta.data
  bridgeRNA_list[[batchName]] <- SingleCellExperiment(
    assays = list(counts = assay(bridgeRNA, "counts")[,rownames(batchColData)]),
    colData = batchColData
  )
}
saveRDS(bridgeRNA_list, paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))
# scran normalise
bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))
library(scran)
for(i in 1:length(bridgeRNA_list)){
  
  cat("Scran normalise", i, "th\n")
  batchName <- names(bridgeRNA_list)[[i]]
  query <- bridgeRNA_list[[batchName]]
  query <- computeSumFactors(query, cluster = quickCluster(query))
  query <- logNormCounts(query)
  bridgeRNA_list[[batchName]] <- query
}
saveRDS(bridgeRNA_list, paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_cor.rds"))









#### Azimuth ref and bridge (they intergrated Neurips2021 as bridge)
library(Seurat)
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/"
ext <- readRDS(paste0(dat_dir, "data/10x_multiome/ext.rds"))
ext <- UpdateSeuratObject(ext)

ext@assays$refAssay$data

sce <- SingleCellExperiment::SingleCellExperiment(
  list(data = ext[["refAssay"]]$data),
  colData = ext@meta.data
)
saveRDS(sce, paste0(dat_dir, "data/10x_multiome/ext_sce.rds"))















