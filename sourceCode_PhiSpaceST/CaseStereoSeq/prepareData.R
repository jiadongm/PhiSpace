suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(ggplot2))
source("/data/gpfs/projects/punim0613/JiaDong/PhiSpace/CaseStereoSeq/utils.R")
dat_dir <- "/data/projects/punim0613/JiaDong/Dawson/MouseSpleen/"

### Stereo-seq mouse spleen with AML barcodes from Holze et al. (2024)
library(reticulate)
library(Matrix)
sc <- import("scanpy")
binsize <- 50
ad <- sc$read_h5ad(
  paste0(
    dat_dir, 
    "stereo-seq_m4_additional_bin_sizes/mouse4_bin",
    binsize,
    "_bc.h5ad"
  )
)
cMat <- ad$X %>%
  `dimnames<-`(
    list(ad$obs_names$format(), ad$var_names$format())
  ) %>%
  t() %>% 
  as('CsparseMatrix')
colDat <- ad$obs %>%
  mutate(
    x = ad$obsm["spatial"][,1],
    y = -ad$obsm["spatial"][,2]
  ) # adjusted x y coordinates to make similar to paper
query <- SingleCellExperiment::SingleCellExperiment(
  list(counts = cMat), colData = colDat
)
VizSpatial(query, ptSize = 0.5)
qsave(
  query,
  paste0(
    dat_dir, "stereo-seq_m4_additional_bin_sizes/mouse4_bin",
    binsize, "_bc_sce.qs"
  )
)

### scRNA-seq from the same mouse as Stereo-seq above, als from Holze et al. (2024)
querySC <- readRDS(
  paste0(dat_dir, "MLL_M2_cc_regressed_normalised_singlets.rds")
)
querySC <- UpdateSeuratObject(querySC)
sce <- SingleCellExperiment(
  list(counts = querySC[["RNA"]]$counts),
  colData = querySC@meta.data
)
qsave(
  sce, paste0(dat_dir, "mouse4_scRNAseq_sce.qs")
)



### Mouse Spleen scRNA-seq from Zhang et al. (2023)
ref_seurat <- readRDS(paste0(dat_dir, "scRNA-seq/consistent_version_SPF.rds"))
reference <- SingleCellExperiment(
  list(
    counts = ref_seurat@assays$RNA@counts
  ),
  colData = ref_seurat@meta.data
)
# scran normalisation
clusts <- quickCluster(reference, assay.type = "counts")
reference <- computeSumFactors(reference, assay.type = "counts", cluster = clusts)
if(min(reference$sizeFactor) <= 0){
  smallestPositive <- reference$sizeFactor[which(sort(reference$sizeFactor) > 0)[1]]
  reference$sizeFactor[reference$sizeFactor <= 0] <- smallestPositive
}
reference <- logNormCounts(reference, assay.type = "counts")
# log1p normalisation
reference <- logTransf(reference, use_log1p = TRUE)
saveRDS(reference, paste0(dat_dir, "scRNA-seq/consistent_version_SPF_unint_sce.rds"))





### Mouse bone marrow scRNA-seq from Harris et al. (2021)
library(Seurat)
library(SeuratDisk)
library(Matrix)
# Original file is a loom obj
loom.fn <- paste0(dat_dir, "scRNA-seq/Harris/processed_droplet_data.loom")
s_cnct <- Connect(filename = loom.fn, mode = "r")
s_cnct |> names()
s_cnct[["matrix"]]
s_cnct[["matrix"]]$dims
ncells <- s_cnct[["matrix"]]$dims[1]
ngenes <- s_cnct[["matrix"]]$dims[2]
s_cnct[["matrix"]][1:10,1:10]
mat <- as.sparse(
  s_cnct[["matrix"]][1:ncells,1:ngenes]
)
mat <- t(mat)

# Gene list
gene_names <- s_cnct[["row_attrs"]][["var_names"]][1:ngenes]
# Cell ID
cell_names <- s_cnct[["col_attrs"]][["obs_names"]][1:ncells]
dimnames(mat) <- list(
  gene_names,
  cell_names
)
# Cell type
s_cnct[["col_attrs"]] |> names()
s_cnct[["col_attrs"]][["scNym_str"]][1:ncells] |> table()
cell_types <- s_cnct[["col_attrs"]][["scNym"]][1:ncells] 
# SCE obj
ref <- SingleCellExperiment::SingleCellExperiment(
  list(data = mat),
  colData = data.frame(
    cellType = cell_types
  )
)
saveRDS(ref, paste0(dat_dir, "scRNA-seq/Harris/processed_droplet_data_sce.rds"))
pathRefBM <- paste0(
  dat_dir, 
  "scRNA-seq/refBM_processed.qs"
)
# Subsampe refBM
subsample <- function(
    reference, 
    key = "cellType",
    maxSize = 5000,
    seed = 510
){
  
  set.seed(seed)
  
  labs <- colData(reference)[,key]
  freqTable <- table(labs)
  labs2subsamp <- names(freqTable[freqTable > maxSize])
  
  if(length(labs2subsamp) >= 1){
    
    idx_mat <- matrix(NA, nrow = length(labs2subsamp), ncol = length(labs))
    for(i in 1:length(labs2subsamp)){
      
      lab <- labs2subsamp[i]
      idx <- (labs == lab)
      
      subIdx <- sample(1:sum(idx), maxSize)
      idx[which(idx)[subIdx]] <- FALSE
      idx_mat[i, ] <- !idx
    }
    
    out <- apply(idx_mat, 2, all)
  } else {
    
    out <- rep(TRUE, length(labs))
  }
  
  
  return(reference[,out])
}
refBM <- subsample(refBM, maxSize = 5000)
table(refBM$cellType)
dim(refBM)
refBM <- RankTransf(refBM, "data")
qsave(refBM, pathRefBM)





### Mouse spleen CITE-seq data from Gayoso et al. (2021)
library(reticulate)
library(Matrix)
library(qs)
sc <- import("scanpy")
pd <- import("pandas")
# Two batches
ref1 <- sc$read_h5ad("/data/gpfs/projects/punim0613/JiaDong/PhiSpace/data/Stereo-CITE/CITE-seq/spleen_lymph_206.h5ad")
ref2 <- sc$read_h5ad("/data/gpfs/projects/punim0613/JiaDong/PhiSpace/data/Stereo-CITE/CITE-seq/spleen_lymph_111.h5ad")
cMat <- t(ref1$obsm['protein_expression'])
protNames <- ref1$uns$protein_names
protNames <- sub("ADT_", "", protNames)
protNames <- sub("_A.*", "", protNames)
dimnames(cMat) <- list(
  protNames, ref1$obs_names$format()
)
ref_206_ADT <- SingleCellExperiment::SingleCellExperiment(
  list(counts = cMat),
  colData = ref1$obs
)
cMat <- t(ref1$X)
cMat <- Matrix(cMat, sparse = T)
dimnames(cMat) <- list(
  ref1$var_names$format(), ref1$obs_names$format()
)
ref_206_RNA <- SingleCellExperiment::SingleCellExperiment(
  list(counts = cMat),
  colData = ref1$obs
)
unique(ref_206_ADT$cell_types)
toDelete <- c(
  "B-CD8 T cell doublets", "B doublets", "T doublets", "B-macrophage doublets",
  "B-CD4 T cell doublets", "Low quality B cells", "Low quality T cells"
)
idx <- colnames(ref_206_ADT)[!(ref_206_ADT$cell_types %in% toDelete)] 
qsave(
  list(
    ADT = ref_206_ADT[,idx],
    RNA = ref_206_RNA[,idx]
  ),
  paste0(
    dat_dir, "data/Stereo-CITE/CITE-seq/ref_206_list.qs"
  )
)
cMat <- t(ref2$obsm['protein_expression'])
protNames <- ref2$uns$protein_names
protNames <- sub("ADT_", "", protNames)
protNames <- sub("_A.*", "", protNames)
dimnames(cMat) <- list(
  protNames, ref2$obs_names$format()
)
ref_111_ADT <- SingleCellExperiment::SingleCellExperiment(
  list(counts = cMat),
  colData = ref2$obs
)
cMat <- t(ref2$X)
cMat <- Matrix(cMat, sparse = T)
dimnames(cMat) <- list(
  ref2$var_names$format(), ref2$obs_names$format()
)
ref_111_RNA <- SingleCellExperiment::SingleCellExperiment(
  list(counts = cMat),
  colData = ref2$obs
)
unique(ref_111_ADT$cell_types)
toDelete <- c(
  "B-CD8 T cell doublets", "B doublets", "T doublets", "B-macrophage doublets",
  "B-CD4 T cell doublets", "Low quality B cells", "Low quality T cells"
)
idx <- colnames(ref_111_ADT)[!(ref_111_ADT$cell_types %in% toDelete)] 
qsave(
  list(
    ADT = ref_111_ADT[,idx],
    RNA = ref_111_RNA[,idx]
  ),
  paste0(
    dat_dir, "data/Stereo-CITE/CITE-seq/ref_111_list.qs"
  )
)
# Combine RNA parts of two CITE-seq (for RNA only annotation)
ref206path <- paste0(
  dat_dir, "data/Stereo-CITE/CITE-seq/ref_206_list.qs"
)
ref111path <- paste0(
  dat_dir, "data/Stereo-CITE/CITE-seq/ref_111_list.qs"
)
ref206 <- qread(ref206path)$RNA
ref111 <- qread(ref111path)$RNA
ref206 <- logTransf(ref206, targetAssay = "log1p", use_log1p = TRUE)
ref111 <- logTransf(ref111, targetAssay = "log1p", use_log1p = TRUE)
reference <- cbind(ref206, ref111)
reference$dataset <- rep(
  c("data206", "data111"),
  c(ncol(ref206), ncol(ref111))
)
reference$cell_types <- as.character(reference$cell_types)
qsave(
  reference,
  paste0(
    dat_dir, "data/Stereo-CITE/CITE-seq/refComboRNA.qs"
  )
)



### Single neutrophil data from Ng et al. (2024)
rawRefPath <- paste0(
  dat_dir, "data/Neutrophils/GSE244536/scRNA-seq"
)
metaRefPath <- paste0(
  dat_dir, "data/Neutrophils/GSE244536/scRNA-seq/GSE243466_Metadata_df.csv"
)
refPath <- paste0(
  dat_dir, "data/Neutrophils/GSE244536/scRNA-seq/reference.qs"
)
seu <- Read10X(rawRefPath)
metaDat <- read.csv(metaRefPath, row.names = 1)
reference <- SingleCellExperiment::SingleCellExperiment(
  list(counts = seu[,rownames(metaDat)]),
  colData = metaDat
)
qsave(reference, refPath)





























