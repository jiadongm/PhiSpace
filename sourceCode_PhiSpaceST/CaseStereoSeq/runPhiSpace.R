## Annotate Mouse4 scRNA-seq (published)
pathRefSpleen <- paste0(
  dat_dir, 
  "scRNA-seq/refSpleen_processed.qs"
)
pathRefBM <- paste0(
  dat_dir, 
  "scRNA-seq/refBM_processed.qs"
)
pathRefNeutro <- "/data/gpfs/projects/punim0613/JiaDong/PhiSpace/data/Neutrophils/GSE244536/scRNA-seq/reference.qs"
pathRefCITE <- "/data/projects/punim0613/JiaDong/PhiSpace/data/Stereo-CITE/CITE-seq/refComboRNA.qs"

if(F){ # If haven't done this already, change F to T
  refSpleen <- qread(pathRefSpleen)
  refSpleen <- scranTransf(refSpleen)
  refBM <- qread(pathRefBM)
  refNeutro <- qread(pathRefNeutro)
  refNeutro <- zeroFeatQC(refNeutro)
  refNeutro <- scranTransf(refNeutro)
  refCITE <- qread(pathRefCITE)
  # Only use spleen immune cells in refCITE
  refCITE <- refCITE[,refCITE$hash_id == "Spleen"]
  refCITE <- scranTransf(refCITE)
  
  
  ## Check Batch effects in references
  # There is mild batch effects in refCITE; all other references seem ok
  pcaRes <- getPC(t(assay(refSpleen, "logcounts")), ncomp = 20)
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refSpleen$orig.ident)
  pcaRes <- getPC(t(assay(refBM, "rank")), ncomp = 20)
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refBM$cellType)
  pcaRes <- getPC(t(assay(refNeutro, "log1p")), ncomp = 20)
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refNeutro$ManuscriptClusters)
  pcaRes <- getPC(t(assay(refCITE, "log1p")), ncomp = 20)
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refCITE$cell_types) 
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refCITE$batch_indices)
  
  ## Batch correction for refCITE
  set.seed(89230)
  refCITEcorct <- batchelor::fastMNN(
    refCITE, batch = refCITE$batch_indices, BSPARAM=BiocSingular::RandomParam()
  ) 
  pcaRes <- getPC(t(assay(refCITEcorct, "reconstructed")), ncomp = 20)
  matrixPlot(pcaRes$scores, comp_idx = 1:4, colBy = refCITEcorct$batch) 
  cMat <- assay(refCITEcorct, "reconstructed")
  cMat <- Matrix(cMat, nrow = nrow(cMat), ncol = ncol(cMat), dimnames = dimnames(cMat))
  assay(refCITE, "logcounts") <- cMat
  
  qsave(refSpleen, pathRefSpleen)
  qsave(refBM, pathRefBM)
  qsave(refNeutro, pathRefNeutro)
  qsave(refCITE, pathRefCITE)
} else {
  
  refSpleen <- qread(pathRefSpleen)
  refBM <- qread(pathRefBM)
  refNeutro <- qread(pathRefNeutro)
  refNeutro <- zeroFeatQC(refNeutro)
  refNeutro <- scranTransf(refNeutro)
  refCITE <- qread(pathRefCITE)
  # Only use spleen immune cells in refCITE
  refCITE <- refCITE[,refCITE$hash_id == "Spleen"]
  refCITE <- scranTransf(refCITE)
}

if(F){
  query <- qread(
    paste0(dat_dir, "mouse4_scRNAseq_sce.qs")
  )
  query <- scranTransf(query)
  query <- RankTransf(query)
  qsave(query, paste0(dat_dir, "mouse4_scRNAseq_sce.qs"))
} else {
  
  query <- qread(
    paste0(dat_dir, "mouse4_scRNAseq_sce.qs")
  )
}

pathPhiSc <- paste0(
  dat_dir,
  "output/Mouse4_scRNA-seq_PhiSc_list.rds"
)
impScPath <- paste0(
  dat_dir,
  "output/ImpScores_for_Mouse4_scRNA-seq.rds"
)

PhiSc_list <- impScores_list <- vector("list", 4) %>%
  `names<-`(c("Spleen", "BM", "Neutro", "CITE"))

# Spleen ref
PhiSpaceAssay <- "logcounts"
YtrainName <- "celltypes"
tuneRes <- tunePhiSpace(
  reference = refSpleen,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[1]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refSpleen), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refSpleen, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$Spleen <- normPhiScores(PhiRes$PhiSpaceScore)

# BM ref
PhiSpaceAssay <- "rank"
YtrainName <- "cellType"
tuneRes <- tunePhiSpace(
  reference = refBM,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[2]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refBM), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refBM, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$BM <- normPhiScores(PhiRes$PhiSpaceScore)

# Neutro
PhiSpaceAssay <- "logcounts"
YtrainName <- "ManuscriptClusters"
tuneRes <- tunePhiSpace(
  reference = refNeutro,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[3]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 600)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refNeutro), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refNeutro, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$Neutro <- normPhiScores(PhiRes$PhiSpaceScore)

# CITE
PhiSpaceAssay <- "logcounts"
YtrainName <- "cell_types"
tuneRes <- tunePhiSpace(
  reference = refCITE,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[4]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refCITE), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refCITE, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list[[4]] <- normPhiScores(PhiRes$PhiSpaceScore)


saveRDS(PhiSc_list, pathPhiSc)
saveRDS(impScores_list, impScPath)







#### Old
pathRefSpleen <- paste0(
  dat_dir, 
  "scRNA-seq/refSpleen_processed.qs"
)
pathRefBM <- paste0(
  dat_dir, 
  "scRNA-seq/refBM_processed.qs"
)
pathRefNeutro <- "/data/gpfs/projects/punim0613/JiaDong/PhiSpace/data/Neutrophils/GSE244536/scRNA-seq/reference.qs"
pathRefCITE <- "/data/projects/punim0613/JiaDong/PhiSpace/data/Stereo-CITE/CITE-seq/refComboRNA.qs"
refSpleen <- qread(pathRefSpleen)
refSpleen <- logTransf(refSpleen, targetAssay = "log1p", use_log1p = TRUE)
refBM <- qread(pathRefBM)
refNeutro <- qread(pathRefNeutro)
refNeutro <- zeroFeatQC(refNeutro)
refNeutro <- logTransf(refNeutro, targetAssay = "log1p", use_log1p = TRUE)
refCITE <- qread(pathRefCITE)
# Only use spleen immune cells in refCITE
refCITE <- refCITE[,refCITE$hash_id == "Spleen"]
refFibro <- readRDS(paste0(
  "/data/projects/punim0613/JiaDong/Dawson/MouseSpleen/",
  "scRNA-seq/Alexandre/reference.rds"
))



## Investigate neutrophils in CITE and Spleen datasets, are they very different?
if(F){
  # Calculate maturation scores
  maturationModule <- readxl::read_excel("PhiSpace/CaseStereoSeq/science.adf6493_table_s1.xlsx")$GeneID # See Ng et al. (2024). Nature.
  library(Seurat)
  Seu_CITE <- CreateAssayObject(
    counts = assay(refCITE, "counts")
  ) %>%
    CreateSeuratObject(
      meta.data = colData(refCITE) %>% as.data.frame()
    ) %>% NormalizeData()
  Seu_Spleen <- CreateSeuratObject(
    counts = assay(refSpleen, "counts"), meta.data = colData(refSpleen) %>% as.data.frame()
  ) %>% NormalizeData()
  Seu_Neutro <- CreateSeuratObject(
    counts = assay(refNeutro, "counts"), meta.data = colData(refNeutro) %>% as.data.frame()
  ) %>% NormalizeData()
  Seu_CITE <- Seurat::AddModuleScore(
    Seu_CITE, features = list(maturationModule), 
    name = "MaturationScore", nbin = 24
  ) 
  Seu_Spleen <- Seurat::AddModuleScore(
    Seu_Spleen, features = list(maturationModule), 
    name = "MaturationScore", nbin = 24
  ) 
  Seu_Neutro <- Seurat::AddModuleScore(
    Seu_Neutro, features = list(maturationModule), 
    name = "MaturationScore", nbin = 24
  ) 
  Seu_CITE_neu <- Seu_CITE[,Seu_CITE$cell_types == "Neutrophils"]
  Seu_Spleen_neu <- Seu_Spleen[,Seu_Spleen$celltypes == "Neutrophil"]
  
  plot_dat <- data.frame(
    MaturationScore = c(
      Seu_CITE_neu$MaturationScore1, Seu_Spleen_neu$MaturationScore1,
      Seu_Neutro$MaturationScore1[Seu_Neutro$ManuscriptClusters == "preNeu"],
      Seu_Neutro$MaturationScore1[Seu_Neutro$ManuscriptClusters == "IMM 1"],
      Seu_Neutro$MaturationScore1[Seu_Neutro$ManuscriptClusters == "MAT 1"],
      Seu_Neutro$MaturationScore1[Seu_Neutro$ManuscriptClusters == "MAT 2"],
      Seu_Neutro$MaturationScore1[Seu_Neutro$ManuscriptClusters == "T3"]
    ),
    subType = rep(
      c("CITE", "Spleen", "preNeu", "IMM 1", "MAT 1", "MAT 2", "T3"),
      c(
        ncol(Seu_CITE_neu), ncol(Seu_Spleen_neu),
        sum(Seu_Neutro$ManuscriptClusters == "preNeu"),
        sum(Seu_Neutro$ManuscriptClusters == "IMM 1"),
        sum(Seu_Neutro$ManuscriptClusters == "MAT 1"),
        sum(Seu_Neutro$ManuscriptClusters == "MAT 2"),
        sum(Seu_Neutro$ManuscriptClusters == "T3")
      )
    )
  )
  
  plot_dat %>%
    ggplot(
      aes(subType, MaturationScore)
    ) +
    geom_violin(
      aes(
        fill = subType
      )
    ) 
  
  summary(Seu_CITE_neu$MaturationScore1)
  summary(Seu_Spleen_neu$MaturationScore1)
}






pathPhiSc <- paste0(
  dat_dir,
  "output/PhiSc_list.rds"
)
impScPath <- paste0(
  dat_dir,
  "output/DawsonImpScores.rds"
)

PhiSc_list <- impScores_list <- vector("list", 5) %>%
  `names<-`(c("Spleen", "BM", "Neutro", "CITE", "Fibro"))

# Spleen ref
PhiSpaceAssay <- "log1p"
YtrainName <- "celltypes"
tuneRes <- tunePhiSpace(
  reference = refSpleen,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[1]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refSpleen), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refSpleen, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$Spleen <- normPhiScores(PhiRes$PhiSpaceScore)

# BM ref
PhiSpaceAssay <- "rank"
YtrainName <- "cellType"
tuneRes <- tunePhiSpace(
  reference = refBM,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[2]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refBM), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refBM, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$BM <- normPhiScores(PhiRes$PhiSpaceScore)

# Neutro
PhiSpaceAssay <- "log1p"
YtrainName <- "ManuscriptClusters"
tuneRes <- tunePhiSpace(
  reference = refNeutro,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[3]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 600)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refNeutro), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refNeutro, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list$Neutro <- normPhiScores(PhiRes$PhiSpaceScore)

# CITE
PhiSpaceAssay <- "log1p"
YtrainName <- "cell_types"
tuneRes <- tunePhiSpace(
  reference = refCITE,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[4]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refCITE), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refCITE, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list[[4]] <- normPhiScores(PhiRes$PhiSpaceScore)

# Fibro 
PhiSpaceAssay <- "rank"
refFibro <- RankTransf(refFibro, "data")
YtrainName <- "seurat_clusters"
tuneRes <- tunePhiSpace(
  reference = refFibro,
  assayName = PhiSpaceAssay,
  phenotypes = YtrainName,
  tune_ncomp = F,
  tune_nfeat = F
)
impScores <- tuneRes$impScores
impScores_list[[5]] <- impScores
selectedFeat <- selectFeat(impScores, nfeat = 500)$selectedFeat
length(selectedFeat)
length(
  intersect(
    intersect(
      rownames(refFibro), rownames(query)
    ),
    selectedFeat
  )
)
PhiRes <- PhiSpaceR_1ref(
  refFibro, 
  query, 
  phenotypes = YtrainName, 
  PhiSpaceAssay = PhiSpaceAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  scale = FALSE
)
PhiSc_list[[5]] <- normPhiScores(PhiRes$PhiSpaceScore)

saveRDS(PhiSc_list, pathPhiSc)
saveRDS(impScores_list, impScPath)




























