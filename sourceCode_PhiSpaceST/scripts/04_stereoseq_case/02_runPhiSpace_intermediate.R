## Annotate Mouse4 scRNA-seq (published)
source("scripts/00_setup/paths_loader.R")

pathRefSpleen <- paste0(
  dat_dir, 
  "scRNA-seq/refSpleen_processed.qs"
)
pathRefBM <- paste0(
  dat_dir, 
  "scRNA-seq/refBM_processed.qs"
)
pathRefNeutro <- file.path(paths$data_root, "Neutrophils/GSE244536/scRNA-seq/reference.qs")
pathRefCITE <- file.path(paths$data_root, "Stereo-CITE/CITE-seq/refComboRNA.qs")

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







