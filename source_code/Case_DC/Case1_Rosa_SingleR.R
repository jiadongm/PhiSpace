
### Query
Rosa0 <- readRDS(
  paste0(
    dat_dir,
    "data/Rosa/query_Rosa.rds"
  )
)
set.seed(9048)
idx_list <- 
  lapply(
    1:length(unique(Rosa0$celltype)),
    function(idx){
      lab <- unique(Rosa0$celltype)[idx]
      typeIdx <- which(Rosa0$celltype == lab)
      sample(typeIdx, 73)
    }
  )
idx <- do.call(c, idx_list)
Rosa0 <- Rosa0[,idx]




## Reference
reference <- readRDS(
  paste0(
    dat_dir,
    "data/stemformatics/ref_dc.rds"
  )
)
colnames(reference@colData)
# phenotypes <- c("Cell Type", "Sample Source", "Tissue Type", "Disease State", "Activation Status")
phenotypes <- c("Cell Type", "Sample Source")
YY <- codeY(reference, phenotypes)


## Common genes and rank transform
c(reference, query) %<-% KeepCommonGenes(reference, Rosa0)
reference <- RankTransf(reference, "data")
query <- RankTransf(query, "counts")


## SingleR 
singleRassay <- "rank"
sr_re_l <- vector("list", length(phenotypes))
names(sr_re_l) <- phenotypes
for(idx in 1:length(phenotypes)){
  
  YtrainName <- phenotypes[idx]
  sr_train <- trainSingleR(assay(reference, singleRassay),
                           labels = colData(reference)[, YtrainName])
  sr_re <- classifySingleR(assay(query, singleRassay),
                           sr_train,
                           fine.tune = T)
  sr_ref_re <- classifySingleR(assay(reference, singleRassay),
                           sr_train,
                           fine.tune = T)
  
  sr_re_l[[idx]] <- list(ref = sr_ref_re, query = sr_re)
}
saveRDS(sr_re_l, "output/Case1/Rosa2dc_SR_score.rds")




