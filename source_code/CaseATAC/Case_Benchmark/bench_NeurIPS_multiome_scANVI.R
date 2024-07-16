## Run on HPC

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
dat_dir <- "~/PhiSpace/"

bridgeRNA_list <- readRDS(paste0(dat_dir, "data/NeurIPS2021/bridgeRNA_list_Seu_cor.rds"))

## Python packages
library(reticulate)
use_condaenv("scarches")
sc <- import("scanpy")
tourch <- import("torch")
sca <- import("scarches")
np <- import("numpy")

################ scANVI ##################
cell_type_key = 'cellType'

batchNames <- names(bridgeRNA_list)
for(refIdx in 1:length(batchNames)){
  
  refBatchName <- batchNames[refIdx]
  reference <- bridgeRNA_list[[refBatchName]]
  queryBatchNames <- setdiff(
    batchNames, refBatchName
  )
  
  # Subsettting HVG
  reference <- NormalizeData(reference)
  reference <- FindVariableFeatures(reference)
  HVGs <- VariableFeatures(reference)
  
  countMat <- t(
    as(
      reference[["RNA"]]$counts,
      "RsparseMatrix"
    )[HVGs, ]
  )
  reference <- sc$AnnData(
    X = countMat,
    obs = reference@meta.data
  )
  reference$var_names <- HVGs
  # set up data
  sca$models$SCVI$setup_anndata(
    reference,
    labels_key=cell_type_key
  )
  
  # Set up model
  vae <- sca$models$SCVI(
    reference,
    n_layers=2L,
    encode_covariates=T,
    deeply_inject_covariates=F,
    use_layer_norm="both",
    use_batch_norm="none"
  )
  # Train
  vae$train()
  
  # Train a scANVI model
  scanvae <- sca$models$SCANVI$from_scvi_model(vae, unlabeled_category = "Unknown")
  scanvae$train(max_epochs=20L)
  
  # Save trained model and HVGs
  ref_path = paste0(
    'PhiSpace/Case_Benchmark/ref_model/',
    refBatchName,
    '/'
  )
  scanvae$save(ref_path, overwrite=T)
  saveRDS(
    HVGs,
    paste0(
      ref_path,
      "HVGs.rds"
    )
  )
  
}


# Map query
scANVI_ann_res <- vector(
  "list", length(batchNames)
) %>%
  `names<-`(batchNames)
for(refIdx in 1:length(batchNames)){
  
  refBatchName <- batchNames[refIdx]
  queryBatchNames <- setdiff(
    batchNames, refBatchName
  )
  ref_path = paste0(
    'PhiSpace/Case_Benchmark/ref_model/',
    refBatchName,
    '/'
  )
  HVGs <- readRDS(
    paste0(
      ref_path, "HVGs.rds"
    )
  )
  
  temp_res <- vector(
    "list", length(queryBatchNames)
  ) %>%
    `names<-`(queryBatchNames)
  for(queryIdx in 1:length(queryBatchNames)){
    
    queryBatchName <- queryBatchNames[queryIdx]
    query <- bridgeRNA_list[[queryBatchName]]
    
    
    countMat <- t(
      as(
        query[["RNA"]]$counts,
        "RsparseMatrix"
      )[HVGs,]
    )
    query <- sc$AnnData(
      X = countMat,
      obs = query@meta.data
    )
    
    # Mapping query
    model = sca$models$SCANVI$load_query_data(
      query,
      ref_path,
      freeze_dropout = T,
    )
    model[["_labeled_indices"]] <- model[["_unlabeled_indices"]]
    model[["_unlabeled_indices"]] <- np$arange(query$n_obs)
    
    model$train(
      max_epochs=100L,
      plan_kwargs=dict(weight_decay=0.0),
      check_val_every_n_epoch=10L,
    )
    
    query_path = paste0(
      'PhiSpace/Case_Benchmark/query_model/',
      refBatchName,
      '/',
      queryBatchName,
      '/'
    )
    model$save(query_path, overwrite = T)
    
    temp_res <- list(
      errSR = PhiSpace:::classErr(
        model$predict(),
        query$obs$cellType %>% as.character()
      )$err
    )
    
  }
  
}


# Calculate errors
scANVI_ann_res <- vector(
  "list", length(batchNames)
) %>%
  `names<-`(batchNames)
for(refIdx in 1:length(batchNames)){
  
  refBatchName <- batchNames[refIdx]
  queryBatchNames <- setdiff(
    batchNames, refBatchName
  )
  ref_path = paste0(
    'PhiSpace/Case_Benchmark/ref_model/',
    refBatchName,
    '/'
  )
  HVGs <- readRDS(
    paste0(
      ref_path, "HVGs.rds"
    )
  )
  
  temp_res <- vector(
    "list", length(queryBatchNames)
  ) %>%
    `names<-`(queryBatchNames)
  for(queryIdx in 1:length(queryBatchNames)){
    
    queryBatchName <- queryBatchNames[queryIdx]
    query <- bridgeRNA_list[[queryBatchName]]
    
    
    countMat <- t(
      as(
        query[["RNA"]]$counts,
        "RsparseMatrix"
      )[HVGs,]
    )
    query <- sc$AnnData(
      X = countMat,
      obs = query@meta.data
    )
    
    query_path = paste0(
      'PhiSpace/Case_Benchmark/query_model/',
      refBatchName,
      '/',
      queryBatchName,
      '/'
    )
    model <- sca$models$SCANVI$load(
      query_path, adata = query
    )
    
    temp_res <- list(
      errSR = PhiSpace:::classErr(
        model$predict(),
        query$obs$cellType %>% as.character()
      )$err
    )
    
  }
  
}


if(F){
  reference_latent = sc$AnnData(scanvae$get_latent_representation())
  reference_latent$obs[["cellType"]] = reference$obs[[cell_type_key]]
  sc$pp$neighbors(reference_latent, n_neighbors=8L)
  sc$tl$leiden(reference_latent)
  sc$tl$umap(reference_latent)
  sc$pl$umap(reference_latent,
             color=c('cellType'),
             frameon=F,
             wspace=0.6,
  )
  PhiSpace:::classErr(
    scanvae$predict(),
    reference$obs$cellType %>% as.character()
  )$err
  
  # viz query
  query_latent = sc$AnnData(model$get_latent_representation())
  query_latent$obs[['cellType']] = query$obs[[cell_type_key]]
  sc$pp$neighbors(query_latent)
  sc$tl$leiden(query_latent)
  sc$tl$umap(query_latent)
  sc$pl$umap(
    query_latent,
    color='cellType',
    frameon=F,
    wspace=0.6,
  )
}

