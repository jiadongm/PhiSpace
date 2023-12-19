### Query
if(F){
  ## Rosa query
  data_dir <- "data/Rosa"
  list.files(data_dir)
  data_mat <- Matrix::readMM("data/Rosa/matrix.mtx")
  barcodes <- read.table("data/Rosa/barcodes.tsv")
  celltypes <- read.table("data/Rosa/celltypes.tsv")
  genes <- read.table("data/Rosa/genes.tsv")
  colnames(data_mat) <- barcodes$V1
  rownames(data_mat) <- genes$V1
  colDat <- data.frame(celltype = celltypes$V1)
  rownames(colDat) <- barcodes$V1
  query <- SingleCellExperiment(
    list(counts = data_mat),
    colData = colDat
  )


  # Add labels without differentiation of donors
  labs <- query$celltype
  table(labs)
  if(T){
    newLabs <- rep("Day3", length(labs))
    idx <- grepl("Day6", labs, fixed = TRUE)
    newLabs[idx] <- "Day6"
    idx <- grepl("Day9_DP", labs, fixed = TRUE)
    newLabs[idx] <- "Day9_DP"
    idx <- grepl("Day9_SP", labs, fixed = TRUE)
    newLabs[idx] <- "Day9_SP"
    idx <- grepl("DC1", labs, fixed = TRUE)
    newLabs[idx] <- "DC1"
    idx <- grepl("DC2", labs, fixed = TRUE)
    newLabs[idx] <- "DC2"
    idx <- grepl("HEF", labs, fixed = TRUE)
    newLabs[idx] <- "HEF"
    idx <- grepl("PDC", labs, fixed = TRUE)
    newLabs[idx] <- "pDC"
  }
  query$mainTypes <- newLabs

  saveRDS(query, "output/Case1/query_Rosa.rds")

  ## subsampling
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
  saveRDS(Rosa0, paste0(dat_dir, "output/Case1/query_Rosa_sub.rds"))
}





## Reference
if(F){
  ## DC atlas
  sample_mat <- read.table("data/dc/samples.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
  data_mat <- read.table("data/dc/expression.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
  geneNames <- read.table("data/dc/genes.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
  # two gene IDs correpond to Mar-02 symbol, hence delete it
  data_mat <- data_mat[geneNames$symbol != c("Mar-02","Mar-01"), ]
  geneNames <- geneNames[geneNames$symbol != c("Mar-02","Mar-01"),]
  rownames(data_mat) <- geneNames$symbol
  geneNames <- data.frame(geneNames, row.names = "symbol")
  ref_dc <- SingleCellExperiment(
    list(data = as.matrix(data_mat)),
    colData = sample_mat,
    rowData = geneNames
  )



  ## Shorter names
  c(ref_dc$`Cell Type`, ref_dc$`Activation Status`, ref_dc$`Sample Source`,
    ref_dc$`Tissue Type`, ref_dc$`Platform Category`, ref_dc$`Disease State`) %>%
    unique()
  lookup <- c(
    "plasmacytoid dendritic cell" = "pDC",
    "monocyte" = "mono",
    "DC precursor" = "DC_prec",
    "dendritic cell" = "DC",
    "in vivo (HuMouse)" = "in_vivo_HuMouse",
    "in vitro derived" = "in_vitro_derived",
    "in vivo" = "in_vivo",
    "in vitro" = "in_vitro",
    "ex vivo" = "ex_vivo",
    "bone marrow" = "BM",
    "synovial fluid" = "synovial",
    "small intestine" = "small_int",
    "Pitt-Hopkins syndrome" = "PittHopkins"
  )
  ref_dc$`Cell Type` <-
    stringr::str_replace_all(
      ref_dc$`Cell Type`,
      lookup
    )
  ref_dc$`Sample Source` <-
    stringr::str_replace_all(
      ref_dc$`Sample Source`,
      lookup
    )
  ref_dc$`Tissue Type` <-
    stringr::str_replace_all(
      ref_dc$`Tissue Type`,
      lookup
    )
  ref_dc$`Disease State` <-
    stringr::str_replace_all(
      ref_dc$`Disease State`,
      lookup
    )
  saveRDS(ref_dc, "output/Case1/ref_dc.rds")
}
