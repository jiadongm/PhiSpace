PhiSpaceAssay <- "rank"
phenotypes <- c("Cell Type", "Sample Source")
PhiMethod <- "PLS"
tune_res <- tunePhiSpace(
  reference = reference,
  assayName = PhiSpaceAssay,
  phenotypes = phenotypes,
  regMethod = PhiMethod
)

# In the paper we selected ncomp=30 and nfeat=265

saveRDS(
  tune_res$selectedFeat,
  paste0(dat_dir, "output/CaseDC/ref_dc_test.rds")
)

