## ============================================================================
## Run PhiSpace on all 40 Visium NSCLC samples
## Same workflow as CaseVisium/runPhiSpace.R
## ============================================================================

source("scripts/00_setup/paths_loader.R")

suppressPackageStartupMessages({
  library(PhiSpace)
  library(qs)
  library(dplyr)
})

dat_dir <- paths$phispace_data_root
pdc_dir <- file.path(paths$output_root, "pDC")

PhiAssay <- "log1p"

# =============================================================================
# 1. Load reference (same as CaseVisium)
# =============================================================================

YtrainName <- refLabName <- "ann_finest_level"
refPath <- paste0(dat_dir, "data/LungRef/AzimuthLung2.0_sce_0.1sub.qs")
reference <- qread(refPath)

# Feature selection (same as CaseVisium: 500 features)
impScPath <- paste0(dat_dir, "output/Case3/refImpScores.rds")
impScores <- readRDS(impScPath)
selectedFeat <- selectFeat(impScores, 500)$selectedFeat
cat(sprintf("Selected %d features\n", length(selectedFeat)))

# =============================================================================
# 2. Load query list (all 40 samples from prepareData.R)
# =============================================================================

allPath <- file.path(pdc_dir, "data", "allSamples.qs")
qu_list <- qread(allPath)
cat(sprintf("Loaded %d query samples\n", length(qu_list)))

# QC: remove zero-count features
qu_list <- lapply(qu_list, zeroFeatQC)

# =============================================================================
# 3. Run PhiSpace (same parameters as CaseVisium)
# =============================================================================

cat("Running PhiSpace...\n")
PhiRes <- PhiSpaceR_1ref(
  reference,
  qu_list,
  phenotypes = YtrainName,
  refAssay = PhiAssay,
  queryAssay = PhiAssay,
  selectedFeat = selectedFeat,
  regMethod = "PLS",
  center = TRUE,
  scale = FALSE
)

# =============================================================================
# 4. Save results
# =============================================================================

outPath <- file.path(pdc_dir, "output", "PhiRes_allSamples.qs")
dir.create(file.path(pdc_dir, "output"), showWarnings = FALSE)
qsave(PhiRes, outPath)
cat(sprintf("PhiSpace results saved to %s\n", outPath))

# Quick check
cat("\nSample names in results:\n")
cat(paste(names(PhiRes$PhiSpaceScore), collapse = ", "), "\n")
cat(sprintf("\nCell types: %d\n", ncol(PhiRes$PhiSpaceScore[[1]])))
cat(paste(head(colnames(PhiRes$PhiSpaceScore[[1]]), 10), collapse = ", "), "...\n")
