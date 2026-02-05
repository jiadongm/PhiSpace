# CLAUDE.md

## Project Overview

PhiSpace is an R package for continuous cell state annotation of single-cell and spatial omics data. It uses partial least squares (PLS) regression to phenotype cells on a continuum rather than assigning discrete labels. The R package source lives in `pkg/`.

## Repository Structure

```
PhiSpace/
├── pkg/                    # R package source
│   ├── R/                  # R source files (41 files)
│   ├── man/                # roxygen2-generated documentation
│   ├── vignettes/          # 7 vignettes (Rmd)
│   ├── DESCRIPTION         # Package metadata
│   └── NAMESPACE           # Exports/imports (auto-generated)
├── Test/                   # Manual test scripts and data (git-ignored)
│   ├── data/               # Test datasets (.qs, .rds)
│   └── test_cellTypeThreshold.R
├── docs/                   # pkgdown website
├── figs/                   # README figures
├── sourceCode_PhiSpaceMultiomics/
└── sourceCode_PhiSpaceST/
```

## Key Source Files

- `pkg/R/PhiSpaceR.R` — Main user-facing `PhiSpace()` wrapper (handles single/multiple references, stores normalised scores in `reducedDim`)
- `pkg/R/PhiSpaceR_1ref.R` — Core single-reference implementation `PhiSpaceR_1ref()`. Returns both raw (`PhiSpaceScore`, `YrefHat`) and normalised (`PhiSpaceNorm`, `YrefHatNorm`) scores. Supports `cellTypeThreshold` for filtering rare cell types.
- `pkg/R/tunePhiSpace.R` — Parameter tuning via cross-validation. Also supports `cellTypeThreshold`.
- `pkg/R/codeY.R` — Converts categorical phenotypes to dummy matrices
- `pkg/R/mvr.R` — Internal PLS/PCA regression fitting
- `pkg/R/superPC.R` — Supervised PCA/PLS model builder
- `pkg/R/phenotype.R` — Prediction function
- `pkg/R/normPhiScores.R` — Score normalization (column-wise median centering or row-wise min-max)
- `pkg/R/rankFeatures.R` — Feature importance ranking
- `pkg/R/clusterPhiSpace.R` — K-means clustering on PhiSpace scores
- `pkg/R/spatialSmoother.R` — KNN-based spatial smoothing
- `pkg/R/findNiches.R` — Spatial niche identification
- `pkg/R/vizSpatial.R` — Spatial visualization (`VizSpatial()`)
- `pkg/R/saveCellTypeMaps.R` — Batch spatial heatmap export
- `pkg/R/utils.R` — Shared utilities (color palette, censoring, scaling, etc.)

## Core Function Architecture

```
PhiSpace()                          # User-facing wrapper
 └─ PhiSpaceR_1ref()                # Per-reference core (called in a loop for multi-ref)
     ├─ codeY()                     # Phenotypes → dummy matrix
     ├─ cellTypeThreshold filtering # Optional: remove rare cell types
     ├─ mvr() / SuperPC()           # PLS/PCA model fitting
     ├─ phenotype()                 # Project query onto model
     └─ normPhiScores()             # Normalise scores (returned as PhiSpaceNorm)
```

`PhiSpaceR_1ref()` returns a list with both raw and normalised scores:
- `PhiSpaceScore` / `PhiSpaceNorm` — query scores (matrix or list of matrices)
- `YrefHat` / `YrefHatNorm` — reference predictions

`PhiSpace()` uses the pre-normalised scores from `PhiSpaceR_1ref()` and stores them in `reducedDim(query, "PhiSpace")`.

## Build & Development Commands

All commands should be run from the repo root. The package source is in `pkg/`.

```bash
# Generate/update roxygen2 documentation (man pages + NAMESPACE)
Rscript -e 'devtools::document("pkg")'

# Run R CMD check
Rscript -e 'devtools::check("pkg")'

# Run tests
Rscript -e 'devtools::test("pkg")'

# Install locally
Rscript -e 'devtools::install("pkg")'

# Build pkgdown site
Rscript -e 'pkgdown::build_site("pkg")'

# Run manual integration test (requires Test/data/)
Rscript Test/test_cellTypeThreshold.R
```

## Code Conventions

- **Documentation**: roxygen2 with markdown enabled (`Roxygen: list(markdown = TRUE)`). Man pages in `man/` are auto-generated — edit roxygen comments in `R/*.R` files, then run `devtools::document()`. Use `PhiSpaceR_1ref.R` as a reference for detailed documentation style (structured `@return` with `\describe`, `@details`, `@seealso`, `@references`).
- **Exports**: Declared via `@export` roxygen tag. NAMESPACE is auto-generated.
- **Core data structures**: `SingleCellExperiment` (SCE) and `SpatialExperiment` (SPE) from Bioconductor. Results stored in `reducedDim()` slots and `colData()`.
- **S4 generics**: `cellTypeCorMat` is an S4 generic with methods for matrix, data.frame, and SCE.
- **S3 class**: `PhiSpaceClustering` with print/summary/plot methods (in `clusterPhiSpace.R`).
- **Parameter validation**: Use `stop()` for errors, `warning()` for warnings, `message()` for informational output. Validate new parameters early in the function body (see `cellTypeThreshold` validation pattern).
- **Input flexibility**: Many functions accept both single objects and lists of objects (e.g., `PhiSpace()` accepts single or multiple references/queries). When a return value may be a matrix or list of matrices, handle both cases explicitly.

## Dependencies

- **R >= 4.5.0**
- **Bioconductor**: SummarizedExperiment, SingleCellExperiment, SpatialExperiment, scran, scuttle, S4Vectors, ComplexHeatmap
- **CRAN**: Matrix, ggplot2, ggpubr, dplyr, plyr, magrittr, irlba, FNN, cluster, umap, kerndwd
- **GitHub**: vizOmics (ByronSyun/vizOmics)

## Testing

- **Unit tests**: `testthat` (edition 3). Test files in `pkg/tests/`. Run with `devtools::test("pkg")`.
- **Integration tests**: `Test/test_cellTypeThreshold.R` exercises `PhiSpaceR_1ref()` and `PhiSpace()` with CosMx lung data at multiple `cellTypeThreshold` values, then saves spatial heatmaps via `saveCellTypeMaps()`. Requires `Test/data/ref_list.qs` and `Test/data/Lung5_Rep1.rds`. The `Test/` directory is git-ignored.

## Known Check Notes

- Hidden `.claude` directory in `pkg/` triggers a NOTE (not an error)
- `plot.PhiSpaceClustering` has ggplot2 NSE binding NOTEs (cosmetic, pre-existing)
