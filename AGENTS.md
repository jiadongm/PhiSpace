# AGENTS.md

## Project Overview

PhiSpace is an R package for continuous cell state annotation of single-cell and spatial omics data. It uses partial least squares (PLS) regression to phenotype cells on a continuum rather than assigning discrete labels. The R package source lives in `pkg/`.

## Repository Structure

```
PhiSpace/
├── pkg/                    # R package source
│   ├── R/                  # R source files
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

## Vignette Overhaul: Active Project

The active maintenance project is to make every vignette reproducible with the
current package and publish it through GitHub Actions. Treat the explicit
analysis currently shown in each vignette as the reference implementation.

### Compatibility First, Wrappers Last

Use two separate passes:

1. **Compatibility pass**: make the existing analysis run with current data,
   dependencies, and PhiSpace APIs. Preserve explicit calculations unless they
   are actually broken.
2. **Wrapper-modernisation pass**: only after all seven vignettes pass, compare
   explicit blocks with current wrappers such as `PhiSpace()`,
   `clusterPhiSpace()`, `findNiches()`, `spatialSmoother()`,
   `saveCellTypeMaps()`, and `scoreCells()`. Establish numerical or structural
   equivalence before replacing a block. Keep explicit code when it teaches the
   scientific method rather than mere boilerplate.

Do not mix compatibility fixes and wrapper refactors in the same commit. The
working explicit vignette is the oracle for evaluating a wrapper.

### Vignette Status and Recommended Order

| Order | Vignette | Status | Main legacy issue |
| --- | --- | --- | --- |
| done | `StereoSeq.Rmd` | Compatibility-complete; automated and live | qs2 migration and current APIs; fresh and cached builds pass |
| done | `Visium.Rmd` | Compatibility-complete; automated and live | qs2 migration; fresh and cached builds pass |
| 1 | `getting_started.Rmd` | Pending | Core tutorial; hard-coded paths, disabled evaluation, qs caches |
| 2 | `PerturbSeq.Rmd` | Pending | Hard-coded paths, disabled evaluation, qs caches |
| 3 | `BridgeAnnotation.Rmd` | Pending | Hard-coded paths, disabled evaluation, multiple RDS results |
| 4 | `CITE-seq.Rmd` | Pending | Hard-coded paths, disabled evaluation, multimodal RDS results |
| 5 | `CosMx.Rmd` | Pending | Largest remainder; multiple references, qs caches, spatial and multi-sample analysis |

The order may change for scientific priority, but `getting_started` should
normally be next because it exercises the public `PhiSpace()` interface and
parameter tuning. Complete and publish one vignette-sized change at a time.

### Standard Procedure for Each Vignette

1. Record its scientific purpose, input objects, helper scripts, cached results,
   expected figures/tables, random seeds, and meaningful output invariants.
2. Inventory every data path and download. Never depend on a developer home
   directory or private Dropbox mount.
3. For legacy `.qs` files, use old `qs` only for one-time local conversion.
   Save with `qs2::qs_save()`, read back with checksum validation, and compare
   the restored R object with the original. Renaming `.qs` to `.qs2` is invalid.
4. Publish data in a stable shared folder, preserve relative paths, and record
   SHA-256 hashes in `scripts/<slug>-inputs.tsv`.
5. Add `PHISPACE_<SLUG>_DATA`, `PHISPACE_<SLUG>_INPUTS`, and, for a folder ZIP,
   `PHISPACE_<SLUG>_ARCHIVE`. Source helpers with `local = TRUE`.
6. Enable evaluation and fail on chunk errors. Use
   `knitr::opts_chunk$set(error = FALSE)` and render with `lazy = FALSE`.
7. Run a **fresh-analysis build** with cached results deliberately absent, so
   computation errors cannot be hidden by old outputs.
8. Run a separate **cached-result build** with all downloadable outputs present.
   Confirm it renders without loading the retired `qs` namespace.
9. Record runtime, peak RAM, package versions, dimensions/names, and meaningful
   numerical differences from established results.
10. Add `scripts/prepare-<slug>.R`, `scripts/build-<slug>.R`, its manifest, and
    the vignette to Actions only after both local builds pass. Resolve the exact
    workflow dependency set locally.
11. Commit only that vignette source, scripts, manifest, workflow update, and
    publishing notes. Never commit `Test/` caches or staged websites.
12. After the user pushes, check both Actions jobs and inspect the live article.
    Local success does not prove runner resources, permissions, or deployment.

### Serialization and Data Rules

- New vignette work uses `qs2`, tested with qs2 0.3.1. Do not add `qs` to
  Actions or new vignette code; it has been removed from active CRAN.
- Legacy manual tests outside this overhaul may still reference `.qs`. Do not
  migrate unrelated test data unless it is explicitly in scope.
- Keep the `PhiSpace/` clone for the package and its repository files only.
  Put local downloads, input caches, converted objects, logs, staged sites,
  helper scripts and other working files outside the clone: under
  `PkgOverhaul/Test/` or `PkgOverhaul/data/<Vignette>/`. Do not add files to
  `PhiSpace/Test/`. The build scripts default to repository-relative `Test/`
  paths because Actions uses them; locally, override them with the
  `PHISPACE_SITE_DIR` and `PHISPACE_<SLUG>_INPUTS` variables.
- Preparation scripts extract only manifest-listed files and verify every hash.
  Upstream changes must fail until reviewed and intentionally accepted.
- Actions caches verified inputs, not computed results. Build scripts create a
  fresh temporary analysis directory so result-loading branches cannot mask
  failures.

### Current Automation

- Workflow: `.github/workflows/vignettes-pages.yaml`, displayed in Actions as
  **Build and publish vignettes**.
- It runs on pushes to `main`, pull requests, and manual dispatch. Pull requests
  build an artifact but do not deploy. Main pushes deploy through GitHub Pages
  with the built-in token; no personal token is stored as a workflow secret.
- It stages committed `docs/`, then rebuilds every vignette already automated.
  This is important: omitting an automated article from a later deployment
  would copy its old committed HTML and revert the live page.
- Each vignette runs in a separate R process. pkgdown is pinned to the CRAN
  archive URL for version 2.2.0; qs2 is resolved from active CRAN.
- Do not edit or commit generated `docs/` pages for this workflow. Actions
  deploys the staged artifact.
- Full local build:

```bash
export PHISPACE_SITE_DIR=../Test/site-vignettes-next
export PHISPACE_STEREOSEQ_INPUTS=../Test/stereoseq-inputs
export PHISPACE_VISIUM_INPUTS=../Test/visium-inputs
Rscript --vanilla scripts/build-stereoseq.R
Rscript --vanilla scripts/build-visium.R
```

See `WEBSITE.md` for variables and deployment details. As the suite grows,
separate fresh validation from publication: validate changed vignettes from
scratch, but eventually publish the full site from verified cached results to
keep runtime and memory bounded.

### Completed Data and Validation

- Stereo-seq archive:
  `https://www.dropbox.com/scl/fo/4w5vweo2ky2vuf7g591fe/AMDVr2OAUL5W1wO7ATzM754?rlkey=aggfds07sjymsd2aomyepbylv&dl=1`
  with hashes in `scripts/stereoseq-inputs.tsv`.
- Visium archive:
  `https://www.dropbox.com/scl/fo/gsaxu5jex7d8ftrwu3pf0/ANN2z_e4abgrON62BD2TEHk?rlkey=fg94wev8y6cn096khmauyaass&dl=1`
  with hashes in `scripts/visium-inputs.tsv`.
- Stereo-seq: all 22 chunks passed fresh and cached with qs2 0.3.1. Fresh
  execution used about 9.3 GiB peak RAM and 2 minutes 40 seconds locally.
- Visium: all 11 chunks passed fresh and cached with qs2 0.3.1. Fresh execution
  used about 14.3 GB peak RAM and 1 minute 28 seconds; the cached build used
  about 7.0 GB and 34 seconds.
- The combined workflow-order build passed locally without loading `qs`.
- Each build script runs `scripts/check-results.R` on the fresh outputs. It
  compares dimensions, names, score summaries, sorted k-means cluster sizes and
  enrichment scores with `scripts/<slug>-results.tsv` (absolute tolerances of
  about 1e-5 relative to each value; counts and names must match exactly). Set
  `PHISPACE_RESULTS_OUT` to write the observed values. Update an expected value
  only after reviewing why it changed.
- Stereo-seq difference from the established Dropbox results: the vignette
  removes the 1% of bins with the lowest total counts before annotation (15,611
  to 15,454 bins), but the established `output/PhiRes.qs2` was computed on all
  15,611 bins. Rerunning on all bins reproduces it (max absolute difference
  3.3e-6). On the shared bins, fresh scores differ only by a per-column
  constant (SD of differences <= 1.4e-6) from query centering, but the
  max-absolute scaling in `normPhiScores()` changes. PhiSpace niche cluster
  sizes therefore differ by up to 17 bins; barcode and gene-expression
  clusterings and enrichment scores match. The expected values follow the
  vignette code. The cached-result build still uses the older Dropbox files.
- Visium fresh results match the established `combo_PhiRes.qs2` within 1.1e-6
  relative difference on all checked metrics.

Relevant commits:

- `a1080bc` — fix the Stereo-seq assay argument.
- `a2d20e6` — add Stereo-seq Actions publishing.
- `d61c249` — migrate Stereo-seq to qs2.
- `797648c` — migrate and automate Visium and combine publishing.

Handoff checkpoint on 2026-10-06: commit `797648c` was pushed; Actions run
`37424544829` completed successfully, and the live Visium page was verified to use
`qs2` and `.qs2` files. Future sessions should still recheck current Actions state.

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

- `plot.PhiSpaceClustering` has ggplot2 NSE binding NOTEs (cosmetic, pre-existing)
