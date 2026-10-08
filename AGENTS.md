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
2. **Wrapper-modernisation pass** (plan below, under "Wrapper-Modernisation Plan"): only after all seven vignettes pass, compare
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
| done | `getting_started.Rmd` | Compatibility-complete; automated and live | qs2 migration; fresh and cached builds pass |
| done | `BridgeAnnotation.Rmd` | Compatibility-complete; automated from cached results; live | Private paths, file-name case; fresh build needs 17-24 GB RAM |
| done | `PerturbSeq.Rmd` | Compatibility-complete; automated and live | qs2 migration and new public folder; private `utils.R`; celldex reference shipped as an input |
| done | `CITE-seq.Rmd` | Compatibility-complete; automated from cached results; live | Private paths; fresh RNA branch saved the wrong object name; fresh build needs about 65 GB RAM |
| done | `CosMx.Rmd` | Compatibility-complete; automated and live | qs2 migration in the existing Dropbox folder; five download links replaced by one folder link |

The order may change for scientific priority, but `getting_started` should
normally be next because it exercises the public `PhiSpace()` interface and
parameter tuning. Complete and publish one vignette-sized change at a time.

### Pending User Actions (TODO)

None. The user pushes each commit; after a push, check both Actions jobs and
the live pages.

Local clean-up was done on 2026-10-07: the user deleted all of
`PkgOverhaul/Test/` (input caches, the qs2 library, logs and scripts) and the
data files under `PkgOverhaul/data/`. Only READMEs, conversion manifests and
sessionInfo files remain in `data/`. The `.qs` to `.qs2` conversion scripts and
`make_stereoseq_replacements.R` were lost with `Test/`. On 2026-10-08 the
READMEs in `data/` were rewritten as records that describe the lost scripts'
method; the manifests still record every hash. For a
local build, download inputs again with `scripts/prepare-<slug>.R` and install
qs2 from CRAN.

### Wrapper-Modernisation Plan (agreed 2026-10-07)

This is pass 2 of "Compatibility First, Wrappers Last". Scope, by expected
value (checked by reading the vignettes and wrapper sources). Status is in the
Decision column.

| Order | Vignette | Explicit block | Wrapper | Decision |
| --- | --- | --- | --- | --- |
| 1 | PerturbSeq | `PhiSpaceR_1ref()` + `normPhiScores()`, one reference | `PhiSpace()` | Done 2026-10-07 (see the PerturbSeq wrapper record below). |
| 2 | CosMx | Loop over 4 lineage references with `PhiSpaceR_1ref()` + `normPhiScores()` and a `(lineage)` name suffix | `PhiSpace(reference = ref_list)` | Done 2026-10-07 (see the CosMx wrapper record below). |
| 2 | CosMx | PCA of scores, `set.seed(94863)`, `kmeans(centers = 9, iter.max = 200, nstart = 50)` | `findNiches()` | Done 2026-10-07, with the vignette's k-means arguments. |
| 3 | StereoSeq | Same PCA + k-means pattern (`centers = 8`) | `findNiches()` | Done 2026-10-07 (see the StereoSeq wrapper record below). |
| 3 | StereoSeq | Single-cell scores as `response` | `PhiSpace(response = ...)` | Done 2026-10-07. |
| 4 | Visium | One reference, list of tissues as query | `PhiSpace(query = query_list)` | Done 2026-10-07 (see the Visium wrapper record below). |
| - | BridgeAnnotation, CITE-seq | Bridge transfer; CITE-seq uses `YrefHat` | - | Keep explicit. The code teaches the method, and new cached results need 17-24 GB and about 65 GB RAM. |
| - | getting_started | Already uses `PhiSpace()` and `tunePhiSpace()` | - | No change. |

`findNiches()` defaults differ from the vignettes: Lloyd algorithm, 20 starts and
500 iterations, against the R default Hartigan-Wong algorithm, 50 starts and
200 iterations. Its default seed (94863) matches. Its PCA is the same `getPC()`
call. With the vignette's arguments it reproduces the CosMx labels exactly; with
its defaults the adjusted Rand index is 0.89. Since 2026-10-08 it reports
progress and niche sizes with `message()` (silence with `verbose = FALSE` or
`suppressMessages()`); the CosMx and Stereo-seq niche chunks set
`message = FALSE`, so fresh and cached pages show the same text. The same
change was made in `zeroFeatQC()`, `spatialSampler()` and `spatialSmoother()`
(whose progress bar still writes to standard output). `clusterPhiSpace.R`
keeps `cat()` because it uses it only in `print()` and `summary()` methods.
Fresh local builds after the change: CosMx 35 checks and 14 of 14 figures
identical, Stereo-seq 48 checks and 23 of 23 figures identical. The kmeans
and irlba warnings stay visible, as in the old explicit vignettes.

`pls.fit()` calls `irlba::partial_eigen()`, which draws a random start vector,
so PhiSpace scores vary between runs by about 1e-8 unless the seed is the same.
Compare explicit and wrapper code from the same seed.

Procedure, one vignette per commit:

1. Run the explicit block and the wrapper on the same inputs in a script
   outside the clone. Compare scores with `all.equal()`, and clusters with an
   adjusted Rand index of 1 (cluster numbers may be permuted).
2. If a wrapper differs, fix the wrapper first, in a separate package commit.
   Testing the wrappers against working code is the main gain of this pass.
3. Replace the block, rebuild fresh, and require the existing
   `scripts/<slug>-results.tsv` checks to pass unchanged.
4. `PhiSpace()` returns the query object, not the `PhiSpaceR_1ref()` list, so
   saved result files change. Regenerate them, update the manifests, and
   confirm each Dropbox upload with the user. Save only what the vignette uses
   (the user chose this for PerturbSeq): usually `reducedDim(query, "PhiSpace")`.
   Result checks that need removed components are dropped; retained checks
   must keep their values.
5. Upload the new result file before the user pushes the commit: the preparation
   script requires every manifest file in the Dropbox folder ZIP. Delete the
   superseded file from Dropbox only after the new page is live.

PerturbSeq wrapper record (2026-10-07):

- `PhiSpace(reference, query, phenotypes = "label.fine", refAssay = "logcounts",
  regMethod = "PLS", nfeat = 500)` matches `PhiSpaceR_1ref()` +
  `normPhiScores()`: normalised query scores within 2.1e-15, raw query scores
  within 3.3e-16, reference predictions within 3.3e-15, identical dimnames.
  No wrapper change was needed.
- The vignette now saves `output/PhiSpaceScores.qs2` (the 28,453 x 15
  normalised score matrix, 3.2 MB) instead of `output/PhiRes.qs2` (9.7 MB). It no
  longer prints the number of selected features (2,821).
- `scripts/perturbseq-results.tsv`: dropped the raw-score, reference, `ncomp`
  and feature-count checks (10 metrics); kept 12 with unchanged expected
  values (observed within 1.5e-15); added `norm_score_mean` and
  `norm_score_sd`. 14 metrics pass in fresh (61 s, 13.3 GB) and cached
  (43 s, 6.5 GB) builds, with identical values.
- All 11 figures are byte-identical to a local render of the previous
  explicit vignette. Against the live page, 10 of 11 are identical; the
  heatmap differs in fine detail, as it does for the explicit vignette built
  locally, so the difference comes from the build machine.
- Uploaded `output/PhiSpaceScores.qs2` to `VignetteData/PerturbSeq/` on
  2026-10-07 (Dropbox hashes match); the public folder download passes the new
  manifest. The new page went live in run `37613768456`; `output/PhiRes.qs2`
  was then already absent from Dropbox when the approved deletion ran.

CosMx wrapper record (2026-10-07):

- `PhiSpace(ref_list, query, ...)`, after `logTransf()` of each reference,
  gives bit-identical scores and names to the explicit four-reference loop
  when both start from the same seed. No wrapper change was needed.
- `findNiches(query, n_niches = 9, ncomp = 25, kmeans_algorithm =
  "Hartigan-Wong", kmeans_nstart = 50L, kmeans_iter = 200L, seed = 94863)`
  gives labels identical to `getPC()` + `kmeans()` and to the established
  cached clustering. The irlba start of the PCA does not change the labels.
- The vignette saves `output/Lung5_Rep1_PhiSpaceScores4Refs.qs2` (98,002 x 43
  normalised scores, 29.7 MB) and `output/Lung5_Rep1_PhiNiches4Refs.qs2`
  (niche labels, 0.15 MB), replacing `CosMxLung5Rep1PhiRes4Refs.qs2` and
  `Lung5_Rep1_PhiClusts4Refs.qs2`. It no longer draws the PCA scree plot
  (only fresh builds drew it) or runs the no-op lineage renaming.
- All 35 checks in `scripts/cosmx-results.tsv` pass with unchanged expected
  values: the check splits the score matrix by its `(lineage)` suffix and
  recomputes the within-cluster sum of squares on 25 principal components.
  Fresh: 111 s, 10.3 GB; cached: 69 s, 9.5 GB. All 14 figures are
  byte-identical between fresh and cached builds and to a local render of
  the previous explicit vignette.
- The multi-sample section stays explicit (saved all-lung results and a
  PLS-DA with `mvr()`).
- Uploaded both new files to `VignetteData/CosMx/output/` on 2026-10-07
  (Dropbox hashes match); the public folder download passes the new manifest
  (16 files). The new page went live in run `37616620258`; the superseded
  `CosMxLung5Rep1PhiRes4Refs.qs2` and `Lung5_Rep1_PhiClusts4Refs.qs2` were
  then deleted from Dropbox (user-approved), and the public download still
  passes the manifest.

StereoSeq wrapper record (2026-10-07):

- The bridge annotation uses `PhiSpace(querySC, query, response = ...,
  refAssay = "log1p", nfeat = 500, regMethod = "PLS", scale = FALSE)`, and the
  PhiSpace niches use `findNiches(query, n_niches = 8, ncomp = 30,
  kmeans_algorithm = "Hartigan-Wong", kmeans_nstart = 50L, kmeans_iter = 200L,
  seed = 94863)`. Rendered fresh from the same inputs, the old explicit
  vignette and the new one give bit-identical normalised scores, identical
  niche labels, and identical barcode clusters, gene-expression clusters and
  niche enrichment scores. No wrapper change was needed.
- The barcode k-means (k = 2 to 10) and the gene-expression k-means stay
  explicit: they do not cluster PhiSpace scores, and the barcode chunk
  explains the k-means settings.
- The vignette saves `output/PhiSpaceScores.qs2` (15,454 x 70 normalised
  scores, 7.7 MB) and `output/PhiNiches.qs2` (niche labels, 65 KB), replacing
  `PhiRes.qs2` (151 MB) and `PhiClustRes.qs2`. The Stereo-seq manifest lists
  inputs only, so Actions does not need these files; they serve readers.
  The PCA scree plot and the printed explained-variance vector (fresh builds
  only) are gone.
- `scripts/stereoseq-results.tsv`: the raw-score checks (mean, SD, six
  per-cell-type SDs) and the selected-feature count (3,517) were replaced by
  the same summaries of the normalised scores (8 metrics). The other 40 checks
  keep their expected values. 48 checks pass in a fresh build (167 s,
  9.8 GB) and in a render from all saved results (68 s, 3.5 GB) that
  rewrote no saved file.
- All 23 figures are byte-identical to the old explicit render; its 24th
  figure was the scree plot. The saved-results render gives the same 23.
- Uploaded both new files to `VignetteData/StereoSeq/output/` on 2026-10-07
  (Dropbox hashes match; the public folder ZIP serves the same SHA-256).
  The new page went live in run `37618603280`; `PhiRes.qs2` and
  `PhiClustRes.qs2` were then deleted from Dropbox (user-approved,
  2026-10-07).

Visium wrapper record (2026-10-07):

- `PhiSpace(reference, query_list, phenotypes = "ann_finest_level",
  refAssay = "log1p", selectedFeat = selectedFeat, regMethod = "PLS")`
  normalises each sample separately, as the explicit code did. Rendered fresh
  from the same inputs, the old and new vignettes give bit-identical
  normalised scores and correlation matrices. No wrapper change was needed.
- The vignette saves `output/PhiSpaceScores.qs2` (list of 18 normalised score
  matrices, 18.9 MB), replacing `combo_PhiRes.qs2` (raw scores, 17.2 MB). The
  Visium manifest lists inputs only, so Actions does not need the file.
- The P11_T3 B-cell map now shows normalised scores (legend 0 to 1 instead of
  raw -1 to 0; same spatial pattern). The other 7 figures are byte-identical
  to the old render.
- `scripts/visium-results.tsv`: the 36 per-sample raw-score checks and the
  raw B-cell SD were replaced by the same summaries of the normalised scores
  (37 metrics); the other 21 keep their expected values. 58 checks pass in a
  fresh build (94 s, 14.4 GB) and in a render from the saved scores (36 s,
  7.0 GB) that rewrote no saved file and gave the same 8 figures.
- Uploaded `PhiSpaceScores.qs2` to `VignetteData/Visium/output/` on
  2026-10-07 (Dropbox hash matches). The new page went live in run
  `37685636487` (commit `4026b0e`); `combo_PhiRes.qs2` was then deleted from
  Dropbox (user-approved, 2026-10-08). `output/` now holds
  `PhiSpaceScores.qs2` and `refImpScores.qs2`.
- The first Actions run (`37620847488`, commit `9b17582`) failed in the
  Visium build: 4 of 18 per-sample normalised means differed from the local
  values by 1.0e-7 to 1.7e-7, just over their 1e-7 tolerance floor. A
  normalised mean (about -0.001 to 0.01) is small relative to the spread of
  the scores, so a tolerance relative to the mean is too tight. Every
  normalised-mean tolerance (Visium, Stereo-seq, PerturbSeq) now equals the
  tolerance of the matching SD (1e-5 x SD, 1.6e-6 to 2.5e-6). The local runs
  still pass 58, 48 and 14 checks.

### Dropbox Access

The user configured an rclone remote named `Dropbox` (capital D) on Spartan on
2026-10-07. Vignette data are under
`Dropbox:Research_projects/PhiSpace/VignetteData/<Folder>/`. Writing to Dropbox
publishes data: confirm each upload with the user, back up replaced files to
`Dropbox:Research_projects/PhiSpace/VignetteData_backups/` (outside the shared
folders), compare `rclone hashsum dropbox` of remote and local files, and then
download through the public link and check the SHA-256 manifest.

Done on 2026-10-07 with this procedure:

- Clean-up after all seven vignettes were live (user-approved): deleted the 15
  old `.qs` files in `CosMx/`, the superseded `DC/` and `Perturb-seq/` folders,
  `VignetteData_backups/`, `Visium/data/LungRef/AzimuthLungMarkers.qs` and
  `StereoSeq/.Rhistory` (about 3.5 GB). Old individual file links no longer
  work. `VignetteData/` now has one folder per vignette (`ATAC`, `CITE`,
  `CosMx`, `PerturbSeq`, `StereoSeq`, `Visium`, `getting_started`) and no `.qs`
  files. The local copies were deleted in the local clean-up; the old `.qs`
  originals are no longer kept anywhere.

- Stereo-seq: `output/PhiRes.qs2`, `PhiClustRes.qs2` and `cloneKDEres.qs2`
  replaced by the files in `PkgOverhaul/data/StereoSeq/replacement-2026-10-06/`;
  the old files are in `VignetteData_backups/StereoSeq-2026-10-06/output/`.
  The existing shared link still serves the folder, the replacements match
  `replacement-manifest.csv`, and the five CI inputs still match
  `scripts/stereoseq-inputs.tsv`.
- Getting Started: new folder `VignetteData/getting_started/` with the three
  qs2 inputs and a new public link; the public download matches
  `scripts/getting-started-inputs.tsv`. The old `VignetteData/DC/` folder and
  its individual file links are unchanged.
- PerturbSeq: new folder `VignetteData/PerturbSeq/` with `utils.R`,
  `data/sceStim.qs2`, `data/ref_sce.qs2` and `output/PhiRes.qs2` (1.16 GB) and
  a new public link; the public download matches
  `scripts/perturbseq-inputs.tsv`. The old `VignetteData/Perturb-seq/` folder
  and its public `sceStim.qs` file link are unchanged.
- CosMx: the 15 qs2 files were added next to the old `.qs` files in the
  existing `VignetteData/CosMx/` folder (same layout; `CosMx_utils.R` was
  already identical), and a public folder link was created; the public download
  matches `scripts/cosmx-inputs.tsv`. The user chose one folder per vignette:
  remove the old `.qs` files once the new page is live.
- 2026-10-08 (user-approved): replaced the public `Visium/README.md` and
  `StereoSeq/README.md` with the rewritten records from `PkgOverhaul/data/`
  (Dropbox hashes match the local files). The old versions, which still gave
  upload instructions and listed `combo_PhiRes.qs2`, are in
  `VignetteData_backups/<Folder>-2026-10-08/README.md`.

### Future Work

- **Sparse-aware `mvr()`.** `mvr()` (and `phenotype()` through `scale()`)
  converts sparse predictor matrices to dense. For BridgeAnnotation peaks this
  allocates 8.6 GiB and fresh builds peaked at 17.3 and 24.4 GB RAM. A sparse-aware fit
  (for example, uncentred cross-products on `dgCMatrix`) could let such
  vignettes build fresh in Actions. This is a package change: validate it
  against the BridgeAnnotation expected results before switching that vignette
  back to fresh builds.

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
  `PhiSpace/Test/`; its earlier contents were moved to `PkgOverhaul/Test/` on
  2026-10-06, and `PkgOverhaul/Test/` was deleted on 2026-10-07 (see Pending
  User Actions). The build scripts default to repository-relative `Test/`
  paths because Actions uses them; locally, override them with the
  `PHISPACE_SITE_DIR` and `PHISPACE_<SLUG>_INPUTS` variables.
- Preparation scripts extract only manifest-listed files and verify every hash.
  Upstream changes must fail until reviewed and intentionally accepted.
- By default, Actions caches verified inputs, not computed results, and build
  scripts create a fresh temporary analysis directory so result-loading
  branches cannot mask failures.
- Exception: a vignette whose fresh build exceeds the GitHub runner (16 GB RAM
  as of 2026-10; recheck) is published from verified cached results. Its
  manifest gives each file a `role` (`input` or `cached`), and its build script
  takes `PHISPACE_<SLUG>_MODE` (`cached` by default, `fresh` locally). A cached
  build stages the hashed results, fails if the render rewrites any of them,
  and still runs the result checks. Run a local fresh build and compare it with
  the cached results before changing a cached file or an expected value.
  Currently this applies to BridgeAnnotation and CITE-seq.

### Current Automation

- Package check: `.github/workflows/R-CMD-check.yaml`, displayed as
  **R CMD check** (added 2026-10-08; see Known Check Notes).
- Workflow: `.github/workflows/vignettes-pages.yaml`, displayed in Actions as
  **Build and publish vignettes**.
- It runs on pushes to `main`, pull requests, and manual dispatch. Pull requests
  build an artifact but do not deploy. Main pushes deploy through GitHub Pages
  with the built-in token; no personal token is stored as a workflow secret.
- It stages committed `docs/`, rebuilds every vignette already automated, and
  then rebuilds the function reference from `pkg/man`
  (`scripts/build-reference.R`, added 2026-10-08). The reference step must run
  after Stereo-seq: `build-stereoseq.R` stages the site itself and refuses an
  existing site folder (run `37714789741` failed when the reference step ran
  first). The reference step deletes the staged `reference/` first,
  so topics removed from `pkg/man` disappear, and rewrites the reference entries
  of the staged `sitemap.xml`. The first local run gave 57 pages: it dropped 13
  pages of functions that moved to vizOmics or were removed (for example
  `align_clusters`, `matrixPlot`, `plotSankey`) and added `piScore` and
  `scoreCells`. Examples are not run; every example is in `\dontrun{}`.
  This is important: omitting an automated article from a later deployment
  would copy its old committed HTML and revert the live page.
- The runner is pinned to `ubuntu-24.04`, because `ubuntu-latest` moves to
  Ubuntu 26 from 2026-10-19 and could change binary packages and numerical
  results. Actions run on Node.js 24: `checkout@v5`, `cache@v5`,
  `upload-artifact@v6`, `upload-pages-artifact@v5` (excludes hidden files such
  as `.nojekyll`, which Actions-deployed Pages do not need) and
  `deploy-pages@v5`.
- Each vignette runs in a separate R process. pkgdown is pinned to the CRAN
  archive URL for version 2.2.0; qs2 is resolved from active CRAN.
- Do not edit or commit generated `docs/` pages for this workflow. Actions
  deploys the staged artifact.
- Full local build:

```bash
export PHISPACE_SITE_DIR=../Test/site-vignettes-next
export PHISPACE_STEREOSEQ_INPUTS=../Test/stereoseq-inputs
export PHISPACE_VISIUM_INPUTS=../Test/visium-inputs
export PHISPACE_GETTING_STARTED_INPUTS=../Test/getting-started-inputs
export PHISPACE_BRIDGE_ANNOTATION_INPUTS=../Test/bridge-annotation-inputs
export PHISPACE_CITE_SEQ_INPUTS=../Test/cite-seq-inputs
export PHISPACE_PERTURBSEQ_INPUTS=../Test/perturbseq-inputs
export PHISPACE_COSMX_INPUTS=../Test/cosmx-inputs
Rscript --vanilla scripts/build-stereoseq.R
Rscript --vanilla scripts/build-visium.R
Rscript --vanilla scripts/build-bridge-annotation.R
Rscript --vanilla scripts/build-getting-started.R
Rscript --vanilla scripts/build-cite-seq.R
Rscript --vanilla scripts/build-perturbseq.R
Rscript --vanilla scripts/build-cosmx.R
Rscript --vanilla scripts/build-reference.R
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
- Sorted k-means cluster sizes use a tolerance of 1% of each size. The first
  Actions run with result checks (run `37538262677`, 2026-10-07) failed on
  Stereo-seq PhiSpace niche sizes: the runner's BLAS moved 1 to 13 of 15,454
  bins between niches, while all continuous metrics, including the k-means
  within-cluster sum of squares, passed at 1e-5. Exact tolerances on counts
  are safe only for quantities that do not depend on floating-point ties.
- Stereo-seq difference from the established Dropbox results: the vignette
  removes the 1% of bins with the lowest total counts before annotation (15,611
  to 15,454 bins), but the established `output/PhiRes.qs2` was computed on all
  15,611 bins. Rerunning on all bins reproduces it (max absolute difference
  3.3e-6). On the shared bins, fresh scores differ only by a per-column
  constant (SD of differences <= 1.4e-6) from query centering, but the
  max-absolute scaling in `normPhiScores()` changes. PhiSpace niche cluster
  sizes therefore differ by up to 17 bins; barcode and gene-expression
  clusterings and enrichment scores match. The expected values follow the
  vignette code.
- Stereo-seq cached results: with current ggplot2 (4.0.3), the established
  `output/cloneKDEres.qs2` renders the 16-clone density figure blank without an
  error, because it stores ggplot objects saved before ggplot2 4.0. Replacement
  `PhiRes.qs2`, `PhiClustRes.qs2` and `cloneKDEres.qs2`, computed with the
  vignette code, are in `PkgOverhaul/data/StereoSeq/replacement-2026-10-06/`
  with hashes and a README. With them and the DWD seed below, a cached build
  passes all 49 result checks and reproduces all 23 shared fresh-build figures
  byte for byte. Pending: the user uploads them to the Dropbox `output/` folder.
- Stereo-seq DWD figure: `cv.kerndwd()` assigns random folds. Without a seed,
  fresh builds depended on earlier `set.seed()` calls inside compute branches,
  which cached builds skip; the selected lambda and loadings then changed (top
  loadings and their scale differed about 2-fold). The vignette now calls
  `set.seed(94871)` before `cv.kerndwd()`. This changed the published DWD
  loading figure once; its top positive loadings are HPC(BM) and
  Neutro(Spleen). The loadings are sensitive to the fold assignment, so the
  figure should be read as illustrative.
- Visium fresh results match the established `combo_PhiRes.qs2` within 1.1e-6
  relative difference on all checked metrics.
- BridgeAnnotation archive (RDS, no qs conversion):
  `https://www.dropbox.com/scl/fo/jeuqzjfyyr2j7doa922ve/AI-2U_wtZBpPGOswxMzReGQ?rlkey=n328yyr2llf81gz3chynjg0r6&dl=1`
  with 8 inputs and 3 cached results in `scripts/bridge-annotation-inputs.tsv`.
  The vignette read `data/CellTypeTable.rds`, but the file is
  `cellTypeTable.rds`; this only worked on case-insensitive file systems.
  Fresh builds: about 5 minutes 18 seconds, 24.4 and 17.3 GB peak RAM in two runs, 28 result metrics within
  1e-5 of the established results, and the printed classification errors
  identical to the published page (GA 0.2291131/0.3536215, peaks
  0.2892031/0.3978839). Cached build: 37 seconds, 5.4 GB. `mvr()` converts
  the predictor matrix to dense (since at least `afcf7c8`, 2025-01), so
  `center = F` does not keep the peak matrix sparse; with Matrix 1.7.5 the
  fresh page shows 2.6, 8.6 and 5.4 GiB coercion warnings.
- CosMx folder:
  `https://www.dropbox.com/scl/fo/z01qvst71kzxqog4bf5po/AKZzi-PgkZCDYZdUwLqfBWg?rlkey=65kg8zqatk2zdk5mnirgn1vmc&dl=1`
  with 14 inputs and 2 cached results in `scripts/cosmx-inputs.tsv`. All 15
  `.qs` objects were converted to qs2 (restored objects identical).
  `output/CosMxAllLungsPhiRes4Refs.qs2` (all eight lungs, four references) is
  loaded unconditionally and no vignette code computes it, so it is an input;
  builds verify only its hash. The DWD step sets its seed immediately before
  `cv.kerndwd()`, so fresh and cached builds agree. Fresh build: 1 minute 54
  seconds, 9.7 GB; cached: 1 minute 10 seconds, 9.5 GB; 35 result metrics pass
  in both. Fresh Lung5_Rep1 results match the established ones to about 1e-9
  with identical niche sizes; printed outputs and the tumour-signature table
  match the published page; 13 of 14 figures are byte-identical between fresh
  and cached builds (the mesothelial map differs in 0.024% of pixels).
- PerturbSeq folder:
  `https://www.dropbox.com/scl/fo/4gm1ef27yl4wgb0f6zio7/AG0bvE9LXKgv61AAKmD5ho8?rlkey=bkd7nlwo6qlltkl21645xdw90&dl=1`
  with 3 inputs and 1 cached result in `scripts/perturbseq-inputs.tsv`. The
  three `.qs` objects were converted to qs2 (restored objects identical). The
  public `sceStim.qs` link served the same file as the private Dropbox copy.
  `utils.R` (defines `tempPvals()`, which uses the vignette's `query` and
  `ctrlReg`) existed only in the private folder. Rebuilding the DICE reference
  with celldex 1.20.0 and the vignette code reproduces `ref_sce` (assays within
  4e-12, identical labels), so `ref_sce.qs2` is shipped as an input and Actions
  needs no celldex or ExperimentHub download; the celldex code runs only when
  the file is missing. Fresh build: 59 seconds, 13.4 GB; cached: 41 seconds,
  6.5 GB; 22 result metrics pass in both. Printed results match the published
  page; all 11 figures are byte-identical between fresh and cached builds, and
  differ from the published page only in fonts.
- CITE-seq archive (RDS, no qs conversion):
  `https://www.dropbox.com/scl/fo/it8uwxd2v4k2lyoo936at/AKws5m_xeNjOOs8Vkho-wP8?rlkey=s4d5a7cfkl5sj9qiogc7mng5t&dl=1`
  with 7 inputs and 5 cached results in `scripts/cite-seq-inputs.tsv`. The
  fresh RNA branch saved `PhiRes` but later code reads `PhiResRNA`; fixed in
  `13ca910`. The ADT and RNA objects share cells, order and metadata, so the
  fresh branch replacing `reference`/`query` with the RNA objects is harmless.
  Fresh build: 7 minutes 42 seconds, 64.9 GB peak RAM (dense coercions of 10.1
  and 15.2 GiB), 34 result metrics within 1e-5 of the established results.
  Cached build: 52 seconds, 3.5 GB. Printed results match the published page;
  heatmaps match the cached build to 1/255 colour value; the UMAP grid differs
  in a fresh build because UMAP is recomputed (checks cover only UMAP cells).
- Getting Started: the three Dropbox `.qs` inputs were converted to `.qs2`
  (restored objects `identical()` to the originals) and republished in one
  shared folder; hashes are in `scripts/getting-started-inputs.tsv`. All chunks
  passed fresh (1 minute 10 seconds, 11.7 GB peak RAM) and cached (12 seconds)
  without loading `qs`. All five printed result blocks, including the 1782
  selected genes, match the previously published page. The CV tuning chunk
  remains disabled with `if(F)`, so builds do not test it. A Matrix
  sparse-to-dense coercion warning now appears in the rendered page.
  The old `tuneRes.rds` file was a qs file with an `.rds` name; it is now
  `tuneRes.qs2`.

Relevant commits:

- `a1080bc` — fix the Stereo-seq assay argument.
- `a2d20e6` — add Stereo-seq Actions publishing.
- `d61c249` — migrate Stereo-seq to qs2.
- `797648c` — migrate and automate Visium and combine publishing.

Handoff checkpoint on 2026-10-07: commit `c77f1f1` was pushed; Actions run
`37580149580` built all seven vignettes (Stereo-seq, Visium, Getting Started,
PerturbSeq and CosMx fresh; BridgeAnnotation and CITE-seq from cached results)
with result checks in about 10 minutes, and deployed. The live CosMx article
was checked for its new content and contains no private Dropbox paths. Job
logs need repository admin rights; the public API gives run and step status.
Future sessions should still recheck current Actions state.

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
  Since 2026-10-08, 40 tests (137 expectations) also cover `findNiches()`,
  `getPC()`, `zeroFeatQC()`, `spatialSampler()`, `spatialSmoother()` and
  `clusterPhiSpace()`, using the simulated 300-cell SpatialExperiment in
  `helper-spatial.R`. Fifteen of the new tests fail on the code before that
  day's fixes (commit `a5dd5e1`). `clusterPhiSpace()` keeps `ncol(x) - 1`
  components; it used irlba, which warned and returned wrong trailing
  components (scores up to 3.4 from `prcomp()` with 6 columns). It now uses
  the internal full-SVD `.getPC_svd()` (same elements as `getPC()`), which
  matches `prcomp()` to 1e-14; 42 tests, 152 expectations, no warnings.
- **Integration tests**: none are currently available. Earlier notes described
  `Test/test_cellTypeThreshold.R` with CosMx lung data, but on 2026-10-06 neither
  the script nor its data existed. The vignette builds and their result checks
  act as integration tests.

## Known Check Notes

- None. On 2026-10-08, `R CMD check --no-manual --ignore-vignettes` (R 4.5.0)
  gave Status OK after `plot.PhiSpaceClustering` switched to `.data$` columns.
  The workflow `.github/workflows/R-CMD-check.yaml` runs the same check on
  pushes to `main` and pull requests and fails on warnings. Vignettes are
  excluded because they need external data; the Pages workflow builds them.
