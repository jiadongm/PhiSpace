# PhiSpace: Comprehensive Usage Guide

> **Purpose**: Reference document for Claude Code sessions working with PhiSpace and its sister package vizOmics.
> **Package source**: `pkg/` directory. Install with `devtools::install("pkg")`.
> **Documentation site**: https://jiadongm.github.io/PhiSpace/
> **Workshop**: https://bhuvad.github.io/CosMxSpatialAnalysisWorkshop/articles/Part_II_PhiSpace.html

---

## 1. Package Overview

PhiSpace annotates query single-cell or spatial omics data against an annotated reference dataset using **partial least squares (PLS) regression**. Instead of discrete cell type labels, it produces **continuous phenotype scores** — one score per phenotype per cell — enabling graded annotation across a phenotypic continuum.

**Key concepts:**
- **Reference**: An annotated `SingleCellExperiment` (SCE) with categorical phenotypes in `colData()` (e.g., cell type, disease state).
- **Query**: An SCE or `SpatialExperiment` (SPE) to be annotated.
- **PhiSpace scores**: A matrix stored in `reducedDim(query, "PhiSpace")` with dimensions cells × phenotype-levels. Values after normalization range in [-1, 1] (column-wise) or [0, 1] (row-wise).
- **No batch correction needed**: Reference and query only need to use the same normalization method (e.g., both rank-transformed).

---

## 2. Installation

```r
# Install from GitHub (Bioconductor dependencies required)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("jiadongm/PhiSpace/pkg")

# Load packages
library(PhiSpace)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
```

---

## 3. Core Workflow

### 3.1 Data Preparation

Reference and query must share gene names and use the **same normalization** in the assay passed to PhiSpace.

```r
# --- Rank transformation (recommended for cross-platform) ---
# Stored as assay "rank"; use refAssay = "rank"
reference <- RankTransf(reference, assayName = "counts")
query     <- RankTransf(query,     assayName = "counts")

# --- Log1p normalization ---
# Stored as assay "log1p"; use refAssay = "log1p"
reference <- logTransf(reference, use_log1p = TRUE, targetAssay = "log1p")
query     <- logTransf(query,     use_log1p = TRUE, targetAssay = "log1p")

# --- scran normalization (for single-cell reference) ---
reference <- scranTransf(reference)   # creates "logcounts" assay

# --- Quality control ---
query <- zeroFeatQC(query)            # removes all-zero genes

# --- Stratified subsampling (to balance cell types in reference) ---
reference <- subsample(
  reference,
  key        = "cell_type",   # colData column to stratify by
  proportion = 0.2,           # fraction to keep
  minCellNum = 50,            # minimum cells per category
  seed       = 5202056
)
```

### 3.2 Basic Annotation with `PhiSpace()`

`PhiSpace()` is the **main user-facing wrapper**. It handles single or multiple references, stores normalized scores in `reducedDim(query, "PhiSpace")`, and returns the updated query object.

```r
# Single reference, single query
query <- PhiSpace(
  reference  = reference_sce,
  query      = query_sce,
  phenotypes = "cell_type",      # column name(s) in colData(reference)
  refAssay   = "logcounts",      # assay in reference
  queryAssay = "logcounts",      # assay in query (defaults to refAssay)
  regMethod  = "PLS",            # "PLS" (default) or "PCA"
  ncomp      = NULL,             # NULL → auto = number of phenotype levels
  center     = TRUE,             # recommended TRUE
  scale      = FALSE             # recommended FALSE (see §10)
)

# Access scores: matrix of shape cells × cell_types
scores <- reducedDim(query, "PhiSpace")
dim(scores)       # e.g., 5000 × 20
colnames(scores)  # cell type names

# Discretize to single label per cell
query$PhiCellType <- getClass(reducedDim(query, "PhiSpace"))
```

**Multiple phenotype columns** (e.g., cell type + disease):

```r
query <- PhiSpace(
  reference  = reference_sce,
  query      = query_sce,
  phenotypes = c("cell_type", "disease"),  # two colData columns
  refAssay   = "rank"
  # ncomp auto = n_cell_types + n_disease_levels
)
```

**Multiple references** (scores are column-concatenated, with `"(RefName)"` suffix):

```r
ref_list <- list(
  immune       = immune_sce,
  epithelial   = epithelial_sce,
  endothelial  = endothelial_sce
)

query <- PhiSpace(
  reference  = ref_list,
  query      = query_sce,
  phenotypes = "cell_type",
  refAssay   = "log1p"
)
# Column names: "CD4(immune)", "AT1(epithelial)", "EC(endothelial)", ...
```

**Multiple queries:**

```r
query_list <- PhiSpace(
  reference  = reference_sce,
  query      = list(query1_sce, query2_sce),
  phenotypes = "cell_type",
  refAssay   = "logcounts"
)
# Returns a list of updated query objects
```

**Retrieve reference scores** (`updateRef = TRUE`):

```r
result <- PhiSpace(
  reference  = reference_sce,
  query      = query_sce,
  phenotypes = "cell_type",
  refAssay   = "logcounts",
  updateRef  = TRUE          # returns list(reference=..., query=...)
)
reference_updated <- result$reference
query_updated     <- result$query

YrefHatNorm <- reducedDim(reference_updated, "PhiSpace")
```

**Store raw (unnormalized) scores** alongside normalized:

```r
query <- PhiSpace(..., storeUnNorm = TRUE)
raw_scores  <- reducedDim(query, "PhiSpace_nonNorm")
norm_scores <- reducedDim(query, "PhiSpace")
```

### 3.3 Single-Reference Core: `PhiSpaceR_1ref()`

Use this when you need the full return value (importance scores, model, reference predictions).

```r
PhiRes <- PhiSpaceR_1ref(
  reference  = reference_sce,
  query      = query_sce,       # or list of SCEs
  phenotypes = "cell_type",
  refAssay   = "log1p",         # default "log1p" (not "rank" like PhiSpace wrapper)
  regMethod  = "PLS",
  ncomp      = NULL,
  center     = TRUE,
  scale      = FALSE
)

# Return value fields:
PhiRes$PhiSpaceScore   # raw query scores (matrix or list of matrices)
PhiRes$PhiSpaceNorm    # normalized query scores
PhiRes$YrefHat         # raw reference predictions
PhiRes$YrefHatNorm     # normalized reference predictions
PhiRes$impScores       # feature importance matrix (genes × phenotypes)
PhiRes$selectedFeat    # genes used in model
PhiRes$ncomp           # number of components
PhiRes$phenoDict       # data.frame mapping column names to phenotype categories
PhiRes$atlas_re        # internal model object (SuperPC output)
```

**Manually normalize scores:**

```r
sc_norm <- normPhiScores(PhiRes$PhiSpaceScore, method = "col")
```

### 3.4 Feature Selection

Feature selection can improve model accuracy and interpretability. Two approaches:

```r
# Approach 1: let PhiSpace select nfeat top genes per phenotype
query <- PhiSpace(
  reference  = reference_sce,
  query      = query_sce,
  phenotypes = "cell_type",
  refAssay   = "logcounts",
  nfeat      = 200        # union of top 200 genes per cell type
)

# Approach 2: provide pre-selected genes
query <- PhiSpace(
  reference     = reference_sce,
  query         = query_sce,
  phenotypes    = "cell_type",
  refAssay      = "logcounts",
  selectedFeat  = my_gene_list  # character vector of gene names
)
```

### 3.5 Filtering Rare Cell Types: `cellTypeThreshold`

```r
query <- PhiSpace(
  reference          = reference_sce,
  query              = query_sce,
  phenotypes         = "cell_type",
  refAssay           = "logcounts",
  cellTypeThreshold  = 50    # remove cell types with < 50 reference cells
)
```

### 3.6 Parameter Tuning: `tunePhiSpace()`

Cross-validates `ncomp` and/or `nfeat` on the reference.

```r
tune_res <- tunePhiSpace(
  reference     = reference_sce,
  assayName     = "logcounts",
  phenotypes    = "cell_type",
  tune_ncomp    = TRUE,
  tune_nfeat    = TRUE,
  ncompLimits   = c(1, 100),
  ncompGridSize = 20,
  nfeatLimits   = c(10, 15000),
  nfeatGridSize = 30,
  Kfolds        = 5,
  regMethod     = "PLS",
  Ncores        = 4,
  seed          = 5202056,
  cellTypeThreshold = 50
)

# Use tuned parameters
query <- PhiSpace(
  reference    = reference_sce,
  query        = query_sce,
  phenotypes   = "cell_type",
  refAssay     = "logcounts",
  ncomp        = tune_res$ncomp,
  selectedFeat = tune_res$selectedFeat
)

# Inspect tuning curves
tune_res$tuneNcompErrs   # CV error by ncomp
tune_res$tuneNfeatErrs   # CV error by nfeat
tune_res$impScores       # feature importance matrix
```

### 3.7 Direct Response Matrix

For continuous or custom response variables:

```r
# Build response matrix manually (cells × phenotypes)
Y <- matrix(...)
rownames(Y) <- colnames(reference_sce)
colnames(Y) <- c("age", "BMI", "treatment_score")

PhiRes <- PhiSpaceR_1ref(
  reference = reference_sce,
  query     = query_sce,
  response  = Y,           # overrides phenotypes
  refAssay  = "logcounts"
)
```

---

## 4. Score Normalization: `normPhiScores()`

```r
normPhiScores(X, method = c("col", "row"))
```

| Method | Formula | Range | Use case |
|--------|---------|-------|---------|
| `"col"` (default) | `(x - median(x)) / max(|x - median(x)|)` per column | [-1, 1] | Compare across cell types; default in `PhiSpaceR_1ref()` |
| `"row"` | `((x - min) / (max - min))^3` per row | [0, 1] | SingleR-style; emphasize dominant cell type per cell |

```r
# Column-wise (default): stored in reducedDim(query, "PhiSpace")
sc_col <- normPhiScores(raw_scores, method = "col")

# Row-wise: SingleR-compatible
sc_row <- normPhiScores(raw_scores, method = "row")
```

**PhiSpace()** stores column-wise normalized scores in `reducedDim`. `PhiSpaceR_1ref()` returns both `$PhiSpaceNorm` (col) and `$PhiSpaceScore` (raw).

---

## 5. Spatial Analysis

### 5.1 Spatial Visualization: `VizSpatial()`

Plots spatial distribution of metadata, gene expression, or PhiSpace scores. Returns a `ggplot2` object.

```r
VizSpatial(
  obj,
  x_coord        = "x",          # colData column for x (SCE only; SPE auto-detects)
  y_coord        = "y",          # colData column for y
  ptSize         = 2,
  ptShape        = 16,
  colBy          = NULL,         # color by a colData column (discrete or continuous)
  feature        = NULL,         # color by gene expression
  assay2use      = NULL,         # assay for feature (default: "counts")
  reducedDim     = NULL,         # name of column in reducedDim slot to color by
  reducedDim2use = "PhiSpace",   # which reducedDim slot
  censor         = FALSE,        # set low values to quantile (recommended for PhiSpace scores)
  quant          = 0.5,          # censoring quantile
  legend.position = "right",
  reOrder        = FALSE,        # reorder points by value (avoids overplotting)
  fsize          = 14
)
```

**Examples:**

```r
# Color by discrete cell type label
VizSpatial(query, colBy = "PhiCellType", ptSize = 1)

# Color by specific PhiSpace score (e.g., "AT1" cell type)
VizSpatial(query, reducedDim = "AT1", reducedDim2use = "PhiSpace",
           censor = TRUE, quant = 0.5, ptSize = 1)

# Color by gene expression
VizSpatial(query, feature = "EPCAM", assay2use = "logcounts", ptSize = 1)

# For SPE: coordinates auto-detected from spatialCoords()
VizSpatial(spe, reducedDim = "CD4_T", censor = TRUE, ptSize = 0.5)
```

**Note**: For `SpatialExperiment` objects, `x_coord`/`y_coord` arguments are ignored — coordinates are read from `spatialCoords(obj)`. For `SingleCellExperiment`, specify `x_coord` and `y_coord` as `colData` column names.

### 5.2 Batch Spatial Export: `saveCellTypeMaps()`

Creates multi-panel heatmaps for all (or selected) cell types and saves to disk.

```r
saveCellTypeMaps(
  sce,
  methodName     = "PhiSpace",                       # reducedDim slot
  tissueName     = "SpatialSample",                  # used in filenames
  coordNames     = c("x", "y"),                      # auto-detected for SPE
  freeColScale   = FALSE,   # FALSE: unified scale (good for comparison)
  quants         = c(0.1, 1),                        # color limits (unified scale)
  psize          = 0.5,
  savePlots      = TRUE,
  returnPlots    = FALSE,   # TRUE: return list of ggplot objects
  legendPosition = "none",
  censQuant      = 1,       # censor scores < this quantile (1 = no censoring)
  ctypes         = NULL,    # NULL = all; or character vector of cell type names
  outputDir      = "figs/spatialCellTypeHeatmaps",
  width          = 10, height = 10,
  fignrow        = 4, figncol = 4,    # grid per output file
  fileFormat     = "png",             # "png", "pdf", "jpeg", "tiff"
  dpi            = 300
)
```

**Examples:**

```r
# Save all cell types with censoring
saveCellTypeMaps(query, tissueName = "Lung5", censQuant = 0.9)

# Return plots for further customization
plots <- saveCellTypeMaps(query, savePlots = FALSE, returnPlots = TRUE)
plots[["AT1"]]  # access individual plot

# Specific cell types with free color scales
saveCellTypeMaps(query,
  ctypes       = c("AT1", "AT2", "Club"),
  freeColScale = TRUE,
  fignrow = 1, figncol = 3)
```

**Output structure**: Files saved as `outputDir/tissueName/tissueName_PhiSpace_1.png`, `..._2.png`, etc. (one file per `fignrow × figncol` panel of cell types).

### 5.3 Spatial Smoothing: `spatialSmoother()`

KNN-based spatial smoothing of PhiSpace scores or gene expression.

```r
# Smooth PhiSpace scores (default)
query <- spatialSmoother(
  object            = query,
  reducedDim2smooth = "PhiSpace",
  smoothReducedDim  = TRUE,          # TRUE: smooth reducedDim; FALSE: smooth assay
  smoothedReducedDim = NULL,         # auto: "PhiSpace_smoothed"
  k                 = 10,            # number of spatial neighbors
  kernel            = "gaussian",    # "gaussian", "uniform", or "linear"
  sigma             = NULL,          # auto-computed from median distance
  include_self      = TRUE,
  verbose           = TRUE
)
# Result stored in reducedDim(query, "PhiSpace_smoothed")

# For SCE objects, specify coordinate columns
query <- spatialSmoother(
  query,
  x_coord = "sdimx", y_coord = "sdimy",
  k = 15, kernel = "gaussian"
)

# Smooth gene expression instead
query <- spatialSmoother(
  query,
  smoothReducedDim = FALSE,
  assay2smooth     = "logcounts",
  smoothedAssay    = "logcounts_smoothed",
  k = 10
)
```

**Kernels:**
- `"gaussian"`: `exp(-d² / (2σ²))`, σ auto-set to median distance
- `"uniform"`: equal weights for all k neighbors
- `"linear"`: `1 - d/max_d`

### 5.4 Spatial Niche Identification: `findNiches()`

K-means clustering on PhiSpace scores (with optional PCA pre-processing). Simpler interface than `clusterPhiSpace()`.

```r
# Single k
spe <- findNiches(
  spe,
  reducedDimName  = "PhiSpace",
  n_niches        = 9,
  use_pca         = TRUE,
  ncomp           = NULL,   # auto: min(30, floor(n_cell_types / 2))
  center          = TRUE, scale = FALSE,
  kmeans_iter     = 500, kmeans_nstart = 20,
  seed            = 94863,
  store_pca       = FALSE,
  pca_name        = "PhiSpace_PCA"
)
# Stored as colData(spe)$spatial_niches (factor: "Niche_1", ..., "Niche_9")

# Multiple k values (for comparison)
spe <- findNiches(spe, n_niches = c(7, 8, 9, 10))
# Stored as colData(spe)$spatial_niches_k7, ...$spatial_niches_k10

# Visualize
VizSpatial(spe, colBy = "spatial_niches", ptSize = 0.5)
table(spe$spatial_niches)
```

### 5.5 Clustering with Diagnostics: `clusterPhiSpace()`

More flexible than `findNiches()`. Returns a `PhiSpaceClustering` S3 object with print/summary/plot methods.

```r
# Fixed k
result <- clusterPhiSpace(
  x              = query,           # SCE/SPE, matrix, data.frame, or sparse matrix
  k              = 9,
  reducedDimName = "PhiSpace",      # ignored if x is a matrix
  ncomp          = NULL,            # PCs to use; auto = min(30, ncol-1)
  nstart         = 20,
  iter.max       = 500,
  algorithm      = "Lloyd",
  center         = TRUE, scale = FALSE,
  seed           = 123,
  return_pca     = TRUE,
  store_in_colData = FALSE,
  cluster_name   = "PhiClust"
)

# Automatic k selection
result <- clusterPhiSpace(
  x              = query,
  k_range        = c(5, 15),
  select_k_method = "silhouette"   # or "elbow"
)

# S3 methods
print(result)            # cluster sizes, quality metrics
summary(result)          # full details + CV metrics
plot(result, type = "pca")        # PC1 vs PC2 colored by cluster
plot(result, type = "elbow")      # elbow curve (k_range only)
plot(result, type = "silhouette") # silhouette curve (k_range only)
plot(result, type = "variance")   # PCA variance explained

# Extract assignments
query$niche <- result$clusters
result$optimal_k         # selected k
result$pca_result        # PCA output (getPC() format)
result$pc_scores         # matrix used for clustering
result$kmeans_result     # raw kmeans object
```

**Cluster on assay instead of PhiSpace:**

```r
result <- clusterPhiSpace(
  x        = query,
  use_assay = "logcounts",
  k        = 9,
  ncomp    = 50
)
```

---

## 6. Feature Analysis

### 6.1 Feature Ranking: `rankFeatures()`

Rank genes or PhiSpace score dimensions by discriminative power.

```r
result <- rankFeatures(
  data           = query,           # matrix, SCE/SPE/SummarizedExperiment
  response       = "spatial_niches", # colData column OR vector
  method         = "PLSDA",         # "PLSDA" (multi-class), "PLS" (continuous), "DWD" (binary)
  source         = "reducedDim",    # "reducedDim" or "assay"
  assay_name     = "logcounts",     # used when source = "assay"
  reducedDim_name = "PhiSpace",     # used when source = "reducedDim"
  ncomp          = NULL,            # auto: n_classes for PLSDA
  center         = TRUE, scale = FALSE,
  seed           = NULL
)

# Results
head(result$feature_ranking, 20)   # data.frame: feature, importance_norm, importance_<class>...
result$importance_scores            # matrix: features × classes
result$model                        # fitted PLS/DWD model
result$response_summary             # distribution of response
```

**Method summary:**
- `"PLSDA"`: Multi-class; ranked by Euclidean norm of coefficients across all classes
- `"PLS"`: Continuous response; ranked by `|coefficient|`
- `"DWD"`: Binary only; requires `kerndwd` package; ranked by discriminant weight magnitude

**Ranking genes in expression assay:**

```r
result <- rankFeatures(
  data      = query,
  response  = "niche",
  method    = "PLSDA",
  source    = "assay",
  assay_name = "logcounts"
)
```

### 6.2 Cell Type Correlation: `cellTypeCorMat()`

S4 generic computing pairwise correlation between cell type scores.

```r
# On a matrix
corMat <- cellTypeCorMat(reducedDim(query, "PhiSpace"))

# On an SCE/SPE (uses reducedDim)
corMat <- cellTypeCorMat(query)
```

### 6.3 PCA on PhiSpace Scores: `getPC()`

General-purpose PCA using partial SVD (irlba). Used internally by `clusterPhiSpace()` and `findNiches()`.

```r
pc_res <- getPC(
  X      = reducedDim(query, "PhiSpace"),
  ncomp  = 20,
  center = TRUE,
  scale  = FALSE
)

pc_res$scores    # cells × ncomp matrix, row names from X
pc_res$loadings  # features × ncomp, row names = colnames(X)
pc_res$props     # proportion variance explained per PC
pc_res$accuProps # cumulative variance explained
pc_res$sdev      # standard deviations
pc_res$Xmeans   # column means used for centering (NULL if center=FALSE)

# Project query using reference loadings (manual transfer)
queryEmbedding <- scale(reducedDim(query, "PhiSpace"),
                        center = TRUE, scale = FALSE) %*% pc_res$loadings
```

---

## 7. Utility Functions

### Data Encoding

```r
# Convert colData phenotypes to dummy matrix (cells × phenotype_levels)
YY <- codeY(sce, phenotypes = c("cell_type", "disease"), method = "-1,1")
# method: "-1,1" (default) or "0,1"

# Convert factor/vector to dummy matrix
YY <- codeY_vec(vec, rowNames = NULL, method = "-1,1")
```

### Data Preprocessing

```r
# Rank transformation → "rank" assay
sce <- RankTransf(sce, assayName = "counts", sparse = FALSE)

# Log1p transformation → target assay
sce <- logTransf(sce, use_log1p = TRUE, targetAssay = "log1p")

# scran normalization → "logcounts" assay
sce <- scranTransf(sce)

# CLR normalization (e.g., for CITE-seq protein)
sce <- CLRnorm(sce)

# Remove all-zero genes
sce <- zeroFeatQC(sce)

# Keep only common genes between reference and query
sce_list <- KeepCommonGenes(list(reference, query))  # or two separate objects
```

### Dimensionality Reduction

```r
# UMAP on PhiSpace scores or assay
query <- computeUMAP(query, reducedDimName = "PhiSpace")
# Stores result in reducedDim(query, "UMAP")
```

### Cell Aggregation & Sampling

```r
# Pseudobulk aggregation
pb_sce <- pseudoBulk(sce, group_by = "sample_id")

# Stratified subsampling
sce_small <- subsample(sce, key = "cell_type", proportion = 0.2,
                       minCellNum = 50, seed = 5202056)
```

### Score Utilities

```r
# Discretize continuous scores to single label
labels <- getClass(reducedDim(query, "PhiSpace"), labPerSample = 1)
# labPerSample > 1: return top-N cell types per cell (matrix output)

# Convert scores to probabilities (softmax-like)
probs <- Score2Prob(reducedDim(query, "PhiSpace"))

# Translate between label schemes
new_labels <- translateLabel(old_labels, mapping_df)
```

### Visualization Heatmap

```r
# Heatmap of PhiSpace scores grouped by query cell type
library(ComplexHeatmap)
p <- plotPhiSpaceHeatMap(
  PhiSpaceScore = reducedDim(query, "PhiSpace"),
  phenoDict     = PhiRes$phenoDict,  # from PhiSpaceR_1ref()
  queryLabs     = query$cell_type,
  name          = "PhiSpace score",
  column_names_rot = 20,
  show_row_dend    = FALSE,
  show_column_dend = TRUE
)
draw(p, heatmap_legend_side = "top")
```

---

## 8. Complete End-to-End Examples

### Example A: scRNA-seq Reference → scRNA-seq Query

```r
library(PhiSpace)
library(SingleCellExperiment)

# 1. Prepare data
reference <- scranTransf(reference)   # logcounts
query     <- scranTransf(query)

# 2. Subsample reference to balance cell types
reference <- subsample(reference, key = "cell_type", proportion = 0.2, minCellNum = 50)

# 3. Run PhiSpace with reference + query scores
result <- PhiSpace(
  reference  = reference,
  query      = query,
  phenotypes = c("cell_type", "sample_source"),
  refAssay   = "logcounts",
  updateRef  = TRUE
)
reference <- result$reference
query     <- result$query

# 4. Access scores
query_scores <- reducedDim(query, "PhiSpace")
ref_scores   <- reducedDim(reference, "PhiSpace")

# 5. PCA on reference scores for 3D visualization
pc_res <- getPC(ref_scores, ncomp = 3)
ref_embedding   <- as.data.frame(pc_res$scores)
query_embedding <- as.data.frame(
  scale(query_scores, center = TRUE, scale = FALSE) %*% pc_res$loadings
)

# 6. Heatmap of query scores
query$PhiCellType <- getClass(query_scores)
```

### Example B: CosMx Spatial Transcriptomics

```r
library(PhiSpace)
library(SpatialExperiment)

# 1. Prepare data
spe <- logTransf(spe, use_log1p = TRUE, targetAssay = "log1p")
spe <- zeroFeatQC(spe)

ref_list <- list(
  immune      = immune_sce,
  epithelial  = epithelial_sce,
  stromal     = stromal_sce
)
ref_list <- lapply(ref_list, logTransf, use_log1p = TRUE, targetAssay = "log1p")

# 2. Multi-reference annotation
spe <- PhiSpace(
  reference  = ref_list,
  query      = spe,
  phenotypes = "cell_type",
  refAssay   = "log1p"
)
# reducedDim columns: "AT1(epithelial)", "CD4(immune)", "Fibroblast(stromal)", ...

# 3. Smooth scores
spe <- spatialSmoother(spe, k = 10, kernel = "gaussian")

# 4. Discretize
spe$PhiCellType <- getClass(reducedDim(spe, "PhiSpace"))

# 5. Visualize individual cell type
VizSpatial(spe, reducedDim = "AT1(epithelial)", censor = TRUE, quant = 0.5, ptSize = 0.5)

# 6. Save all cell type maps
saveCellTypeMaps(spe, tissueName = "Lung5_Rep1", censQuant = 0.9, psize = 0.3)

# 7. Identify spatial niches
spe <- findNiches(spe, n_niches = 9, use_pca = TRUE)
VizSpatial(spe, colBy = "spatial_niches", ptSize = 0.5)

# 8. Rank cell types defining each niche
result <- rankFeatures(spe, response = "spatial_niches", method = "PLSDA",
                       source = "reducedDim", reducedDim_name = "PhiSpace")
head(result$feature_ranking, 10)
```

### Example C: Workshop-Style Simplified Workflow

```r
library(PhiSpace)

# Built-in example data
data("lung5_norm")   # pre-normalized CosMx SPE
data("ref_luca")     # reference scRNA-seq

# Normalize reference
ref_luca <- scranTransf(ref_luca)

# Annotate
lung5_norm <- PhiSpace(
  reference  = ref_luca,
  query      = lung5_norm,
  phenotypes = "cell_type",
  refAssay   = "logcounts"
)

# View scores
scores <- reducedDim(lung5_norm, "PhiSpace")
scores[1:5, 1:5]

# Discrete label
lung5_norm$PhiCellType <- getClass(scores)

# Interactive visualization
library(plotly)
p <- VizSpatial(lung5_norm, colBy = "PhiCellType", ptSize = 1)
plotly::ggplotly(p)

# Niche analysis
lung5_norm <- findNiches(lung5_norm, n_niches = 8)
VizSpatial(lung5_norm, colBy = "spatial_niches", ptSize = 0.8)
```

---

## 9. Key Parameters Quick Reference

### `PhiSpace()` / `PhiSpaceR_1ref()`

| Parameter | Default | Notes |
|-----------|---------|-------|
| `reference` | required | SCE or named list of SCEs |
| `query` | required | SCE, SPE, or list |
| `phenotypes` | `NULL` | colData column names; mutually exclusive with `response` |
| `response` | `NULL` | Direct response matrix (rows=cells, cols=phenotypes); for continuous phenotypes |
| `refAssay` | `"rank"` (PhiSpace) / `"log1p"` (PhiSpaceR_1ref) | Must match `queryAssay` normalization |
| `queryAssay` | `NULL` → `refAssay` | |
| `regMethod` | `"PLS"` | `"PLS"` or `"PCA"` |
| `ncomp` | `NULL` → n_phenotype_levels | Number of PLS components |
| `nfeat` | `NULL` → all | Top genes per phenotype (union used) |
| `selectedFeat` | `NULL` → all | Pre-selected gene names; overrides `nfeat` |
| `center` | `TRUE` | Recommended TRUE |
| `scale` | `FALSE` | Recommended FALSE with many features |
| `cellTypeThreshold` | `NULL` | Remove cell types with fewer than N cells |
| `updateRef` | `FALSE` | If TRUE, returns `list(reference=..., query=...)` |
| `storeUnNorm` | `FALSE` | Store raw scores in `"PhiSpace_nonNorm"` |
| `reducedDimName` | `"PhiSpace"` | Where to store scores (PhiSpace wrapper only) |

### `tunePhiSpace()`

| Parameter | Default | Notes |
|-----------|---------|-------|
| `reference` | required | Reference SCE |
| `assayName` | `"logcounts"` | Assay to use |
| `phenotypes` | `NULL` | colData columns |
| `tune_ncomp` | `TRUE` | Cross-validate ncomp |
| `tune_nfeat` | `TRUE` | Cross-validate nfeat |
| `ncompLimits` | `c(1, 100)` | ncomp search range |
| `nfeatLimits` | `c(10, 15000)` | nfeat search range |
| `Kfolds` | `5` | CV folds |
| `Ncores` | `1` | Parallel cores |
| `cellTypeThreshold` | `NULL` | Filter rare cell types |

### `VizSpatial()`

| Parameter | Notes |
|-----------|-------|
| `colBy` | colData column (discrete → categorical colors, numeric → gradient) |
| `reducedDim` | Column name within the reducedDim slot (e.g., a specific cell type name) |
| `reducedDim2use` | Which reducedDim slot (default `"PhiSpace"`) |
| `feature` | Gene name from `rownames(obj)` |
| `censor = TRUE` | Recommended for PhiSpace scores to improve contrast |

### `clusterPhiSpace()`

| Parameter | Default | Notes |
|-----------|---------|-------|
| `k` | `NULL` | Fixed k; mutually exclusive with `k_range` |
| `k_range` | `NULL` | Auto-select k via silhouette or elbow |
| `select_k_method` | `"silhouette"` | `"silhouette"` or `"elbow"` |
| `ncomp` | `NULL` → min(30, n-1) | PCs for clustering |
| `store_in_colData` | `FALSE` | Store result in `colData(x)[[cluster_name]]` |

---

## 10. Important Notes & Caveats

1. **Normalization must match**: Reference and query must use the **same assay normalization**. No batch correction is needed beyond this — rank transformation inherently handles different platforms.

2. **`ncomp` default**: Automatically set to the total number of phenotype categories. If `phenotypes = c("cell_type", "disease")` with 20 cell types and 3 disease states, `ncomp = 23`. Tune with `tunePhiSpace()` for better accuracy.

3. **`center=TRUE, scale=FALSE`**: PhiSpace uses all genes (not just 2,000 HVGs), so many low-variance genes are included. Scaling inflates these genes' contribution; centering only is the recommended default.

4. **Multi-reference column naming**: When using a named list of references, output columns are `"CellType(RefName)"`. Use this naming when calling `VizSpatial(reducedDim = "AT1(epithelial)")`.

5. **`refAssay` defaults differ**: `PhiSpace()` defaults to `refAssay = "rank"`, while `PhiSpaceR_1ref()` defaults to `"log1p"`. Always set explicitly to avoid surprises.

6. **Score interpretation**:
   - Column-wise normalized scores: high positive = confidently that cell type; near zero = average; negative = confidently not that cell type
   - Row-wise normalized scores (SingleR-style): 0–1 scale, 1 = highest scoring cell type for that cell

7. **`cellTypeThreshold`**: Applied before model fitting. Removes cells of rare types from the reference, preventing them from creating noisy prediction columns. Useful when reference has long-tail cell type distributions.

8. **SPE vs SCE for VizSpatial**: `SpatialExperiment` objects auto-detect coordinates from `spatialCoords()`. `SingleCellExperiment` objects require `x_coord` and `y_coord` arguments pointing to `colData` columns.

9. **`findNiches` vs `clusterPhiSpace`**: `findNiches()` is simpler and stores results directly in the object's `colData`. `clusterPhiSpace()` returns a `PhiSpaceClustering` object with diagnostics (elbow/silhouette plots, variance explained) and can auto-select k.

10. **Feature selection order**: If both `nfeat` and `selectedFeat` are given, `selectedFeat` takes priority. Feature selection is a two-step process internally: first fit PLS on all features to get importance scores, then refit on the selected union.

---

## 11. Object Slots Summary

After running `PhiSpace()`:

```
reducedDim(query, "PhiSpace")          # normalized scores: cells × phenotypes
reducedDim(query, "PhiSpace_nonNorm")  # raw scores (if storeUnNorm = TRUE)
colData(query)$PhiCellType             # discrete labels (after getClass())
colData(query)$spatial_niches          # niche assignments (after findNiches())
colData(query)$PhiClust                # cluster assignments (after clusterPhiSpace with store_in_colData)
```

After `spatialSmoother()`:

```
reducedDim(query, "PhiSpace_smoothed") # smoothed PhiSpace scores
```

After `findNiches()` with multiple k:

```
colData(spe)$spatial_niches_k7
colData(spe)$spatial_niches_k8
colData(spe)$spatial_niches_k9
S4Vectors::metadata(spe)$nicheAnalysis  # clustering metadata
```

---

## 12. Exported Functions Reference (38 exports)

| Function | Category | Description |
|----------|----------|-------------|
| `PhiSpace()` | Core | Main wrapper: single/multi reference annotation |
| `PhiSpaceR_1ref()` | Core | Single-reference core; returns full result list |
| `tunePhiSpace()` | Core | Cross-validate ncomp/nfeat |
| `normPhiScores()` | Core | Column-wise or row-wise score normalization |
| `mvr()` | Core | Internal PLS/PCA regression (exported) |
| `codeY()` | Encoding | SCE colData → dummy matrix |
| `codeY_vec()` | Encoding | Factor/vector → dummy matrix |
| `getClass()` | Scores | Discretize continuous scores → labels |
| `Score2Prob()` | Scores | Scores → probabilities |
| `VizSpatial()` | Spatial viz | Single-panel spatial plot |
| `saveCellTypeMaps()` | Spatial viz | Batch export all cell type heatmaps |
| `spatialSmoother()` | Spatial | KNN spatial smoothing |
| `findNiches()` | Spatial | Niche identification (simple interface) |
| `clusterPhiSpace()` | Spatial | Niche identification with diagnostics |
| `print.PhiSpaceClustering()` | S3 | Print method |
| `summary.PhiSpaceClustering()` | S3 | Summary method |
| `plot.PhiSpaceClustering()` | S3 | Plot method (pca/elbow/silhouette/variance) |
| `rankFeatures()` | Analysis | Feature importance ranking (PLSDA/PLS/DWD) |
| `cellTypeCorMat()` | Analysis | Cell type correlation matrix (S4 generic) |
| `plotPhiSpaceHeatMap()` | Viz | ComplexHeatmap of PhiSpace scores |
| `getPC()` | Utilities | PCA via partial SVD (irlba) |
| `computeUMAP()` | Utilities | UMAP computation |
| `selectFeat()` | Utilities | Select top features from importance scores |
| `RankTransf()` | Normalization | Rank transformation → "rank" assay |
| `logTransf()` | Normalization | Log/log1p transformation |
| `scranTransf()` | Normalization | scran normalization → "logcounts" |
| `CLRnorm()` | Normalization | Centered log-ratio normalization |
| `zeroFeatQC()` | QC | Remove all-zero genes |
| `KeepCommonGenes()` | QC | Intersect genes between datasets |
| `subsample()` | Sampling | Stratified subsampling of SCE |
| `pseudoBulk()` | Aggregation | Aggregate to pseudobulk |
| `spatialSampler()` | Sampling | Spatial sampling utilities |
| `translateLabel()` | Utilities | Translate between label schemes |
| `reNameCols()` | Utilities | Rename matrix columns with prefix |
| `doubleCent()` | Utilities | Double-center matrix by row/column means |
| `Score2Prob()` | Scores | Convert PhiSpace scores to probabilities |

---

## 13. Build & Development

```bash
# Regenerate documentation (man/ + NAMESPACE)
Rscript -e 'devtools::document("pkg")'

# Run R CMD check
Rscript -e 'devtools::check("pkg")'

# Run unit tests
Rscript -e 'devtools::test("pkg")'

# Install locally
Rscript -e 'devtools::install("pkg")'

# Build pkgdown site
Rscript -e 'pkgdown::build_site("pkg")'
```

**Package source layout:**
```
pkg/
├── R/           # 41 R source files (edit these)
├── man/         # auto-generated docs (never edit directly)
├── vignettes/   # 7 Rmd vignettes
├── tests/       # testthat unit tests
├── DESCRIPTION  # package metadata
└── NAMESPACE    # auto-generated exports/imports
```
