# Extra PhiSpace Functions

Plans for new PhiSpace utilities motivated by the CaseCRCLM analyses. Currently scoped to two functions:

1. `scoreCells()` — correlation- and signature-based per-cell-per-class scoring (PhiSpace fallback for few/single-class references).
2. `piScore()` — combine log fold change and p-values into a single significance score for ranking DE genes.

---

## `scoreCells()`

### Motivation

PhiSpace's discriminative phenotype-space model needs a reference with multiple classes to be meaningful. Two regimes break this assumption:

1. **Single-class reference** — e.g. Juan's CD/HFD cancer-cell-only references. PhiSpace cannot run at all.
2. **Very few classes** (2–3), or very few cells per class — PhiSpace can technically run but is unstable.

In both cases we still want to ask the same biological question PhiSpace answers: *for each query cell/bin, how strongly does it resemble each reference class?* `scoreCells()` provides a lightweight, correlation- and signature-based fallback that produces the same shape of output (a `[n_query × n_classes]` score matrix stored in a `reducedDim`) so downstream PhiSpace tooling (visualisation, smoothing, clustering) works unchanged.

The name is `scoreCells` (not `scoreSpatialBins`) because the function applies equally to spatial bins and to ordinary single cells — spatial smoothing is a separate, optional follow-up step.

### Origin: what we did for Juan

In `score_cancerArea_CD.R` and the generalised `score_cancerArea_all.R`, we scored Stereo-seq spatial bins against single-class (cancer-cell-only) references using two complementary methods:

- **Correlation score**: Spearman correlation of each bin's expression vector with the reference centroid (mean expression vector across all cancer cells).
- **Signature score**: mean expression of a curated 100-gene cancer signature, z-scored.

The signature was built by ranking genes by mean expression in the reference, taking the top 300, stripping housekeeping genes (regex patterns plus a curated list), and keeping the top 100. Both scores were spatially smoothed with k-NN averaging.

Empirically the two scoring methods agreed strongly (r ≈ 0.94) and were stable across reference conditions (CD vs HFD vs merged). `scoreCells()` is the multi-class generalisation of this pipeline, packaged as a reusable PhiSpace function.

### Behaviour

`scoreCells()` produces, for each query cell/bin and each reference class, two scores:

- A **correlation score** (Spearman or Pearson) between the cell's expression and the class's centroid.
- A **signature score**: the mean expression of a class-specific gene signature, z-scored across classes within each cell so rows are directly comparable.

Both are written into the query object as `reducedDim` slots, mirroring how PhiSpace stores its own annotation scores.

### Design decisions (settled)

| Question | Decision |
|---|---|
| Multi-class signature selection | `scran::findMarkers` (test type configurable; default `t`, `pval.type="any"`, `direction="up"`). |
| Single-class signature selection | Per-class mean-rank of all genes, minus housekeeping (the existing CRCLM approach). `findMarkers` falls back to this automatically when only one class is present. |
| Shared genes across class signatures | **Allowed** — a gene can be in multiple classes' signatures. (`pval.type="any"` rather than `"all"`.) |
| Z-scoring of signature scores | **Across classes per cell** (within each row), so values within a cell are directly comparable regardless of overall expression magnitude. |
| Spatial smoothing | **Not done inside `scoreCells()`.** Reuse `PhiSpace::spatialSmoother()` as a separate, optional step on the resulting `reducedDim` slot. |
| Output container | The query object itself (SCE/SPE), with new `reducedDim` slots and auxiliary outputs in `metadata()`. Mirrors `spatialSmoother()`'s in-place style. |

### Function signature (sketch)

```r
scoreCells(
  reference, query,
  class_col = NULL,                 # colData(reference) col; NULL/1-class -> single-class mode
  assay_ref   = "logcounts",
  assay_query = "logcounts",

  ## Signature selection
  signature_method   = c("mean", "findMarkers"),
  marker_test        = c("t", "wilcox", "binom"),
  marker_pval_type   = "any",
  marker_direction   = "up",
  n_top_genes        = 300,         # candidate pool per class
  n_signature_genes  = 100,         # final signature size per class
  remove_housekeeping   = TRUE,
  housekeeping_patterns = c("^Rpl","^Rps","^mt-","^Mrpl","^Mrps"),
  housekeeping_genes    = NULL,     # extra curated list; built-in default
  shared_genes_only     = TRUE,     # intersect ref & query feature sets

  ## Scoring
  scoring          = c("both", "correlation", "signature"),
  cor_method       = c("spearman", "pearson"),
  zscore_signature = TRUE,          # across classes per cell

  ## Output / misc
  score_prefix = "scoreCells",
  BPPARAM      = SerialParam(),
  verbose      = TRUE
)
```

Returns the (modified) query object with:

- `reducedDim(query, "<prefix>_correlation")` — `[n_query × n_classes]` matrix.
- `reducedDim(query, "<prefix>_signature")`   — `[n_query × n_classes]` matrix, z-scored across classes per row.
- `metadata(query)$scoreCells` — list with: `signatures` (named list of gene vectors per class, possibly overlapping), `marker_stats` (data.frame: gene, class, summary.logFC, FDR, rank, is_housekeeping, kept), `centroids` (`[n_genes × n_classes]`), `classes`, `params`.

### Internal helpers (exported for composability)

- `.buildClassSignatures()` — builds per-class signatures via `findMarkers` (multi-class) or per-class mean rank (single-class), strips housekeeping, returns signatures + `marker_stats` + centroids.
- `.scoreCorrelation()` — vectorised rank-based Spearman (or Pearson) of each cell against each class centroid; returns `[n_query × n_classes]`.
- `.scoreSignature()` — per-class mean expression of signature genes; z-score across classes per cell.

### Typical usage

```r
# Spatial bins, multi-class reference
query <- scoreCells(reference, query, class_col = "celltype")
query <- spatialSmoother(query,
                         reducedDim2smooth = "scoreCells_signature",
                         k = 10, kernel = "linear")
# Adds reducedDim "scoreCells_signature_smoothed"

# Single-class cancer reference (the original Juan use case)
query <- scoreCells(reference, query, class_col = NULL)
query <- spatialSmoother(query, reducedDim2smooth = "scoreCells_correlation", k = 10)

# Plain single cells (non-spatial) — skip smoothing
query <- scoreCells(reference, query, class_col = "celltype")
```

### PhiSpace integration

- Public function exported alongside the existing PhiSpace annotation entry points.
- PhiSpace's main annotation function can dispatch automatically: if `n_classes < threshold` (or per-class sample size is too small), call `scoreCells()` instead and wrap its `reducedDim` output in the same score container PhiSpace returns. Downstream visualisation / smoothing / clustering code is then unchanged.
- The `score_prefix` argument lets callers rename the slot to match PhiSpace's existing reducedDim conventions when used as a fallback.

### Open / deferred items

- **Performance**: `spatialSmoother()` rebuilds the kNN graph on every call, so smoothing both `_correlation` and `_signature` pays that cost twice. Minor for now; if it becomes a bottleneck, add an optional precomputed-kNN argument to `spatialSmoother()` itself rather than to `scoreCells()`.
- **Marker test default**: `t` is fast and a reasonable default but `wilcox` may be more robust for sparse data. Revisit after testing on Juan's CD/HFD references.
- **Built-in housekeeping list**: ship the curated list from `score_cancerArea_CD.R` as the default for `housekeeping_genes`.
- **Validation plan**: re-run `scoreCells()` on the CD/HFD/All references in single-class mode and confirm the resulting scores match (within tolerance) the existing `score_cancerArea_*.R` outputs before considering the function ready for PhiSpace.

---

## `piScore()`

### Motivation

After a DE test (Seurat `FindMarkers`, `scran::findMarkers`, limma, edgeR, DESeq2, MAST, …) we routinely need a single number per gene to rank by — neither raw p-value nor fold change alone is right. The π-score of Xiao et al. (2014) combines them into one statistic and is what we have used informally across PhiSpace case studies (e.g. `CaseActivationAtlas/utils.R`, `CaseFibrosis/utils.R`, `CaseVisium/utils.R` all compute `-log10(p_val_adj) * avg_log2FC` inline). Promoting it to a small, reusable PhiSpace function removes the copy-paste and centralises the choice of conventions.

### Definition (Xiao et al. 2014, *Bioinformatics* 30(6):801–807, [doi:10.1093/bioinformatics/btr671](https://doi.org/10.1093/bioinformatics/btr671))

The paper's definition is:

```
π = |log2(FC)| × (-log10(p))
```

i.e. **unsigned**. The paper notes:

- p can be raw or adjusted; the score is not tied to any particular DE test.
- Suggested ad-hoc thresholds: π > 1.3 (≈ 2-fold change at p = 0.05) or π > 2.0 (≈ 2-fold change at p = 0.01).

### Convention adopted in PhiSpace usage

Across the existing case-study utilities we have actually used the **signed** version:

```
π_signed = log2(FC) × (-log10(p))
```

This preserves direction (positive = up-regulated, negative = down-regulated) so that `rank(-π_signed)` returns the top up-regulated genes directly. That is the more useful form when picking class-specific markers.

`piScore()` should therefore expose both via a `signed` argument (default `TRUE` to match how we already use it; `FALSE` for the strict Xiao et al. definition).

### Function signature (sketch)

```r
piScore(
  de_results,                       # data.frame from a DE test
  logfc_col   = "avg_log2FC",       # Seurat default; auto-detect common alternatives
  pval_col    = "p_val_adj",        # adjusted by default; user may pass raw
  signed      = TRUE,               # TRUE = signed (PhiSpace convention),
                                    # FALSE = |log2FC| (Xiao et al. 2014)
  pval_floor  = .Machine$double.xmin,  # clamp p before -log10 to avoid Inf
  add_column  = TRUE,               # append "pi_score" column to de_results
  rank        = FALSE               # if TRUE also append "pi_rank" (descending)
)
```

Returns the input `data.frame` with `pi_score` (and optionally `pi_rank`) added, or — if `add_column = FALSE` — a numeric vector aligned with `rownames(de_results)`.

### Auto-detection of column names

Common DE tools use different conventions; the function should fall back gracefully if defaults aren't found:

| Tool | logFC column | p-value column |
|---|---|---|
| Seurat `FindMarkers` | `avg_log2FC` | `p_val_adj` (raw: `p_val`) |
| `scran::findMarkers` | `summary.logFC` | `FDR` |
| limma `topTable` | `logFC` | `adj.P.Val` |
| edgeR `topTags` | `logFC` | `FDR` |
| DESeq2 `results` | `log2FoldChange` | `padj` |

A small `.detect_de_columns()` helper picks the first matching pair present in the input, with explicit user-supplied `logfc_col` / `pval_col` overriding detection.

### Edge cases

- **`p = 0`** (common with very small p-values that underflow): clamp to `pval_floor` to avoid `Inf`. Warn if any clamping happens.
- **`p = NA`** or **`logFC = NA`**: pass through as `NA` in the output.
- **Adjusted vs raw p-values**: leave the choice to the caller via `pval_col`. Document the trade-off (adjusted is the conservative default; raw preserves more dynamic range when ranking is the only goal).

### Typical usage

```r
markers <- Seurat::FindMarkers(seurat_obj, ident.1 = "tumor")
markers <- piScore(markers, rank = TRUE)
top_up   <- head(markers[order(-markers$pi_score), ], 50)
top_down <- head(markers[order( markers$pi_score), ], 50)
```

Companion use inside `scoreCells()`: when `signature_method = "findMarkers"`, after the `findMarkers` call we can rank candidate genes per class by `piScore()` rather than by raw p-value, before taking the top `n_top_genes`. This unifies the two functions' gene-selection logic.

### Citation requirement

The final function **must cite Xiao et al. (2014) within the code itself** — both in the roxygen2 `@references` block and as an inline comment next to the formula. Required citation:

> Xiao Y, Hsiao T-H, Suresh U, et al. A novel significance score for gene selection and ranking. *Bioinformatics* 30(6):801–807 (2014). doi:10.1093/bioinformatics/btr671

Roxygen skeleton to use verbatim:

```r
#' Compute pi-score for ranking differentially expressed genes
#'
#' Combines log fold change and p-value into a single significance score
#' (Xiao et al. 2014). By default returns the signed variant
#' \code{log2FC * -log10(p)} for directional ranking; set \code{signed = FALSE}
#' for the paper's original \code{|log2FC| * -log10(p)}.
#'
#' @references
#' Xiao Y, Hsiao T-H, Suresh U, et al. (2014) A novel significance score for
#'   gene selection and ranking. \emph{Bioinformatics} 30(6):801-807.
#'   \doi{10.1093/bioinformatics/btr671}
```

### Open / deferred items

- **Default `signed`**: TRUE matches existing PhiSpace usage but diverges from the paper's definition. Worth flagging in the function docs.
- **Threshold helper**: optional convenience function `piThreshold(score, cutoff = 1.3)` returning a logical mask, mirroring the paper's suggested cutoffs.
