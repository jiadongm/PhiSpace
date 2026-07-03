# Handoff: revise `PhiSpace::scoreCells()`

## Target repository

Open a new Codex project/session at:

```text
/Users/jmao1/Documents/GitHub/PhiSpace/Pkg
```

This handoff was written from:

```text
/Users/jmao1/Library/CloudStorage/Dropbox/Research_projects/CaseMacroFragment
```

The installed PhiSpace version inspected during planning was PhiSpace 1.1.0.
Do not patch the installed package under the system R library; make all changes
in the source repository above.

## Requested changes

Revise `PhiSpace::scoreCells()` so that:

1. Correlation-based scoring uses only the signature genes for the relevant
   reference class, rather than all genes shared between the reference and
   query.
2. Users can provide their own signature genes instead of requiring
   `scoreCells()` to derive signatures from the reference.
3. When `scoring = "both"`, correlation and signature-expression scoring use
   exactly the same effective signature genes.

## Current behavior

The current control flow is approximately:

```text
scoreCells()
    -> define classes with .get_score_cells_classes()
    -> restrict reference and query to shared genes
    -> derive centroids and signatures with .buildClassSignatures()
    -> correlation:
           .scoreCorrelation(query_expr, centroids, method)
           uses all shared genes
    -> signature:
           .scoreSignature(query_expr, signatures, zscore)
           uses signature genes
```

Important current details:

- `class_col = NULL` assigns every reference column to a synthetic class named
  `"reference"`.
- A one-class reference always uses mean-expression signature selection,
  because `scran::findMarkers()` is unavailable with only one class.
- Mean selection ranks genes by class centroid expression, truncates to
  `n_top_genes`, applies the housekeeping exclusion, and retains up to
  `n_signature_genes`.
- Default housekeeping patterns are:
  `^Rpl`, `^Rps`, `^mt-`, `^MT-`, `^Mrpl`, and `^Mrps`.
- Signature scoring calculates the mean query expression over the signature.
  With one class, it z-scores across query cells. With multiple classes, it
  z-scores across classes within each query cell.
- Correlation scoring currently calls `stats::cor()` between each query column
  and each class centroid across all shared genes.

## Proposed public API

Add an argument:

```r
signature_genes = NULL
```

Supported forms:

```r
# Automatically derive signatures: preserve existing selection behavior.
signature_genes = NULL

# One-class reference.
signature_genes = c("C1qa", "C1qb", "C1qc", "Csf1r")

# Multi-class reference.
signature_genes = list(
  macrophage = c("C1qa", "C1qb", "Csf1r"),
  monocyte = c("Lyz2", "Ccr2", "Ly6c2")
)
```

Recommended validation rules:

- A character vector is accepted only when the reference has one class. Convert
  it internally to a named list using the sole class name.
- A multi-class reference requires a named list.
- List names must correspond exactly to the reference class names. Fail on
  missing, unknown, or duplicated class names rather than silently guessing.
- Every signature must be a non-empty character vector without `NA` or blank
  values.
- Remove duplicated genes within a signature while preserving order.
- Intersect each signature with the genes available in both the selected
  reference assay and selected query assay.
- Report supplied genes that are unavailable. Prefer a warning containing the
  affected class and gene count; retain the effective intersection.
- Correlation requires at least two effective genes and nonzero centroid
  variance. Fail clearly when a class has too few genes. Allow ordinary
  per-cell zero-variance cases to produce `NA` with an informative warning or
  documented behavior.
- A user-provided signature is authoritative: do not remove housekeeping genes
  and do not re-rank or truncate it using `n_top_genes` or
  `n_signature_genes`.
- Document that marker-selection arguments are ignored when
  `signature_genes` is supplied. Avoid noisy warnings unless that matches the
  package's established style.

Before coding, check the package's existing API conventions and tests. Adjust
the argument name only if the repository already has a clearer established
name, but retain the behavior above.

## Internal design

Normalize both generated and supplied signatures into the same representation:

```r
signatures <- list(
  class_1 = c("geneA", "geneB"),
  class_2 = c("geneC", "geneD")
)
```

Suggested flow:

```text
classes <- get_reference_classes(reference, class_col)
shared_genes <- intersect(rownames(reference), rownames(query))
reference <- reference[shared_genes, ]
query <- query[shared_genes, ]

IF signature_genes is NULL:
    build <- .buildClassSignatures(...)
    signatures <- build$signatures
    signature_source <- "generated"
ELSE:
    centroids <- .compute_class_centroids(ref_expr, classes)
    signatures <- .prepare_user_signatures(
        signature_genes,
        classes,
        shared_genes
    )
    marker_stats <- suitable empty/NA result consistent with package API
    signature_source <- "user"

IF scoring includes "correlation":
    correlation_scores <- .scoreCorrelation(
        query_expr,
        centroids,
        signatures,
        method
    )

IF scoring includes "signature":
    signature_scores <- .scoreSignature(
        query_expr,
        signatures,
        zscore
    )
```

Refactor `.scoreCorrelation()` to score each class independently:

```text
FOR each class:
    genes <- intersect(signatures[[class]], rownames(query_expr))

    correlation_scores[, class] <-
        correlation(
            each query cell's expression across genes,
            centroids[genes, class],
            method = cor_method
        )
```

For multiple classes, different class scores may use different gene sets. This
is intentional. Document that cross-class correlation scores can therefore be
based on class-specific features.

Do not inadvertently transpose the output. Preserve:

```text
rows    = query cells/samples
columns = reference classes
```

and preserve query column names and class names.

## Metadata

Continue storing the effective signatures in:

```r
S4Vectors::metadata(query)$scoreCells$signatures
```

Also record enough provenance to reproduce the run. Suggested additions:

```r
metadata(query)$scoreCells$signature_source
# "generated" or "user"

metadata(query)$scoreCells$correlation_genes
# named list of effective genes per class

metadata(query)$scoreCells$params$signature_genes_supplied
# logical, or another compact nonduplicative provenance field
```

Do not store a misleading generated-marker table for supplied signatures.
Either use an empty data frame with stable columns or document `NULL`, following
the package's existing metadata conventions.

## Required tests

First run the existing test suite unchanged to establish the baseline. Then add
focused synthetic tests covering the following.

### Generated signatures

- Single-class generated signature correlation equals a manual
  `stats::cor()` calculation restricted to the generated signature genes.
- Multi-class generated correlation uses each class's own generated signature.
- Adding or dramatically changing a non-signature gene does not change the
  correlation score.
- Existing signature-expression scores are numerically unchanged.

### User signatures

- A character vector works for a one-class reference.
- A named list works for a multi-class reference.
- Each multi-class correlation agrees with a manual correlation using that
  class's supplied genes.
- `scoring = "signature"` uses the supplied genes.
- `scoring = "both"` uses identical effective genes for both scoring branches.
- User-provided housekeeping genes remain in the signature.
- User signatures bypass mean ranking and `scran::findMarkers()`.

### Validation and edge cases

- Reject an unnamed list for a multi-class reference.
- Reject missing, extra, or duplicated class names.
- Reject non-character entries, `NA`, blank strings, and empty signatures.
- Deduplicate repeated genes deterministically.
- Warn about genes missing from the reference/query.
- Fail when fewer than two effective genes remain for correlation.
- Test constant centroid values and constant query profiles.
- Preserve score dimensions and dimnames.
- Preserve behavior for `SingleCellExperiment` and, if already tested by the
  package, `SpatialExperiment`.
- Preserve custom `score_prefix`.
- Ensure `scoring = "correlation"`, `"signature"`, and `"both"` create only the
  expected reduced dimensions.

Use small deterministic fixtures. Tests should compare numerical results
directly rather than merely checking dimensions.

## Documentation

Update the roxygen documentation and regenerate `man/` files using the
repository's established workflow.

Clearly explain:

```text
correlation score:
    similarity of the relative expression pattern across signature genes

signature score:
    average expression level across signature genes, optionally z-scored
```

Add examples for:

- an automatically generated one-class signature;
- a supplied one-class character vector;
- supplied multi-class named signatures;
- `scoring = "both"`.

Update `NEWS.md` and the package version if required by repository policy.
Call out the correlation change as user-visible: prior versions used all shared
genes, whereas the revised function uses signature genes.

## Acceptance criteria

The work is complete when:

1. Correlation scores are calculated only from each class's effective signature
   genes.
2. Users can supply signatures for one-class and multi-class references.
3. Generated-signature behavior remains available with
   `signature_genes = NULL`.
4. Signature-expression results remain backward-compatible when no custom
   signature is supplied.
5. Metadata identifies the signature source and effective genes.
6. New unit tests pass.
7. The full package test suite and package check pass without new warnings,
   notes, or errors attributable to the change.

## Downstream CaseMacroFragment sensitivity analysis

After the package change is complete and installed from source, return to:

```text
/Users/jmao1/Library/CloudStorage/Dropbox/Research_projects/CaseMacroFragment
```

The current reporter-BM-to-TMS analysis is:

```text
atlas/score_tms_cells_pure_bm_phispace_null.R
```

Its locked 50-gene signature is:

```text
atlas/figs/TMS_signature_propagation_PhiSpace/
pure_bm_signature_locked.csv
```

Run a downstream sensitivity comparison using:

1. The existing single-class mean-expression signature score.
2. Correlation using the automatically generated 50-gene signature.
3. Correlation using the supplied locked 50-gene signature.

Do not overwrite the completed corrected analysis. Write sensitivity outputs to
a new directory and record the installed development version/commit of
PhiSpace.

The biological interpretation should keep the score types distinct:

- Signature expression asks whether the reporter programme is elevated.
- Correlation asks whether the relative pattern among those genes resembles
  the reporter centroid.

