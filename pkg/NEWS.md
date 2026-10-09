# PhiSpace 1.1.0.9000

## Changes that can alter results

* `getPC()` now computes `totVar`, `props` and `accuProps` from the centred
  (and, if requested, scaled) data that it decomposes. Previous versions divided
  by the total variance of the raw data, so the proportions were too small
  whenever column means were not zero. Scores and loadings are unchanged. The
  corrected proportions also change the `variance_explained` value that
  `findNiches()` stores in its metadata.
* `clusterPhiSpace()` now computes its principal components with a full
  singular value decomposition instead of irlba. irlba could return inaccurate
  trailing components, and the default `ncomp` (up to 30 components) could
  include them in the clustering. Clusters can therefore differ from previous
  versions, mainly for score matrices with few columns.
* PLS fits (`mvr(method = "PLS")`, and therefore `PhiSpace()`,
  `PhiSpaceR_1ref()` and `tunePhiSpace()`) now find each weight vector with an
  exact eigendecomposition (`eigen()`) instead of `irlba::partial_eigen()`.
  irlba stopped at a tolerance of about 1e-5 and used a random start, so
  coefficients differed from `pls::kernelpls.fit` by up to about 2e-5
  (relative) and varied with the random seed. They now agree with
  `pls::kernelpls.fit` to about 1e-14 and do not depend on the seed. Scores
  change by a similarly small amount.
* With a `"rank"` assay, `PhiSpace()`, `PhiSpaceR_1ref()` and
  `tunePhiSpace()` rank the selected features again after feature selection.
  This second ranking now ranks the genes within each cell, as
  `RankTransf()` does. Previous versions ranked each gene across cells, so a
  cell's ranked values depended on the other cells in the same reference or
  query. Scores computed from a `"rank"` assay change.
* `mvr(method = "PCA")` now uses a full singular value decomposition when
  `ncomp` is at least half of the smaller dimension of `X`, where irlba can be
  inaccurate.

## Memory and speed

* `mvr()` and `phenotype()` no longer convert a sparse predictor matrix to a
  dense one. Centring and scaling are applied implicitly, inside the matrix
  products. On a simulated 20,000 x 8,000 sparse matrix (10% non-zero, 20
  responses), peak R memory (including the 185 MB input) fell from 3,574 MB
  to 906 MB for the PLS fit and from 4,646 MB to 988 MB for prediction; run
  times fell from 9.3 s to 0.9 s and from 4.0 s to 0.4 s. In fresh vignette
  builds, peak memory fell from 17-24 GB to 7.2 GB for BridgeAnnotation and
  from 64.9 GB to 16.1 GB for CITE-seq.
* `mvr()` gains a `keepComps` argument that selects the numbers of components
  whose coefficients are returned. The default, `1:ncomp`, returns all of them
  as before. `PhiSpace()`, `PhiSpaceR_1ref()`, `tunePhiSpace()` and
  `rankFeatures()` keep only the final coefficients. As a result, the
  `atlas_re$reg_re$coefficients` array returned by `PhiSpaceR_1ref()` now has
  one slice, named `"<ncomp> comps"`; select it by that name, not by
  position. `phenotype()` accepts both forms, so objects saved by earlier
  versions still work.
* `RankTransf()` no longer converts a sparse count matrix with no negative
  values to a dense one. The ranks are unchanged.
* `PhiSpaceR_1ref()` (and therefore `PhiSpace()`) no longer subsets the
  reference and query objects to their shared genes, which copied every
  assay. It extracts only the assay it uses and keeps the selected genes
  before it transposes the matrix. Results are unchanged. On simulated data
  (3.3 GB of input: a 10,000-cell reference and a 40,000-cell query, each
  with three assays), peak R memory fell from 7.4 GB to 3.7 GB.
* `spatialSmoother()` smooths an assay with one sparse cell-by-cell weight
  matrix instead of a loop over cells, and `pseudoBulk()` aggregates with a
  sparse indicator matrix instead of a loop over pseudo-bulks.
  `scoreCells()` computes correlation scores for blocks of cells at once and
  no longer holds the class's whole query matrix as a dense matrix. Results
  agree with the previous versions to about 2e-15. On simulated data (5,000
  genes, 30,000 cells), smoothing an assay took 8 s instead of 332 s,
  `pseudoBulk()` 0.7 s instead of 22 s, and Spearman correlation scoring
  (25,000 query cells, 20 classes) 5 s instead of 26 s.

## Bug fixes

* `pseudoBulk()` drew the pseudo-bulks of a cluster with only one cell from
  the wrong cells: for the cell at position i, `sample()` drew from cells 1
  to i. It now draws from that cell only. Clusters with two or more cells
  give the same pseudo-bulks as before, for the same `seed`.
* `print()`, `summary()` and `plot()` for `clusterPhiSpace()` results now read
  the variance proportions that `getPC()` returns. Previously `print()`
  reported 0% variance explained, `summary()` and `plot(type = "variance")`
  failed, and the PCA plot axis labels showed no percentages.

## Progress messages

* `findNiches()`, `zeroFeatQC()`, `spatialSampler()` and `spatialSmoother()`
  now report progress with `message()` instead of `cat()`, so
  `suppressMessages()` and the knitr chunk option `message = FALSE` silence
  them. `findNiches()`, `zeroFeatQC()` and `spatialSampler()` gain a
  `verbose` argument (default `TRUE`); `spatialSmoother()` already had one.
* `findNiches()` reports niche sizes on one line instead of printing a table.

## Other changes

* The `License` field changed from `AGPL (>= 3)` to `AGPL-3`.
* Raw scores from `phenotype()` (for example `PhiSpaceScore` and `YrefHat`)
  no longer carry a `"scaled:center"` attribute.
* `pls` is now a suggested package; the unit tests compare `mvr()` with
  `pls::kernelpls.fit`.
* The `PhiSpaceR_1ref()` help page now documents every value that the
  `scoreCells` fallback changes.

## scoreCells()

* `scoreCells()` now calculates each class's correlation score using that
  class's effective signature genes. Previous versions used every gene shared
  by the reference and query.
* `scoreCells()` gains a `signature_genes` argument for supplying a character
  vector to a one-class reference or an exactly named list to a multi-class
  reference.
* Runs now record whether signatures were generated or supplied and which
  effective genes were used for correlation scoring.
