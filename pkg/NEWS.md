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

## Bug fixes

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
