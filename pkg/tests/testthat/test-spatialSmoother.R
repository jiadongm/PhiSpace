test_that("spatialSmoother adds a smoothed reduced dimension", {
  spe <- make_spatial_fixture()
  out <- spatialSmoother(spe, k = 5, verbose = FALSE)

  smoothed <- SingleCellExperiment::reducedDim(out, "PhiSpace_smoothed")
  original <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  expect_equal(dim(smoothed), dim(original))
  expect_equal(dimnames(smoothed), dimnames(original))
  # Averaging over spatial neighbours reduces within-group noise.
  expect_lt(sd(smoothed[spe$group == "g1", 4]), sd(original[spe$group == "g1", 4]))
})

test_that("spatialSmoother keeps a constant assay constant", {
  spe <- make_spatial_fixture()
  SummarizedExperiment::assay(spe, "logcounts")[] <- 3
  out <- spatialSmoother(spe, smoothReducedDim = FALSE, k = 5, verbose = FALSE)

  expect_equal(unname(SummarizedExperiment::assay(out, "logcounts_smoothed")),
               matrix(3, nrow(spe), ncol(spe)))
})

test_that("spatialSmoother reports with messages and verbose = FALSE is silent", {
  spe <- make_spatial_fixture()
  loud <- capture_all(spatialSmoother(spe, k = 5))
  quiet <- capture_all(spatialSmoother(spe, k = 5, verbose = FALSE))

  expect_true("Added smoothed reduced dimension: PhiSpace_smoothed" %in% sub("\n$", "", loud$messages))
  expect_length(quiet$output, 0)
  expect_length(quiet$messages, 0)
})

test_that("spatialSmoother uses colData coordinates for a SingleCellExperiment", {
  spe <- make_spatial_fixture()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = SummarizedExperiment::assay(spe, "logcounts")),
    colData = cbind(SummarizedExperiment::colData(spe), SpatialExperiment::spatialCoords(spe)),
    reducedDims = list(PhiSpace = SingleCellExperiment::reducedDim(spe, "PhiSpace"))
  )
  from_sce <- spatialSmoother(sce, x_coord = "x", y_coord = "y", k = 5, verbose = FALSE)
  from_spe <- spatialSmoother(spe, k = 5, verbose = FALSE)

  expect_equal(SingleCellExperiment::reducedDim(from_sce, "PhiSpace_smoothed"),
               SingleCellExperiment::reducedDim(from_spe, "PhiSpace_smoothed"))
  expect_error(spatialSmoother(sce, k = 5), "x_coord and y_coord must be specified")
})

test_that("spatial_smooth_expression matches a weighted average over neighbours", {
  spe <- make_spatial_fixture()
  coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
  X <- Matrix::Matrix(SummarizedExperiment::assay(spe, "counts"), sparse = TRUE)
  for (kernel in c("linear", "gaussian", "uniform")) {
    for (self in c(TRUE, FALSE)) {
      knn <- FNN::get.knnx(coords, coords, k = if (self) 6 else 7)
      idx <- knn$nn.index
      dst <- knn$nn.dist
      if (!self) {
        idx <- idx[, -1]
        dst <- dst[, -1]
      }
      w <- compute_kernel_weights(dst, kernel)
      expected <- sapply(seq_len(ncol(X)), function(i) {
        as.vector(as.matrix(X[, idx[i, ]]) %*% (w[i, ] / sum(w[i, ])))
      })
      out <- spatial_smooth_expression(X, coords, k = 6, kernel = kernel,
                                       include_self = self, verbose = FALSE)
      expect_true(is.matrix(out))
      expect_equal(dimnames(out), dimnames(X))
      expect_equal(out, expected, ignore_attr = TRUE, tolerance = 1e-12)
    }
  }
})
