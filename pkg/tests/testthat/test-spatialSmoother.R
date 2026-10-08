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
