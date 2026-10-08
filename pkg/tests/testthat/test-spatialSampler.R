test_that("spatialSampler random method returns the target number of cells", {
  spe <- make_spatial_fixture()
  run <- capture_all(spatialSampler(spe, prop = 0.2, method = "random"))

  expect_s4_class(run$value, "SpatialExperiment")
  expect_equal(ncol(run$value), 60)
  expect_true(all(colnames(run$value) %in% colnames(spe)))
  meta <- S4Vectors::metadata(run$value)$spatialSampling[[1]]
  expect_equal(meta$method, "random")
  expect_equal(meta$sampled_n_cells, 60)
})

test_that("spatialSampler grid and kmeans methods sample every region", {
  spe <- make_spatial_fixture()
  for (method in c("grid", "kmeans")) {
    out <- spatialSampler(spe, prop = 0.2, method = method, verbose = FALSE)
    expect_gt(ncol(out), 0)
    expect_setequal(unique(out$group), c("g1", "g2", "g3"))
  }
})

test_that("spatialSampler is reproducible with the same seed", {
  spe <- make_spatial_fixture()
  a <- spatialSampler(spe, prop = 0.2, seed = 7, verbose = FALSE)
  b <- spatialSampler(spe, prop = 0.2, seed = 7, verbose = FALSE)
  expect_identical(colnames(a), colnames(b))
})

test_that("spatialSampler reports with messages and verbose = FALSE is silent", {
  spe <- make_spatial_fixture()
  loud <- capture_all(spatialSampler(spe, prop = 0.2))
  quiet <- capture_all(spatialSampler(spe, prop = 0.2, verbose = FALSE))

  expect_length(loud$output, 0)
  expect_true(any(grepl("^Spatial sampling: 60 cells from 300 total cells", loud$messages)))
  expect_length(quiet$output, 0)
  expect_length(quiet$messages, 0)
  expect_identical(colnames(quiet$value), colnames(loud$value))
})

test_that("spatialSampler validates its inputs", {
  spe <- make_spatial_fixture()
  expect_error(spatialSampler(spe, prop = 0), "prop must be between 0 and 1")
  expect_error(spatialSampler(spe, method = "hex"), "method must be")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1, 2, 2)))
  expect_error(spatialSampler(sce), "must be a SpatialExperiment")
})
