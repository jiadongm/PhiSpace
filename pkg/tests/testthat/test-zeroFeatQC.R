test_that("zeroFeatQC removes features that are zero in every cell", {
  spe <- make_spatial_fixture()
  SummarizedExperiment::assay(spe, "counts")[c(2, 5), ] <- 0
  run <- capture_all(zeroFeatQC(spe))

  expect_equal(nrow(run$value), nrow(spe) - 2)
  expect_false(any(c("gene2", "gene5") %in% rownames(run$value)))
  expect_length(run$output, 0)
  expect_identical(run$messages, "Deleted features with all zeros.\n")
})

test_that("zeroFeatQC keeps every feature when none is all zero", {
  spe <- make_spatial_fixture()
  SummarizedExperiment::assay(spe, "counts")[, 1] <- 1
  run <- capture_all(zeroFeatQC(spe))

  expect_equal(dim(run$value), dim(spe))
  expect_identical(run$messages, "All features have at least 1 nonzero value.\n")
})

test_that("zeroFeatQC verbose = FALSE is silent", {
  spe <- make_spatial_fixture()
  run <- capture_all(zeroFeatQC(spe, verbose = FALSE))

  expect_length(run$output, 0)
  expect_length(run$messages, 0)
})
