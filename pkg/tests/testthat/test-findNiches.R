test_that("findNiches recovers well-separated niches", {
  spe <- make_spatial_fixture()
  res <- capture_all(findNiches(spe, n_niches = 3, ncomp = 3))$value

  expect_s3_class(res$spatial_niches, "factor")
  expect_equal(levels(res$spatial_niches), paste0("Niche_", 1:3))
  expect_true(same_partition(res$spatial_niches, res$group))
})

test_that("findNiches reports progress with messages, not standard output", {
  spe <- make_spatial_fixture()
  run <- capture_all(findNiches(spe, n_niches = 3, ncomp = 3))

  expect_length(run$output, 0)
  expect_true(any(grepl("^Performing k-means clustering with k = 3", run$messages)))
  expect_true(any(grepl("^Niche sizes for k = 3: Niche_1 = ", run$messages)))
})

test_that("findNiches verbose = FALSE is silent and gives the same niches", {
  spe <- make_spatial_fixture()
  loud <- capture_all(findNiches(spe, n_niches = 3, ncomp = 3))$value
  quiet <- capture_all(findNiches(spe, n_niches = 3, ncomp = 3, verbose = FALSE))

  expect_length(quiet$output, 0)
  expect_length(quiet$messages, 0)
  expect_identical(quiet$value$spatial_niches, loud$spatial_niches)
})

test_that("findNiches stores one column per k and records metadata", {
  spe <- make_spatial_fixture()
  res <- findNiches(spe, n_niches = c(3, 2), use_pca = FALSE, verbose = FALSE)

  expect_true(all(c("spatial_niches_k2", "spatial_niches_k3") %in% colnames(SummarizedExperiment::colData(res))))
  meta <- S4Vectors::metadata(res)$nicheAnalysis[[1]]
  expect_equal(meta$n_niches_tested, c(2, 3))
  expect_named(meta$clustering_results, c("2", "3"))
})

test_that("findNiches records the variance explained by the PCs it used", {
  spe <- make_spatial_fixture()
  res <- findNiches(spe, n_niches = 3, ncomp = 3, verbose = FALSE)

  pca <- stats::prcomp(SingleCellExperiment::reducedDim(spe, "PhiSpace"))
  expected <- sum(pca$sdev[1:3]^2) / sum(pca$sdev^2)
  expect_equal(S4Vectors::metadata(res)$nicheAnalysis[[1]]$variance_explained, expected,
               tolerance = 1e-6)
})

test_that("findNiches validates its inputs", {
  spe <- make_spatial_fixture()
  expect_error(findNiches(spe, n_niches = 3, reducedDimName = "missing"), "not found")
  expect_error(findNiches(spe), "n_niches must be")
  expect_error(findNiches(matrix(1:4, 2), n_niches = 2), "SpatialExperiment or SingleCellExperiment")
})
