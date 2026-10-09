fit_clusters <- function(X) {
  clusterPhiSpace(X, k_range = c(2, 5), ncomp = 5, algorithm = "Hartigan-Wong", seed = 1)
}

test_that("clusterPhiSpace selects k and recovers well-separated groups", {
  spe <- make_spatial_fixture()
  X <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  res <- fit_clusters(X)

  expect_s3_class(res, "PhiSpaceClustering")
  expect_equal(res$optimal_k, 3)
  expect_true(same_partition(res$clusters, spe$group))
})

test_that("clusterPhiSpace print and summary report the variance explained", {
  spe <- make_spatial_fixture()
  X <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  res <- fit_clusters(X)

  printed <- utils::capture.output(print(res))
  line <- grep("Variance explained by selected PCs:", printed, value = TRUE)
  expect_length(line, 1)
  expect_gt(as.numeric(sub(".*: *([0-9.]+) %.*", "\\1", line)), 0)
  expect_no_error(utils::capture.output(summary(res)))
})

test_that("clusterPhiSpace plots build for every type", {
  spe <- make_spatial_fixture()
  X <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  res <- fit_clusters(X)

  for (type in c("pca", "silhouette", "variance")) {
    p <- plot(res, type = type)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_match(plot(res, type = "pca")$labels$x, "^PC1 \\([0-9.]+%\\)$")
})

test_that("clusterPhiSpace computes all components without warnings", {
  spe <- make_spatial_fixture()
  X <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  expect_no_warning(res <- fit_clusters(X))

  ref <- stats::prcomp(X)
  expect_equal(res$pca_result$ncomp, ncol(X) - 1)
  expect_equal(res$pca_result$props, (ref$sdev^2 / sum(ref$sdev^2))[1:(ncol(X) - 1)])
  expect_equal(abs(unname(res$pc_scores)), abs(unname(ref$x[, 1:5])))
})

test_that("clusterPhiSpace computes silhouette widths on a subsample above silhouette_max_cells", {
  spe <- make_spatial_fixture()
  X <- SingleCellExperiment::reducedDim(spe, "PhiSpace")
  full <- clusterPhiSpace(X, k_range = c(2, 5), ncomp = 5, algorithm = "Hartigan-Wong", seed = 1)
  same <- clusterPhiSpace(X, k_range = c(2, 5), ncomp = 5, algorithm = "Hartigan-Wong", seed = 1,
                          silhouette_max_cells = nrow(X))
  expect_identical(full$k_selection, same$k_selection)
  expect_equal(full$k_selection$silhouette_cells, nrow(X))

  expect_message(
    sub <- clusterPhiSpace(X, k_range = c(2, 5), ncomp = 5, algorithm = "Hartigan-Wong", seed = 1,
                           silhouette_max_cells = 150),
    "random subsample of 150 of 300 cells"
  )
  expect_equal(sub$k_selection$silhouette_cells, 150)
  expect_length(sub$clusters, nrow(X))
  expect_equal(sub$optimal_k, 3)
  expect_output(summary(sub), "subsample of 150 cells")
})
