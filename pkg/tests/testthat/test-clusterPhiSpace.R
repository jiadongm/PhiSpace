# clusterPhiSpace() computes ncol(x) - 1 components with irlba whatever ncomp is,
# so irlba warns that it computes too many singular values and may not converge
# on the last ones. Only those two known warnings are muffled here.
fit_clusters <- function(X) {
  withCallingHandlers(
    clusterPhiSpace(X, k_range = c(2, 5), ncomp = 5, algorithm = "Hartigan-Wong", seed = 1),
    warning = function(w) {
      if (grepl("too large a percentage|did not converge", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
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
