test_that("getPC variance proportions match prcomp for centred data", {
  set.seed(1)
  X <- matrix(rnorm(500 * 40, mean = 2), 500, 40) %*% diag(1:40)
  pc <- getPC(X, ncomp = 5, center = TRUE, scale = FALSE)
  ref <- stats::prcomp(X, center = TRUE, scale. = FALSE)

  expect_equal(pc$totVar, sum(ref$sdev^2))
  expect_equal(pc$props, (ref$sdev^2 / sum(ref$sdev^2))[1:5])
  expect_equal(pc$accuProps, cumsum(pc$props))
  expect_equal(abs(unname(pc$scores)), abs(unname(ref$x[, 1:5])), tolerance = 1e-6)
})

test_that("getPC variance proportions match prcomp for scaled data", {
  set.seed(2)
  X <- matrix(rnorm(500 * 40, mean = 2), 500, 40) %*% diag(1:40)
  pc <- getPC(X, ncomp = 5, center = TRUE, scale = TRUE)
  ref <- stats::prcomp(X, center = TRUE, scale. = TRUE)

  expect_equal(pc$totVar, ncol(X))
  expect_equal(pc$props, (ref$sdev^2 / sum(ref$sdev^2))[1:5])
})

test_that("getPC gives the same total variance for sparse and dense input", {
  set.seed(3)
  Xs <- Matrix::rsparsematrix(300, 30, density = 0.2)
  sparse <- getPC(Xs, ncomp = 5)
  dense <- getPC(as.matrix(Xs), ncomp = 5)

  expect_equal(sparse$totVar, dense$totVar)
  expect_equal(sparse$props, dense$props, tolerance = 1e-6)
})
