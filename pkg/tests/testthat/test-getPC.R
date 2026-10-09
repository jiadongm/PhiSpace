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

test_that("the SVD fallback of getPC matches prcomp and getPC", {
  set.seed(4)
  X <- matrix(rnorm(200 * 12, mean = 1), 200, 12) %*% diag(1:12)
  for (scale in c(FALSE, TRUE)) {
    full <- .getPC_svd(X, ncomp = 11, center = TRUE, scale = scale)
    ref <- stats::prcomp(X, center = TRUE, scale. = scale)
    expect_equal(full$props, (ref$sdev^2 / sum(ref$sdev^2))[1:11])
    expect_equal(abs(unname(full$scores)), abs(unname(ref$x[, 1:11])))

    few <- getPC(X, ncomp = 3, center = TRUE, scale = scale)
    expect_equal(full$totVar, few$totVar)
    expect_equal(abs(full$scores[, 1:3]), abs(few$scores), tolerance = 1e-6)
    expect_named(full, names(few))
  }
  expect_equal(.getPC_svd(X, ncomp = 50)$ncomp, 12)
})

test_that("getPC uses a full SVD when ncomp is at least half of min(dim(X))", {
  set.seed(5)
  X <- matrix(rnorm(300 * 20, mean = 1), 300, 20) %*% diag(1:20)
  ref <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  expect_no_warning(pc <- getPC(X, ncomp = 10))
  expect_equal(pc$props, (ref$sdev^2 / sum(ref$sdev^2))[1:10])
  expect_equal(abs(unname(pc$scores)), abs(unname(ref$x[, 1:10])))
  expect_equal(pc$ncomp, 10)
})

test_that("getPC scales sparse and dense input the same way", {
  set.seed(6)
  Xs <- Matrix::rsparsematrix(300, 30, density = 0.2)
  sparse <- getPC(Xs, ncomp = 5, scale = TRUE)
  dense <- getPC(as.matrix(Xs), ncomp = 5, scale = TRUE)
  expect_equal(sparse$Xscals, dense$Xscals)
  expect_equal(sparse$totVar, ncol(Xs))
  expect_equal(abs(sparse$scores), abs(dense$scores), tolerance = 1e-6)
})
