make_dwd_data <- function(n, p = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p, dimnames = list(NULL, paste0("g", 1:p)))
  y <- factor(ifelse(X[, 1] + stats::rnorm(n) > 0, "a", "b"))
  list(X = X, y = y)
}

test_that("DWD with a non-linear kernel stops above max_cells before fitting", {
  d <- make_dwd_data(300)
  expect_error(
    rankFeatures(d$X, response = d$y, method = "DWD",
                 dwd_params = list(lambda = 0.1, kernel = kerndwd::rbfdot(sigma = 0.05),
                                   max_cells = 200)),
    "300 x 300 kernel matrix"
  )
  res <- rankFeatures(d$X, response = d$y, method = "DWD",
                      dwd_params = list(lambda = 0.1, kernel = kerndwd::rbfdot(sigma = 0.05)))
  expect_equal(nrow(res$importance_scores), ncol(d$X))
})

test_that("linear DWD works with at least as many cells as genes", {
  for (n in c(300, 20)) {                       # cells > genes, cells == genes
    d <- make_dwd_data(n)
    res <- rankFeatures(d$X, response = d$y, method = "DWD",
                        dwd_params = list(lambda = 0.1))
    beta <- res$importance_scores
    expect_equal(dim(beta), c(ncol(d$X), 1))
    expect_equal(rownames(beta), colnames(d$X))
    Xc <- scale(d$X, center = TRUE, scale = FALSE)
    expect_equal(unname(res$scores[, 1]),
                 unname(drop(Xc %*% beta) + res$model$alpha[1, 1]))
  }
  # g1 drives the class labels
  d <- make_dwd_data(300)
  res <- rankFeatures(d$X, response = d$y, method = "DWD",
                      dwd_params = list(lambda = 0.1))
  expect_equal(res$feature_ranking$feature[1], "g1")
})

test_that("linear DWD weights agree between the primal and dual fits", {
  # With fewer cells than genes kerndwd fits the dual; X' alpha must then
  # give the same kind of per-gene weight as the primal coefficients
  d <- make_dwd_data(15)
  res <- rankFeatures(d$X, response = d$y, method = "DWD",
                      dwd_params = list(lambda = 0.1))
  Xc <- scale(d$X, center = TRUE, scale = FALSE)
  expect_equal(unname(res$scores[, 1]),
               unname(drop(Xc %*% res$importance_scores) + res$model$alpha[1, 1]),
               tolerance = 1e-8)
})

test_that("PLS and PLSDA return the component scores", {
  skip_if_not_installed("pls")
  d <- make_dwd_data(60)
  y3 <- factor(rep(c("a", "b", "c"), 20))
  for (method in c("PLSDA", "PLS")) {
    resp <- if (method == "PLS") d$X[, 1] + stats::rnorm(60) else y3
    res <- rankFeatures(d$X[, -1], response = resp, method = method, ncomp = 3)
    expect_equal(dim(res$scores), c(60, 3))
    Y <- if (method == "PLS") resp else codeY_vec(resp)
    ref <- pls::kernelpls.fit(scale(d$X[, -1], scale = FALSE),
                              scale(as.matrix(Y), scale = FALSE), ncomp = 3)
    # Component signs are arbitrary
    expect_equal(abs(unname(res$scores)), abs(unname(unclass(ref$scores))),
                 tolerance = 1e-10)
  }
})
