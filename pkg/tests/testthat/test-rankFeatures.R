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
