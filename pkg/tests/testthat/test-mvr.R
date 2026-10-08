sim_xy <- function(n = 300, p = 60, K = 5, seed = 1) {
  set.seed(seed)
  X <- Matrix::rsparsematrix(n, p, density = 0.2,
                             rand.x = function(m) log1p(stats::rpois(m, 3) + 1))
  dimnames(X) <- list(paste0("cell", seq_len(n)), paste0("gene", seq_len(p)))
  lab <- sample(paste0("type", seq_len(K)), n, replace = TRUE)
  list(X = X, Y = codeY_vec(lab, rownames(X)))
}

sd_scale <- function(M, center, scale) {
  base::scale(M, center = center,
              scale = if (scale) apply(M, 2, stats::sd) else FALSE)
}

test_that("PLS coefficients match pls::kernelpls.fit", {
  skip_if_not_installed("pls")
  d <- sim_xy()
  K <- ncol(d$Y)
  for (center in c(TRUE, FALSE)) for (scale in c(FALSE, TRUE)) {
    fit <- mvr(d$X, d$Y, ncomp = K, method = "PLS", center = center, scale = scale)
    ref <- pls::kernelpls.fit(sd_scale(as.matrix(d$X), center, scale),
                              sd_scale(d$Y, center, scale),
                              ncomp = K, center = center)
    expect_equal(unname(fit$coefficients), unname(ref$coefficients),
                 tolerance = 1e-10)
  }
})

test_that("PLS with one response matches pls::kernelpls.fit", {
  skip_if_not_installed("pls")
  d <- sim_xy()
  y <- d$Y[, 1, drop = FALSE]
  fit <- mvr(d$X, y, ncomp = 4, method = "PLS")
  ref <- pls::kernelpls.fit(sd_scale(as.matrix(d$X), TRUE, FALSE),
                            sd_scale(y, TRUE, FALSE), ncomp = 4)
  expect_equal(unname(fit$coefficients), unname(ref$coefficients), tolerance = 1e-10)
})

test_that("PLS with more responses than features matches pls::kernelpls.fit", {
  skip_if_not_installed("pls")
  d <- sim_xy(p = 6, K = 8)
  fit <- mvr(d$X, d$Y, ncomp = 4, method = "PLS")
  ref <- pls::kernelpls.fit(sd_scale(as.matrix(d$X), TRUE, FALSE),
                            sd_scale(d$Y, TRUE, FALSE), ncomp = 4)
  expect_equal(unname(fit$coefficients), unname(ref$coefficients), tolerance = 1e-10)
})

test_that("PLS is deterministic and gives the same result for sparse and dense X", {
  d <- sim_xy()
  set.seed(1)
  a <- mvr(d$X, d$Y, ncomp = 5, method = "PLS", scale = TRUE)
  set.seed(2)
  b <- mvr(d$X, d$Y, ncomp = 5, method = "PLS", scale = TRUE)
  dense <- mvr(as.matrix(d$X), d$Y, ncomp = 5, method = "PLS", scale = TRUE)
  expect_identical(a$coefficients, b$coefficients)
  expect_equal(a$coefficients, dense$coefficients, tolerance = 1e-12)
  expect_equal(a$Xscals, apply(as.matrix(d$X), 2, stats::sd), tolerance = 1e-12)
})

test_that("PCA regression matches a regression on prcomp scores", {
  d <- sim_xy()
  for (ncomp in c(5, 40)) {   # 40 uses the full SVD branch (>= half of 60 features)
    fit <- mvr(d$X, d$Y, ncomp = ncomp, method = "PCA")
    Xc <- base::scale(as.matrix(d$X), scale = FALSE)
    pc <- stats::prcomp(Xc, center = FALSE)
    TT <- pc$x[, seq_len(ncomp)]
    ref <- pc$rotation[, seq_len(ncomp)] %*% solve(crossprod(TT), crossprod(TT, d$Y))
    expect_equal(unname(fit$coefficients[, , ncomp]), unname(ref), tolerance = 1e-6)
  }
})

test_that("keepComps keeps the requested coefficient slices", {
  d <- sim_xy()
  for (method in c("PLS", "PCA")) {
    set.seed(1)
    full <- mvr(d$X, d$Y, ncomp = 5, method = method)
    set.seed(1)
    part <- mvr(d$X, d$Y, ncomp = 5, method = method, keepComps = c(5, 2))
    expect_equal(dim(full$coefficients), c(60, 5, 5))
    expect_equal(dim(part$coefficients), c(60, 5, 2))
    expect_equal(dimnames(part$coefficients)[[3]], c("2 comps", "5 comps"))
    expect_equal(part$coefficients[, , "5 comps"], full$coefficients[, , 5])
    expect_equal(part$coefficients[, , "2 comps"], full$coefficients[, , 2])
  }
  expect_error(mvr(d$X, d$Y, ncomp = 5, method = "PLS", keepComps = 6), "keepComps")
})

test_that("phenotype accepts full and final-slice coefficient arrays", {
  d <- sim_xy()
  full <- mvr(d$X, d$Y, ncomp = 5, method = "PLS")
  final <- mvr(d$X, d$Y, ncomp = 5, method = "PLS", keepComps = 5)
  atlas <- function(fit) list(ncomp = 5, selectedFeat = colnames(d$X),
                              center = TRUE, scale = FALSE, reg_re = fit)
  a <- phenotype(d$X, atlas(full), assayName = "log1p")
  b <- phenotype(d$X, atlas(final), assayName = "log1p")
  expect_equal(a, b)
  expect_equal(dimnames(a), list(rownames(d$X), colnames(d$Y)))
})

test_that("phenotype matches explicit centring and scaling", {
  d <- sim_xy()
  Xd <- as.matrix(d$X)
  Xq <- d$X[1:120, ]
  Xqd <- as.matrix(Xq)
  for (center in c(TRUE, FALSE)) for (scale in c(FALSE, TRUE)) {
    fit <- mvr(d$X, d$Y, ncomp = 5, method = "PLS", center = center, scale = scale)
    atlas <- list(ncomp = 5, selectedFeat = colnames(d$X),
                  center = center, scale = scale, reg_re = fit)
    B <- fit$coefficients[, , 5]
    Ym <- if (center) colMeans(d$Y) else 0

    byQuery <- base::scale(Xqd, center = center, scale = scale) %*% B
    expect_equal(phenotype(Xq, atlas, "log1p", "byQuery"),
                 sweep(byQuery, 2, Ym, "+"), tolerance = 1e-12,
                 ignore_attr = TRUE)

    byRef <- base::scale(Xqd, center = if (center) colMeans(Xd) else FALSE,
                         scale = if (scale) fit$Xscals else FALSE) %*% B
    expect_equal(phenotype(Xq, atlas, "log1p", "byRef"),
                 sweep(byRef, 2, Ym, "+"), tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
})
