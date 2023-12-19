#' Principle component regression
#'
#' Same as `pls::svdpc.fit`, except using `rARPACK::svds` instead of `La.svd`.
#'
#' @param X
#' @param Y
#' @param ncomp
#' @param center
#' @param stripped
#'
#' @return
svdspc.fit <-
  function (X, Y, ncomp, center = TRUE, stripped = FALSE) {
    Y <- as.matrix(Y)
    if (!stripped) {
      dnX <- dimnames(X)
      dnY <- dimnames(Y)
    }
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
    B <- array(0, dim = c(npred, nresp, ncomp))
    if (!stripped)
      fitted <- array(0, dim = c(nobj, nresp, ncomp))
    if (center) {
      Xmeans <- colMeans(X)
      X <- X - rep(Xmeans, each = nobj)
      Ymeans <- colMeans(Y)
      Y <- Y - rep(Ymeans, each = nobj)
    }
    else {
      Xmeans <- rep_len(0, npred)
      Ymeans <- rep_len(0, nresp)
    }
    # This step may incur warnings if compute all singular values
    huhn <- suppressWarnings(rARPACK::svds(X, k = ncomp))
    D <- huhn$d[1:ncomp]
    TT <- huhn$u[, 1:ncomp, drop = FALSE] %*% diag(D, nrow = ncomp)
    P <- huhn$v[, 1:ncomp, drop = FALSE]
    tQ <- crossprod(TT, Y)/D^2
    for (a in 1:ncomp) {
      B[, , a] <- P[, 1:a, drop = FALSE] %*% tQ[1:a, ]
      if (!stripped)
        fitted[, , a] <- TT[, 1:a, drop = FALSE] %*% tQ[1:a,
        ]
    }
    if (stripped) {
      list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
    }
    else {
      residuals <- c(Y) - fitted
      fitted <- fitted + rep(Ymeans, each = nobj)
      objnames <- dnX[[1]]
      if (is.null(objnames))
        objnames <- dnY[[1]]
      prednames <- dnX[[2]]
      respnames <- dnY[[2]]
      compnames <- paste0("comp", 1:ncomp)
      nCompnames <- paste(1:ncomp, "comps")
      dimnames(TT) <- list(objnames, compnames)
      dimnames(P) <- list(prednames, compnames)
      dimnames(tQ) <- list(compnames, respnames)
      dimnames(B) <- list(prednames, respnames, nCompnames)
      dimnames(fitted) <- dimnames(residuals) <- list(objnames,
                                                      respnames, nCompnames)
      names(D) <- compnames
      # class(TT) <- "scores"
      R <- P
      # class(P) <- class(tQ) <- "loadings"
      list(coefficients = B, scores = TT, loadings = P, Yloadings = t(tQ),
           projection = R, Xmeans = Xmeans, Ymeans = Ymeans,
           fitted.values = fitted, residuals = residuals, Xvar = D^2,
           Xtotvar = sum(X * X))
    }
  }
