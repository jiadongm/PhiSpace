#' Principle component regression
#'
#' Same as `pls::svdpc.fit`, except using `rARPACK::svds` instead of `La.svd`.
#'
#' @param X Matrix. Scaled predictor matrix
#' @param Y Matrix. Scaled response matrix
#' @param ncomp Integer.
#'
#' @return A list containing
#' \item{coefficients}{Regression coefficient matrices.}
svdspc.fit <-
  function (X, Y, ncomp) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
    B <- array(0, dim = c(npred, nresp, ncomp))
    # This step may incur warnings if compute all singular values
    huhn <- suppressWarnings(rARPACK::svds(X, k = ncomp))
    D <- huhn$d[1:ncomp]
    TT <- huhn$u[, 1:ncomp, drop = FALSE] %*% diag(D, nrow = ncomp)
    P <- huhn$v[, 1:ncomp, drop = FALSE]
    tQ <- crossprod(TT, Y)/D^2
    for (a in 1:ncomp) {
      B[, , a] <- P[, 1:a, drop = FALSE] %*% tQ[1:a, ]
    }

    # Dimnames
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste0("comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(B) <- list(prednames, respnames, nCompnames)


    return(
      list(coefficients = B)
    )

  }
