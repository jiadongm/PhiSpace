#' Title
#'
#' Same as `pls::kernelpls.fit`, but replacing a full eigen decomposition by a partial one via `rARPACK::eigs_sym`.
#'
#' @param X Matrix. Scaled predictor matrix
#' @param Y Matrix. Scaled response matrix
#' @param ncomp Integer.
#'
#' @return A list containing
#' \item{coefficients}{Regression coefficient matrices.}
pls.fit <-
  function (X, Y, ncomp)
  {

    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]

    R <- P <- matrix(0, ncol = ncomp, nrow = npred)
    tQ <- matrix(0, ncol = nresp, nrow = ncomp)
    B <- array(0, c(npred, nresp, ncomp))
    XtY <- crossprod(X, Y)
    for (a in 1:ncomp) {
      if (nresp == 1) {
        w.a <- XtY/sqrt(c(crossprod(XtY)))
      }
      else {
        if (nresp < npred) {
          q <- rARPACK::eigs_sym(crossprod(XtY), k = 1)$vectors[,1]
          w.a <- XtY %*% q
          w.a <- w.a/sqrt(c(crossprod(w.a)))
        }
        else {
          w.a <- rARPACK::eigs_sym(crossprod(XtY), k = 1)$vectors[,1]
        }
      }
      r.a <- w.a
      if (a > 5) {
        r.a <- r.a - colSums(crossprod(w.a, P[, 1:(a - 1),
                                              drop = FALSE]) %*% t(R[, 1:(a - 1), drop = FALSE]))
      }
      else if (a > 1) {
        for (j in 1:(a - 1)) r.a <- r.a - c(P[, j] %*% w.a) *
            R[, j]
      }
      t.a <- X %*% r.a
      tsq <- c(crossprod(t.a))
      p.a <- crossprod(X, t.a)/tsq
      q.a <- crossprod(XtY, r.a)/tsq
      XtY <- XtY - (tsq * p.a) %*% t(q.a)
      R[, a] <- r.a
      P[, a] <- p.a
      tQ[a, ] <- q.a

      B[, , a] <- R[, 1:a, drop = FALSE] %*% tQ[1:a, , drop = FALSE]

    }

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
