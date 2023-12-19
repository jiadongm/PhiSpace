#' Title
#'
#' Same as `pls::kernelpls.fit`, but replacing a full eigen decomposition by a partial one via `rARPACK::eigs_sym`.
#'
#' @param X
#' @param Y
#' @param ncomp
#' @param center
#' @param stripped
#'
#' @return
pls.fit <-
  function (X, Y, ncomp, center = TRUE, stripped = FALSE)
  {
    Y <- as.matrix(Y)
    if (!stripped) {
      dnX <- dimnames(X)
      dnY <- dimnames(Y)
    }
    dimnames(X) <- dimnames(Y) <- NULL
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
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
    R <- P <- matrix(0, ncol = ncomp, nrow = npred)
    tQ <- matrix(0, ncol = nresp, nrow = ncomp)
    B <- array(0, c(npred, nresp, ncomp))
    if (!stripped) {
      W <- P
      U <- TT <- matrix(0, ncol = ncomp, nrow = nobj)
      tsqs <- rep.int(1, ncomp)
      fitted <- array(0, c(nobj, nresp, ncomp))
    }
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
      if (!stripped) {
        tsqs[a] <- tsq
        u.a <- Y %*% q.a/c(crossprod(q.a))
        if (a > 1)
          u.a <- u.a - TT %*% (crossprod(TT, u.a)/tsqs)
        U[, a] <- u.a
        TT[, a] <- t.a
        W[, a] <- w.a
        fitted[, , a] <- TT[, 1:a] %*% tQ[1:a, , drop = FALSE]
      }
    }
    if (stripped) {
      list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
    }
    else {
      residuals <- -fitted + c(Y)
      if (center) {
        fitted <- fitted + rep(Ymeans, each = nobj)
      }
      objnames <- dnX[[1]]
      if (is.null(objnames))
        objnames <- dnY[[1]]
      prednames <- dnX[[2]]
      respnames <- dnY[[2]]
      compnames <- paste0("comp", 1:ncomp)
      nCompnames <- paste(1:ncomp, "comps")
      dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
      dimnames(R) <- dimnames(W) <- dimnames(P) <- list(prednames,
                                                        compnames)
      dimnames(tQ) <- list(compnames, respnames)
      dimnames(B) <- list(prednames, respnames, nCompnames)
      dimnames(fitted) <- dimnames(residuals) <- list(objnames,
                                                      respnames, nCompnames)
      # class(TT) <- class(U) <- "scores"
      # class(P) <- class(W) <- class(tQ) <- "loadings"
      list(coefficients = B, scores = TT, loadings = P, loading.weights = W,
           Yscores = U, Yloadings = t(tQ), projection = R, Xmeans = Xmeans,
           Ymeans = Ymeans, fitted.values = fitted, residuals = residuals,
           Xvar = colSums(P * P) * tsqs, Xtotvar = sum(X * X))
    }
  }
