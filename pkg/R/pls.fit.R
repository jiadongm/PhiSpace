pls.fit <-
  function (X, Y, ncomp, center = TRUE, scale = FALSE, DRinfo = FALSE)
  {

    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]

    if(center){

      Xmeans <- colMeans(X)
      Ymeans <- colMeans(Y)
    } else {

      Xmeans <- NULL
      Ymeans <- NULL
    }

    if(scale){

      Xscals <- apply(X, 2, stats::sd)
      Yscals <- apply(Y, 2, stats::sd)
    } else {

      Xscals <- NULL
      Yscals <- NULL
    }

    X <- scal(
      X,
      center = Xmeans,
      scale = Xscals
    )

    Y <- scal(
      Y,
      center = Ymeans,
      scale = Yscals
    )

    if(DRinfo){
      TT <- matrix(0, ncol = ncomp, nrow = nobj)
      dimnames(TT) <- list(
        dnX[[1]],
        paste0("comp", 1:ncomp)
      )
    }

    R <- P <- matrix(0, ncol = ncomp, nrow = npred)
    tQ <- matrix(0, ncol = nresp, nrow = ncomp)
    B <- array(0, c(npred, nresp, ncomp))
    XtY <- crossprod(X, Y)
    for (a in 1:ncomp) {

      if (nresp == 1) {

        w.a <- XtY/sqrt(as.numeric(crossprod(XtY)))
      } else {

        if (nresp < npred) {
          q <- irlba::partial_eigen(crossprod(XtY), n = 1)$vectors[,1]
          w.a <- XtY %*% q
          w.a <- w.a/sqrt(as.numeric(crossprod(w.a)))
        } else {
          w.a <- irlba::partial_eigen(tcrossprod(XtY), n = 1)$vectors[,1]
        }

      }
      r.a <- w.a
      if (a > 5) {

        r.a <- r.a - colSums(crossprod(w.a, P[, 1:(a - 1), drop = FALSE]) %*% t(R[, 1:(a - 1), drop = FALSE]))

      } else if (a > 1) {

        for (j in 1:(a - 1)) r.a <- r.a - as.numeric(P[, j] %*% w.a) * R[, j]

      }

      t.a <- X %*% r.a
      tsq <- as.numeric(crossprod(t.a))
      p.a <- crossprod(X, t.a)/tsq
      q.a <- crossprod(XtY, r.a)/tsq
      XtY <- XtY - (tsq * p.a) %*% t(q.a)
      R[, a] <- as.vector(r.a)
      P[, a] <- as.vector(p.a)
      tQ[a, ] <- as.vector(q.a)

      B[, , a] <- R[, 1:a, drop = FALSE] %*% tQ[1:a, , drop = FALSE]


      if (DRinfo) {
        TT[, a] <- t.a
        dimnames(P) <- list(
          dnX[[2]],
          paste0("comp", 1:ncomp)
        )
      }

    }

    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste0("comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(B) <- list(prednames, respnames, nCompnames)


    if(DRinfo){

      return(
        list(
          coefficients = B,
          scores = TT,
          loadings = P,
          Xmeans = Xmeans,
          Ymeans = Ymeans,
          Xscals = Xscals,
          Yscals = Yscals
        )
      )

    } else {

      return(
        list(
          coefficients = B,
          Xmeans = Xmeans,
          Ymeans = Ymeans,
          Xscals = Xscals,
          Yscals = Yscals
        )
      )

    }

  }
