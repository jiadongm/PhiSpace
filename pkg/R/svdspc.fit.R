svdspc.fit <- function (X, Y, ncomp, center = TRUE, scale = FALSE,sparse = FALSE, DRinfo = FALSE) {

    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
    B <- array(0, dim = c(npred, nresp, ncomp))

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

    # This step may incur warnings if compute all singular values
    huhn <- suppressWarnings(irlba::irlba(X, nv = ncomp, center = Xmeans, scale = Xscals))
    D <- huhn$d

    # if(sparse){
    #   huhn$u[huhn$u < 1e-15] <- 0
    #   huhn$v[huhn$v < 1e-15] <- 0
    # }

    TT <- huhn$u %*% diag(D, nrow = ncomp)
    P <- huhn$v
    dimnames(TT) <- list(
      dnX[[1]],
      paste0("comp", 1:ncomp)
    )
    dimnames(P) <- list(
      dnX[[2]],
      paste0("comp", 1:ncomp)
    )
    tQ <- crossprod(TT, Y)/D^2
    for (a in 1:ncomp) {
      B[, , a] <- as.matrix(P[, 1:a, drop = FALSE] %*% tQ[1:a, ])
    }

    # Dimnames
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
