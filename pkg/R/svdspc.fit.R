svdspc.fit <- function (X, Y, ncomp, center = TRUE, scale = FALSE, sparse = FALSE, DRinfo = FALSE,
                        keepComps = seq_len(ncomp)) {

    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
    B <- array(0, dim = c(npred, nresp, length(keepComps)))

    if(center){

      Xmeans <- colMeans(X)
      Ymeans <- colMeans(Y)
    } else {

      Xmeans <- NULL
      Ymeans <- NULL
    }

    if(scale){

      Xscals <- .colSds(X)
      Yscals <- apply(Y, 2, stats::sd)
    } else {

      Xscals <- NULL
      Yscals <- NULL
    }

    if(ncomp < 0.5 * min(nobj, npred)){

      # irlba centres and scales a sparse X implicitly
      huhn <- suppressWarnings(irlba::irlba(X, nv = ncomp, center = Xmeans, scale = Xscals))
    } else {

      # irlba is not accurate for most of the singular values. Here one
      # dimension of X is at most 2 * ncomp, so a dense copy is small.
      huhn <- svd(scal(as.matrix(X), center = Xmeans, scale = Xscals), nu = ncomp, nv = ncomp)
      huhn$d <- huhn$d[seq_len(ncomp)]
    }
    D <- huhn$d

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
    Bcum <- matrix(0, nrow = npred, ncol = nresp)
    for (a in 1:ncomp) {
      Bcum <- Bcum + P[, a] %o% tQ[a, ]
      k <- match(a, keepComps)
      if (!is.na(k)) B[, , k] <- Bcum
    }

    # Dimnames
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    dimnames(B) <- list(prednames, respnames, paste(keepComps, "comps"))


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
