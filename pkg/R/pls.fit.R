# Kernel PLS (as in pls::kernelpls.fit). X is centred and scaled implicitly:
# the code computes products with (X - 1 Xmeans') diag(1/Xscals) without
# forming that matrix, so a sparse X stays sparse.
pls.fit <-
  function (X, Y, ncomp, center = TRUE, scale = FALSE, DRinfo = FALSE,
            keepComps = seq_len(ncomp))
  {

    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]

    if (center) {
      Xmeans <- colMeans(X)
      Ymeans <- colMeans(Y)
    } else {
      Xmeans <- NULL
      Ymeans <- NULL
    }

    if (scale) {
      Xscals <- .colSds(X)
      Yscals <- apply(Y, 2, stats::sd)
    } else {
      Xscals <- NULL
      Yscals <- NULL
    }

    Y <- scal(Y, center = Ymeans, scale = Yscals)

    mu <- if (center) Xmeans else rep(0, npred)
    s <- if (scale) Xscals else rep(1, npred)
    # Products with the centred and scaled X
    Xr <- function(r) {
      r <- r/s
      as.numeric(X %*% r) - sum(mu * r)
    }
    Xtt <- function(tt) (as.numeric(crossprod(X, tt)) - mu * sum(tt))/s

    XtY <- (as.matrix(crossprod(X, Y)) - outer(mu, colSums(Y)))/s

    # Scores t = X r and loadings numerator X't for a weight vector r
    project <- function(r) {
      t <- Xr(r)
      list(t = t, tsq = sum(t^2), Xtt = Xtt(t))
    }
    fit <- .kernelpls(XtY, project, ncomp, keepComps, nobj = if (DRinfo) nobj)
    B <- fit$B
    P <- fit$P
    TT <- fit$TT
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste0("comp", 1:ncomp)
    dimnames(B) <- list(prednames, respnames, paste(keepComps, "comps"))

    if (DRinfo) {
      dimnames(TT) <- list(dnX[[1]], compnames)
      dimnames(P) <- list(prednames, compnames)
      return(list(coefficients = B, scores = TT, loadings = P,
                  Xmeans = Xmeans, Ymeans = Ymeans, Xscals = Xscals,
                  Yscals = Yscals))
    }
    else {
      return(list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans,
                  Xscals = Xscals, Yscals = Yscals))
    }

  }


# The kernel PLS loop. XtY is X'Y for the centred and scaled X and Y.
# project(r) returns list(t, tsq, Xtt): the scores t = X r (or NULL), their
# sum of squares t't = r'X'X r, and X't = X'X r. Scores are kept when nobj is
# not NULL. Returns the coefficients B (npred x nresp x length(keepComps)),
# the loadings P and the scores TT.
.kernelpls <- function(XtY, project, ncomp, keepComps, nobj = NULL)
{
  npred <- nrow(XtY)
  nresp <- ncol(XtY)
  keepTT <- !is.null(nobj)
  if (keepTT) TT <- matrix(0, nrow = nobj, ncol = ncomp)
  R <- P <- matrix(0, nrow = npred, ncol = ncomp)
  B <- array(0, c(npred, nresp, length(keepComps)))
  Bcum <- matrix(0, nrow = npred, ncol = nresp)

  for (a in 1:ncomp) {
    if (nresp == 1) {
      w.a <- XtY/sqrt(sum(XtY^2))
    } else {
      if (nresp < npred) {
        q <- eigen(crossprod(XtY), symmetric = TRUE)$vectors[, 1]
        w.a <- XtY %*% q
        w.a <- w.a/sqrt(sum(w.a^2))
      } else {
        w.a <- eigen(tcrossprod(XtY), symmetric = TRUE)$vectors[, 1]
      }
    }
    w.a <- as.numeric(w.a)
    r.a <- w.a
    if (a > 1) {
      r.a <- r.a - as.numeric(R[, 1:(a - 1), drop = FALSE] %*%
                                crossprod(P[, 1:(a - 1), drop = FALSE], w.a))
    }
    pr <- project(r.a)
    tsq <- pr$tsq
    p.a <- pr$Xtt/tsq
    q.a <- as.numeric(crossprod(XtY, r.a))/tsq
    XtY <- XtY - (tsq * p.a) %o% q.a
    R[, a] <- r.a
    P[, a] <- p.a
    Bcum <- Bcum + r.a %o% q.a
    k <- match(a, keepComps)
    if (!is.na(k)) B[, , k] <- Bcum
    if (keepTT) TT[, a] <- pr$t
  }

  list(B = B, P = P, TT = if (keepTT) TT)
}
