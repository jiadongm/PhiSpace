#' Principal component analysis (PCA) based on partial singular value decomposition (SVD).
#'
#' @param X Matrix.
#' @param ncomp Integer.
#' @param center Logic.
#' @param scale Logic.
#'
#' @return A list containing
#' \item{scores}{Score matrix for X.}
#' \item{loadings}{Laoding matrix for X.}
#' \item{sdev}{}
#' \item{totVar}{}
#' \item{props}{}
#' \item{accuProps}{}
#' \item{ncomp}{}
#' \item{selectFeat}{}
#' \item{reg_re}{}
#'
#' @export
getPC <- function(X, ncomp, center = TRUE, scale = FALSE){

  X <- scale(X, center = center, scale = scale)

  Xmeans <- attr(X, "scaled:center")
  Xscals <- attr(X, "scaled:scale")

  huhn <- suppressWarnings(rARPACK::svds(X, k = ncomp))
  D <- huhn$d[1:ncomp]
  TT <- huhn$u[, 1:ncomp, drop = FALSE] %*% diag(D, nrow = ncomp)
  P <- huhn$v[, 1:ncomp, drop = FALSE]
  sdev <- D/sqrt( nrow(X) - 1 )
  totVar <- sum(colSums(X^2)/(nrow(X)-1))
  props <- sdev^2/totVar
  accuProps <- cumsum(props)

  scores <- TT
  rownames(scores) <- rownames(X)
  colnames(scores) <- paste0('comp', 1:ncol(scores))
  loadings <- P
  rownames(loadings) <- colnames(X)
  colnames(loadings) <- paste0('comp', 1:ncol(scores))

  coefficients <- array(NA, c(dim(P), ncomp))
  coefficients[,,ncomp] <- P

  # Store results in reg_re, mimicking supervised mode; see SuperPC.R
  reg_re <- list(Xmeans = Xmeans,
                 Xscals = Xscals,
                 Ymeans = rep(0, ncol(P)),
                 coefficients = coefficients,
                 loadings = P
  )

  return(list(scores = scores,
              loadings = loadings,
              sdev = sdev,
              totVar = totVar,
              props = props,
              accuProps = accuProps,
              ncomp = ncomp,
              selectFeat = colnames(X),
              reg_re = reg_re,
              center = center,
              scale = scale))
}
