#' Principal component analysis (PCA) based on partial singular value decomposition (SVD).
#'
#' @param X
#' @param ncomp
#' @param center
#' @param scale
#'
#' @return
#' @export
getPC <- function(X, ncomp, center = TRUE, scale = FALSE){
  Xmeans <- colMeans(X)

  X <- scale(X, center = Xmeans, scale = scale)
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
              reg_re = reg_re))
}
