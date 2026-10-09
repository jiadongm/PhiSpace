#' Principal component analysis (PCA) based on partial singular value decomposition (SVD).
#'
#' `getPC()` uses irlba's partial SVD. When `ncomp` is at least half of the
#' smaller dimension of `X`, where irlba can be inaccurate, it uses a full SVD
#' of a dense copy of `X` instead.
#'
#' @param X Matrix.
#' @param ncomp Integer.
#' @param center Logic.
#' @param scale Logic.
#' @param sparse Use sparse matrix or not.
#'
#' @return A list containing
#' \item{scores}{Score matrix for X.}
#' \item{loadings}{Laoding matrix for X.}
#' \item{sdev}{Standard deviations of the principal components.}
#' \item{totVar}{Total variance of `X` after the requested centring and scaling.}
#' \item{props}{Proportion of `totVar` explained by each component.}
#' \item{accuProps}{Cumulative sum of `props`.}
#' \item{ncomp}{Number of components.}
#' \item{Xmeans}{Column means used for centring, or `NULL`.}
#' \item{Xscals}{Column standard deviations used for scaling, or `NULL`.}
#'
#' @export
getPC <- function(X, ncomp, center = TRUE, scale = FALSE, sparse = FALSE){

  if(ncomp >= 0.5 * min(dim(X))){
    # irlba warns here and may be inaccurate. One dimension of X is at most
    # 2 * ncomp, so a dense copy is small.
    return(.getPC_svd(X, ncomp, center = center, scale = scale))
  }

  if(center){
    Xmeans <- colMeans(X)
  } else {
    Xmeans <- NULL
  }

  if(scale){
    Xscals <- apply(X, 2, stats::sd)
  } else {
    Xscals <- NULL
  }

  irlba_res <- irlba::irlba(
    X,
    nv = ncomp,
    center = Xmeans,
    scale = Xscals
  )

  loadings <- irlba_res$v
  scores <- irlba_res$u %*% diag(irlba_res$d)

  sdev <- irlba_res$d/sqrt( nrow(X) - 1 )
  # Total variance of the matrix that irlba decomposed (after centring and
  # scaling), computed without densifying a sparse X.
  colSS <- colSums(X^2)
  if(center) colSS <- colSS - nrow(X) * Xmeans^2
  if(scale) colSS <- colSS / Xscals^2
  totVar <- sum(colSS)/(nrow(X)-1)
  props <- sdev^2/totVar
  accuProps <- cumsum(props)


  rownames(scores) <- rownames(X)
  colnames(scores) <- paste0('comp', 1:ncol(scores))
  rownames(loadings) <- colnames(X)
  colnames(loadings) <- paste0('comp', 1:ncol(scores))

  return(list(
    scores = scores,
    loadings = loadings,
    sdev = sdev,
    totVar = totVar,
    props = props,
    accuProps = accuProps,
    ncomp = ncomp,
    Xmeans = Xmeans,
    Xscals = Xscals
  ))
}


# PCA by a full SVD of the dense matrix. It returns the same elements as getPC().
# Use it when nearly all components are needed: irlba is designed for a few
# leading components, and for most of them it warns and may not converge.
.getPC_svd <- function(X, ncomp, center = TRUE, scale = FALSE){

  X <- as.matrix(X)
  ncomp <- min(ncomp, dim(X))
  Xmeans <- if(center) colMeans(X) else NULL
  Xscals <- if(scale) apply(X, 2, stats::sd) else NULL
  Xc <- base::scale(X,
                    center = if(center) Xmeans else FALSE,
                    scale = if(scale) Xscals else FALSE)

  svd_res <- svd(Xc, nu = ncomp, nv = ncomp)
  d <- svd_res$d[seq_len(ncomp)]

  loadings <- svd_res$v
  scores <- svd_res$u %*% diag(d, nrow = ncomp)

  sdev <- d/sqrt( nrow(X) - 1 )
  totVar <- sum(Xc^2)/(nrow(X)-1)
  props <- sdev^2/totVar
  accuProps <- cumsum(props)

  rownames(scores) <- rownames(X)
  colnames(scores) <- paste0('comp', 1:ncomp)
  rownames(loadings) <- colnames(X)
  colnames(loadings) <- paste0('comp', 1:ncomp)

  return(list(
    scores = scores,
    loadings = loadings,
    sdev = sdev,
    totVar = totVar,
    props = props,
    accuProps = accuProps,
    ncomp = ncomp,
    Xmeans = Xmeans,
    Xscals = Xscals
  ))
}
