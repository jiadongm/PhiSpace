#' Principal component analysis (PCA) based on partial singular value decomposition (SVD).
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
