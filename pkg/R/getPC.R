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
#' \item{sdev}{}
#' \item{totVar}{}
#' \item{props}{}
#' \item{accuProps}{}
#' \item{ncomp}{}
#' \item{selectFeat}{}
#' \item{reg_re}{}
#'
#' @export
getPC <- function(X, ncomp, center = TRUE, scale = FALSE, sparse = FALSE){

  if(center){
    Xmeans <- colMeans(X)
  } else {
    Xmeans <- FALSE
  }

  if(scale){
    Xscals <- apply(X, 2, stats::sd)
  } else {
    Xscals <- FALSE
  }


  X <- Matrix(
    scale(X, center = Xmeans, scale = Xscals),
    sparse = sparse
  )


  huhn <- suppressWarnings(rARPACK::svds(X, k = ncomp))
  D <- huhn$d
  if(sparse){
    huhn$u[huhn$u < 1e-15] <- 0
    huhn$v[huhn$v < 1e-15] <- 0
  }
  scores <- Matrix(
    huhn$u %*% diag(D, nrow = ncomp),
    sparse = sparse
  )
  loadings <- Matrix(
    huhn$v,
    sparse = sparse
  )
  sdev <- D/sqrt( nrow(X) - 1 )
  totVar <- sum(colSums(X^2)/(nrow(X)-1))
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
