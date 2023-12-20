#' Multivariate regression via principal component analysis (PCA) or partial least squares (PLS).
#'
#' Simplified version of `pls::mvr`, using computationally faster versions of PCA and PLS.
#'
#' @param X
#' @param Y
#' @param ncomp
#' @param method
#' @param center
#' @param scale
#'
#' @return A list containing
#' \item{coefficients}{Regression coefficient matrices.}
#' \item{scores}{Score matrix for X.}
#' \item{loadings}{Laoding matrix for X.}
#' \item{Yloadings}{Loading matrix for Y.}
#' \item{projection}{}
#' \item{Xmeans}{}
#' \item{Ymeans}{}
#' \item{fitted.values}{}
#' \item{residuals}{}
#' \item{Xvar}{}
#' \item{Xtotvar}{}
#' \item{ncomp}{}
#' \item{method}{}
#'
#' @export
mvr <- function(X, Y, ncomp, method = c("PCA", "PLS"), center = T, scale = F){

  if(scale){
    X <- scale(X, center = F, scale = T)
  }

  fitFunc <- switch(method,
                    PCA = svdspc.fit,
                    PLS = pls.fit
  )

  out <- fitFunc(X, Y, ncomp, center = center)
  out$ncomp <- ncomp
  out$method <- method

  return(out)
}
