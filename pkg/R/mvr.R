#' Multivariate regression via principal component analysis (PCA) or partial least squares (PLS).
#'
#' Simplified version of `pls::mvr`, using computationally faster versions of PCA and PLS.
#'
#' @param X Matrix.
#' @param Y Matrix.
#' @param ncomp Integer.
#' @param method Character.
#' @param center Logic.
#' @param sparse Convert X to sparse matrix or not.
#' @param scale Logic.
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.
#'
#' @return A list containing
#' \item{coefficients}{Regression coefficient matrices.}
#' \item{Xmeans}{}
#' \item{Ymeans}{}
#' \item{ncomp}{}
#' \item{method}{}
#'
#' @export
mvr <- function(
    X,
    Y,
    ncomp,
    method = c("PCA", "PLS"),
    center = TRUE,
    sparse = FALSE,
    scale = FALSE,
    DRinfo = FALSE
  ){


  Y <- scale(Y, center = TRUE, scale = FALSE)

  if(sparse){

    X <- Matrix(X, sparse = sparse)

  } else {

    X <- scale(X, center = center, scale = scale)
  }



  if(method == "PCA"){

    out <- svdspc.fit(X, Y, ncomp, sparse = sparse, DRinfo = DRinfo)
  } else {

    out <- pls.fit(X, Y, ncomp, DRinfo = DRinfo)
  }


  out$ncomp <- ncomp
  out$method <- method
  out$center <- center
  out$scale <- scale
  out$Xmeans <- attr(X, "scaled:center")
  out$Xscals <- attr(X, "scaled:scale")
  out$Ymeans <- attr(Y, "scaled:center")

  return(out)
}
