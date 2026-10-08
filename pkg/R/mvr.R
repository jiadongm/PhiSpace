#' Multivariate regression via principal component analysis (PCA) or partial least squares (PLS).
#'
#' Simplified version of `pls::mvr`, using computationally faster versions of PCA and PLS.
#'
#' A sparse `X` (for example a `dgCMatrix`) stays sparse: centring and
#' scaling are applied implicitly, without a dense copy of `X`. PLS uses the
#' kernel algorithm of `pls::kernelpls.fit`.
#'
#' @param X Matrix (dense or sparse).
#' @param Y Matrix.
#' @param ncomp Integer.
#' @param method Character.
#' @param center Logic.
#' @param sparse Not used; kept for backward compatibility.
#' @param scale Logic.
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.
#' @param keepComps Integer vector. Numbers of components for which to return
#'   regression coefficients. The default returns all, `1:ncomp`. Use
#'   `keepComps = ncomp` to return only the final coefficients and save memory.
#'
#' @return A list containing
#' \item{coefficients}{Array of regression coefficients with dimensions
#'   features x responses x `length(keepComps)`. The third dimension is named
#'   `"<k> comps"`, so `coefficients[, , paste(ncomp, "comps")]` selects the
#'   coefficients for `ncomp` components.}
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
    DRinfo = FALSE,
    keepComps = seq_len(ncomp)
  ){

  method <- match.arg(method)

  keepComps <- sort(unique(as.integer(keepComps)))
  if(length(keepComps) == 0 || any(is.na(keepComps)) ||
     any(keepComps < 1) || any(keepComps > ncomp)){
    stop("keepComps must contain integers between 1 and ncomp.")
  }

  # Make sure X and Y have dims, can be used as argument of eg colMeans.
  # A sparse X stays sparse.
  X <- .fit_matrix(X)
  Y <- as.matrix(Y)

  if(method == "PCA"){

    out <- svdspc.fit(X, Y, ncomp, center = center, scale = scale, sparse = sparse, DRinfo = DRinfo,
                      keepComps = keepComps)
  } else {

    out <- pls.fit(X, Y, ncomp, center = center, scale = scale, DRinfo = DRinfo,
                   keepComps = keepComps)
  }


  out$ncomp <- ncomp
  out$method <- method
  out$center <- center
  out$scale <- scale
  out$Xmeans <- out$Xmeans
  out$Xscals <- out$Xscals
  out$Ymeans <- out$Ymeans
  out$Yscals <- out$Yscals

  return(out)
}
