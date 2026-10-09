#' Soft phenotyping of query assay.
#'
#' @param phenoAssay Matrix. Query assay to phenotype.
#' @param atlas_re List.
#' @param assayName Character.
#' @param scaleMethod Character.
#'
#' @return A list containing
#' \item{Yhat}{}
#' \item{Proj}{}
phenotype <- function(phenoAssay,
                      atlas_re,
                      assayName = 'rank',
                      scaleMethod = c("byQuery", "byRef")){


  scaleMethod <- match.arg(scaleMethod)

  ncomp <- atlas_re$ncomp
  selectedFeat <- atlas_re$selectedFeat
  Bhat <- .coefSlice(atlas_re$reg_re$coefficients, ncomp)

  # If use rank transformed data, do rank transformation again after feature selection
  XX <- .selectCols(phenoAssay, selectedFeat)
  if(assayName == 'rank') XX <- .rankWithinCells(XX)
  XX <- .fit_matrix(XX)

  # Centre and scale XX implicitly: (XX - 1 Xmeans') diag(1/Xscals) %*% Bhat
  # equals XX %*% (Bhat/Xscals) - 1 Xmeans' (Bhat/Xscals), so a sparse XX
  # stays sparse.
  Xmeans <- rep(0, ncol(XX))
  Xscals <- rep(1, ncol(XX))
  if(scaleMethod == "byQuery"){
    if(atlas_re$center) Xmeans <- colMeans(XX)
    if(atlas_re$scale){
      if(atlas_re$center){
        Xscals <- .colSds(XX, Xmeans)
      } else {
        # As base::scale() without centring: root mean square
        Xscals <- sqrt(colSums(XX^2)/(nrow(XX) - 1))
      }
    }
  } else {
    if(!is.null(atlas_re$reg_re$Xmeans)) Xmeans <- atlas_re$reg_re$Xmeans
    if(atlas_re$scale) Xscals <- atlas_re$reg_re$Xscals
  }

  Bscal <- Bhat/Xscals
  offset <- as.numeric(crossprod(Xmeans, Bscal))
  if(!is.null(atlas_re$reg_re$Ymeans)) offset <- offset - atlas_re$reg_re$Ymeans

  Yhat <- as.matrix(XX %*% Bscal)
  Yhat <- sweep(Yhat, 2, offset)

  return(Yhat)
}
