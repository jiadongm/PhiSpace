#' Soft phenotyping of query assay.
#'
#' @param phenoAssay Matrix. Query assay to phenotype.
#' @param atlas_re List.
#' @param assayName Character.
#' @param scale Logic.
#' @param scaleMethod Character.
#'
#' @return A list containing
#' \item{Yhat}{}
#' \item{Proj}{}
phenotype <- function(phenoAssay,
                      atlas_re,
                      assayName = 'rank',
                      scaleMethod = c("byQuery", "byRef"),
                      scale = FALSE){

  scaleMethod <- match.arg(scaleMethod)

  ncomp <- atlas_re$ncomp
  selectedFeat <- atlas_re$selectedFeat
  Bhat <- atlas_re$reg_re$coefficients[,,ncomp]

  # If use rank transformed data, do rank transformation again after feature selection
  if(assayName == 'rank'){
    XX <- RTassay(phenoAssay[ , selectedFeat])
  } else {
    XX <- phenoAssay[, selectedFeat]
  }


  if(scaleMethod == "byQuery"){
    XX <- scale(XX, center = TRUE, scale = scale)
  } else {
    if(scale){
      XX <- scale(XX,
                  scale = atlas_re$reg_re$Xscals,
                  center = atlas_re$reg_re$Xmeans)
    } else {
      XX_cent <- scale(XX_cent,
                       scale = F,
                       center = atlas_re$reg_re$Xmeans)

    }
  }

  Yhat <- scale(XX %*% Bhat,
                scale = F,
                center = -atlas_re$reg_re$Ymeans)

  return(Yhat)
}
