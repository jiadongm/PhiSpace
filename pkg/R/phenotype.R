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
  Bhat <- atlas_re$reg_re$coefficients[,,ncomp]

  # If use rank transformed data, do rank transformation again after feature selection
  if(assayName == 'rank'){
    XX <- RTassay(phenoAssay[ , selectedFeat])
  } else {
    XX <- phenoAssay[, selectedFeat]
  }


  if(scaleMethod == "byQuery"){
    XX_cent <- scale(XX, center = atlas_re$center, scale = atlas_re$scale)
  } else {
    if(atlas_re$scale){
      XX_cent <- scale(XX,
                       scale = atlas_re$reg_re$Xscals,
                       center = atlas_re$reg_re$Xmeans)
    } else {
      XX_cent <- scale(XX,
                       scale = F,
                       center = atlas_re$reg_re$Xmeans)

    }
  }

  if(is.null(atlas_re$reg_re$Ymeans)){
    toCent <- FALSE
  } else {

    toCent <- -atlas_re$reg_re$Ymeans
  }
  Yhat <- scale(
    XX_cent %*% Bhat,
    scale = F,
    center = toCent
  )

  return(Yhat)
}
