#' Soft phenotyping of query assay.
#'
#' @param assayToProj
#' @param atlas_re
#' @param assayName
#' @param scale
#'
#' @return
phenotype <- function(assayToProj, atlas_re, assayName = 'rank', scale = FALSE){

  ncomp <- atlas_re$ncomp
  selectFeat <- atlas_re$selectFeat
  Bhat <- atlas_re$reg_re$coefficients[,,ncomp]

  # If use rank transformed data, do rank transformation again after feature selection
  if(assayName == 'rank'){
    XX_cent <- RTassay(assayToProj[ , selectFeat])
  } else {
    XX_cent <- assayToProj[, selectFeat]
  }

  if(scale){
    XX_cent <- scale(XX_cent,
                     scale = atlas_re$reg_re$Xscals,
                     center = atlas_re$reg_re$Xmeans)
    Yhat <- scale(XX_cent %*% Bhat,
                  scale = atlas_re$reg_re$Xscals,
                  center = -atlas_re$reg_re$Ymeans)
  } else {
    XX_cent <- scale(XX_cent,
                     scale = F,
                     center = atlas_re$reg_re$Xmeans)
    Yhat <- scale(XX_cent %*% Bhat,
                  scale = F,
                  center = -atlas_re$reg_re$Ymeans)
  }

  Proj <- XX_cent %*% atlas_re$reg_re$loadings

  return(list(Yhat = Yhat, Proj = Proj))
}
