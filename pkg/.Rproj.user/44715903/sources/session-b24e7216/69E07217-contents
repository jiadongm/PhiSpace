#' Compute regression model.
#'
#' @param re
#' @param reference
#' @param YY
#' @param labelName
#' @param ncomp
#' @param Nselect
#' @param impScores
#' @param selectFeat
#' @param assayName
#' @param regMethod
#' @param labelCode
#'
#' @return
SuperPC <- function(re = NULL,
                    reference = NULL,
                    YY = NULL,
                    labelName = NULL,
                    ncomp = NULL,
                    Nselect = NULL,
                    impScores = NULL,
                    selectFeat = NULL,
                    # scale = F,
                    assayName = 'logcounts',
                    regMethod = "PCA",
                    labelCode = "-1,1"){
  # re is a list containing parameter tuning results
  # here we allow user to use differnet ncomp, Nselect from those in re

  if(is.null(re)){

    if(is.null(ncomp)) stop("Have to provide ncomp.")
    if(is.null(reference)) stop("Have to provide reference object.")
    if(is.null(YY)){
      if(is.null(labelName)) stop("Have to provide either labelName or response matrix YY.")
    }
    if(is.null(selectFeat)){
      if(is.null(Nselect) | is.null(impScores)) stop("Have to provide either selectFeat or both impScores & Nselect.")
    }

    re <- list()
    re$ncomp <- ncomp
    XX <- as.matrix(t(assay(reference, assayName)))
    re$XX <- XX
    if(is.null(YY)){
      Ytrain <- colData(reference)[,labelName]
      classLabels <- names(table(Ytrain))
      YY <- sapply(1:length(Ytrain),
                   function(x){

                     if(labelCode == "-1/+1"){

                       out <- as.numeric(classLabels == Ytrain[x])
                       out[out == 0] <- -1

                     } else {

                       out <- as.numeric(classLabels == Ytrain[x])

                     }

                     return(out)
                   })
      YY <- t(YY)
      dimnames(YY) <- list(rownames(XX), classLabels)
      re$YY <- YY
    } else {
      re$YY <- YY
    }
    if(is.null(selectFeat)){
      re$Nselect <- Nselect
      re$impScores <- impScores

      # Select Nselect top features from each column (label) of impScores
      orderByCol <-
        apply(impScores, 2,
              function(x){
                names(x) <- rownames(impScores)
                names(sort(abs(x), decreasing = T))
              })
      re$selectFeat <- unique(
        as.vector(
          orderByCol[1:Nselect, ]
        )
      )

    } else {
      re$selectFeat <- selectFeat
    }
    re$regMethod <- regMethod
  } else {

    if(!is.null(ncomp)) re$ncomp <- ncomp
    if(!is.null(Nselect)) re$Nselect <- Nselect

  }


  ## Prepare predictor matrix
  if(assayName == 'rank'){
    # Rank transform
    XX_RT <- RTassay(re$XX[, re$selectFeat])
  } else {
    XX_RT <- re$XX[, re$selectFeat]
  }


  reg_re <- mvr(XX_RT, re$YY,
                ncomp = re$ncomp,
                method = re$regMethod)
  re$reg_re <- reg_re
  return(re)
}
