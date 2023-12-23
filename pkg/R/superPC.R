#' Compute regression model using selected features.
#'
#' @param reference SingleCellExperiment object. The reference
#' @param YY Resonse matrix
#' @param ncomp Integer. Number of components
#' @param selectedFeat Character. Vector containing selected feature names
#' @param assayName Which assay in reference to use for prediction
#' @param scale Logic. Scale the predictor matrix or not
#' @param regMethod Character. Regression method to use, either PCA or PLS.
#'
#' @return A list containing regression input and output.
SuperPC <- function(reference,
                    YY,
                    ncomp,
                    selectedFeat = NULL,
                    assayName = 'logcounts',
                    regMethod = c("PCA", "PLS"),
                    scale = FALSE)
{
  regMethod <- match.arg(regMethod)

  XX <- as.matrix(t(assay(reference, assayName)))


  ## Prepare predictor matrix
  if(assayName == 'rank'){
    XX <- RTassay(XX[, selectedFeat])
  } else {
    XX <- XX[, selectedFeat]
  }


  reg_re <- mvr(XX, YY,
                ncomp = ncomp,
                method = regMethod,
                scale = scale)


  return(list(
    ncomp = ncomp,
    selectedFeat = selectedFeat,
    assayName = assayName,
    regMethod = regMethod,
    reg_re = reg_re,
    scale = scale
  ))
}
