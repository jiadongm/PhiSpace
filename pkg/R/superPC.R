#' Compute regression model using selected features.
#'
#' @param reference SingleCellExperiment object. The reference
#' @param YY Resonse matrix
#' @param ncomp Integer. Number of components
#' @param selectedFeat Character. Vector containing selected feature names
#' @param assayName Which assay in reference to use for prediction
#' @param scale Logic. Scale the predictor matrix or not
#' @param regMethod Character. Regression method to use, either PCA or PLS.
#' @param center Logic.
#' @param sparse Use sparse matrices or not.
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.
#'
#' @return A list containing regression input and output.
SuperPC <- function(
    reference,
    YY,
    ncomp,
    selectedFeat = NULL,
    assayName = 'logcounts',
    regMethod = c("PCA", "PLS"),
    center = TRUE,
    scale = FALSE,
    sparse = TRUE,
    DRinfo = FALSE
  )
{
  regMethod <- match.arg(regMethod)

  XX <- t(assay(reference, assayName))


  ## Prepare predictor matrix
  if(assayName == 'rank'){
    XX <- RTassay(XX[, selectedFeat])
  } else {
    XX <- XX[, selectedFeat]
  }


  reg_re <- mvr(
    XX,
    YY,
    ncomp = ncomp,
    method = regMethod,
    center = center,
    scale = scale,
    DRinfo = DRinfo
  )


  return(list(
    ncomp = ncomp,
    selectedFeat = selectedFeat,
    assayName = assayName,
    regMethod = regMethod,
    reg_re = reg_re,
    center = center,
    scale = scale
  ))
}
