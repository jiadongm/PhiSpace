#' Title
#'
#' @param reference `SingleCellExperiment` object.
#' @param assayName Character.
#' @param phenotypes Character.
#' @param YY Matrix.
#' @param ncomp Integer.
#' @param impScores Matrix.
#' @param center Logical.
#' @param scale Logical.
#' @param tune_ncomp Logic.
#' @param tune_nfeat Logic.
#' @param ncompLimits Vector.
#' @param ncompGridSize Integer.
#' @param ncompGrid Vector.
#' @param selectedFeat Vector.
#' @param nfeatLimits Vector.
#' @param nfeatGridSize Integer.
#' @param nfeatV Vector.
#' @param Kfolds Integer.
#' @param regMethod Character.
#' @param labelCode Character.
#' @param Ncores Integer.
#' @param seed Integer.
#'
#' @return A list
#' @export
tunePhiSpace <- function(reference,
                         assayName = "logcounts",
                         phenotypes = NULL,
                         YY = NULL,
                         ncomp = NULL,
                         impScores = NULL,
                         #
                         center = TRUE,
                         scale = FALSE,
                         #
                         tune_ncomp = TRUE,
                         tune_nfeat = TRUE,
                         #
                         ncompLimits = c(1,100),
                         ncompGridSize = 20,
                         ncompGrid = NULL,
                         selectedFeat = NULL,
                         #
                         nfeatLimits = c(10, 15000),
                         nfeatGridSize = 30,
                         nfeatV = NULL,
                         #
                         Kfolds = 5,
                         regMethod = 'PLS',
                         labelCode = "-1,1",
                         Ncores = 1,
                         seed = 5202056){

  if(tune_ncomp){

    tune_ncomp_res <- CVTune_ncomp(reference = reference,
                                   assayName = assayName,
                                   phenotypes = phenotypes,
                                   YY = YY,
                                   #
                                   center = center,
                                   scale = scale,
                                   #
                                   ncompLimits = ncompLimits,
                                   gridSize = ncompGridSize,
                                   ncompGrid = ncompGrid,
                                   selectedFeat = selectedFeat,
                                   #
                                   Kfolds = Kfolds,
                                   regMethod = regMethod,
                                   labelCode = labelCode,
                                   Ncores = Ncores,
                                   seed = seed)
    ncomp <- tune_ncomp_res$ncomp
    tuneNcompErrs <- tune_ncomp_res$tuneNcompErrs
  } else {

    tuneNcompErrs <- NULL

    if(is.null(phenotypes)){

      if(is.null(YY)) stop("ncomp and YY cannot be both NULL.")
    } else {

      YY <- codeY(reference, phenotypes)
    }

    if(is.null(ncomp)){

      ncomp <- ncol(YY)
    }
  }

  if(tune_nfeat){
    tune_nfeat_res <- CVTune_nfeat(reference = reference,
                                   assayName = assayName,
                                   phenotypes = phenotypes,
                                   YY = YY,
                                   ncomp = ncomp,
                                   impScores = impScores,
                                   #
                                   center = center,
                                   scale = scale,
                                   #
                                   nfeatLimits = nfeatLimits,
                                   gridSize = nfeatGridSize,
                                   nfeatV = nfeatV,
                                   #
                                   Kfolds = Kfolds,
                                   regMethod = regMethod,
                                   labelCode = labelCode,
                                   seed = seed,
                                   Ncores = Ncores)
    nfeat <- tune_nfeat_res$nfeat
    tuneNfeatErrs <- tune_nfeat_res$tuneNfeatErrs
    impScores <- tune_nfeat_res$impScores
    selectedFeat <- tune_nfeat_res$selectedFeat
  } else {
    tuneNfeatErrs <- NULL
    nfeat <- NULL
  }

  # If do not tune ncomp nor nfeat, simply calculate importance scores
  if(!tune_ncomp & !tune_nfeat){

    selectedFeat <- rownames(reference)

    impScores <- mvr(
      X = t(assay(reference, assayName)),
      Y = YY,
      ncomp = ncomp,
      method = regMethod,
      center = center,
      scale = scale
    )$coefficients[,,ncomp]
  }



  return(list(
    ncomp = ncomp,
    nfeat = nfeat,
    impScores = impScores,
    selectedFeat = selectedFeat,
    tuneNfeatErrs = tuneNfeatErrs,
    tuneNcompErrs = tuneNcompErrs
  ))
}
