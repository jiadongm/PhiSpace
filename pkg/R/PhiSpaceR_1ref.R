#' PhiSpace using a single reference
#'
#' @param reference SingleCellExperiment object. The reference
#' @param query SingleCellExperiment object. The query
#' @param phenotypes Charater. Which types of phenotypes to predict
#' @param PhiSpaceAssay Character. Which assay to use as predictors
#' @param regMethod Character. Regression method: one of PLS and PCA
#' @param ncomp Integer.
#' @param nfeat Integer.
#' @param selectedFeat Character.
#' @param center Logic.
#' @param scale Logic.
#'
#' @return A list
#'
#' @export
PhiSpaceR_1ref <- function(reference,
                           query,
                           phenotypes,
                           PhiSpaceAssay = "rank",
                           regMethod = c("PLS", "PCA"),
                           ncomp = NULL,
                           nfeat = NULL,
                           selectedFeat = NULL,
                           center = TRUE,
                           scale = FALSE
){

  regMethod <- match.arg(regMethod)

  ## Build Y matrix
  YY <- codeY(reference, phenotypes)
  phenoDict <-
    data.frame(
      labs = colnames(YY),
      phenotypeCategory =
        rep(phenotypes,
            apply(as.data.frame(colData(reference)[,phenotypes]), 2,
                  function(x) length(unique(x)))
        )
    )

  ## Common genes and rank transform
  c(reference, query) %<-% KeepCommonGenes(reference, query)
  if(PhiSpaceAssay == "rank"){
    reference <- RankTransf(reference, PhiSpaceAssay)
    query <- RankTransf(query, PhiSpaceAssay)
  }

  ## Build atlas
  if(is.null(ncomp)) ncomp <- ncol(YY)

  # Define selectedFeat
  if(!is.null(selectedFeat)){ # if selectedFeat provided

    selectedFeat <- intersect(selectedFeat, rownames(reference))
    if(length(selectedFeat) == 0) stop("Reference and query don't share any selected genes.")

    impScores <- NULL

  } else {

    if(!is.null(nfeat)){ # if nfeat has been specified
      impScores <- mvr(t(assay(reference, PhiSpaceAssay)),
                       YY,
                       ncomp,
                       method = regMethod,
                       center = center, scale = scale)$coefficients[,,ncomp]
      selectedFeat <- selectFeat(impScores, nfeat)
    } else {

      impScores <- NULL
      selectedFeat <- rownames(reference)
    }

  }



  atlas_re <- SuperPC(reference = reference,
                      YY = YY,
                      ncomp = ncomp,
                      selectedFeat = selectedFeat,
                      assayName = PhiSpaceAssay,
                      regMethod = regMethod,
                      center = center,
                      scale = scale)
  YrefHat <- phenotype(phenoAssay = t(assay(reference, PhiSpaceAssay)),
                       atlas_re = atlas_re,
                       assayName = PhiSpaceAssay)
  ## Project query
  PhiSpaceScore <- phenotype(phenoAssay = t(assay(query, PhiSpaceAssay)),
                             atlas_re = atlas_re,
                             assayName = PhiSpaceAssay)


  return(
    list(
      ncomp = ncomp,
      impScores = impScores,
      phenoDict = phenoDict,
      selectedFeat = selectedFeat,
      YrefHat = YrefHat,
      PhiSpaceScore = PhiSpaceScore,
      center = center,
      scale = scale
    )
  )
}