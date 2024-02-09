#' PhiSpace using a single reference
#'
#' @param reference SCE. The reference.
#' @param query SCE or a list of SCE. The query.
#' @param phenotypes Charater. Which types of phenotypes to predict. If `NULL`, then have to specify `response`.
#' @param response Named matrix. Rows correpond to cells (columns) in reference; columns correspond to phenotypes. If not `NULL`, then will override `phenotypes`.
#' @param PhiSpaceAssay Character. Which assay to use to train
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
                           phenotypes = NULL,
                           response = NULL,
                           PhiSpaceAssay = "rank",
                           regMethod = c("PLS", "PCA"),
                           ncomp = NULL,
                           nfeat = NULL,
                           selectedFeat = NULL,
                           center = TRUE,
                           scale = FALSE
){

  if(!inherits(query, "list")) query <- list(query)

  if(!(PhiSpaceAssay %in% assayNames(reference))) stop("PhiSpaceAssay is not present in reference.")

  # Intersection of names of assays in all queries
  allAssayNames <- lapply(query, assayNames)
  allAssayNames <- Reduce(intersect, allAssayNames)
  if(!(PhiSpaceAssay %in% allAssayNames)) stop("PhiSpaceAssay needs to be present in every query.")

  regMethod <- match.arg(regMethod)

  ## Build Y matrix
  if(!is.null(response)){ # If response is provided

    YY <- response
    phenoDict <- NULL

  } else {

    # Phenotypes has to be specified
    if(is.null(phenotypes)) stop("phenotypes and response cannot both be NULL.")

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
  }

  ## Common genes and rank transform
  featNames <- lapply(query, rownames)
  featNames <- Reduce(intersect, featNames)
  featNames <- intersect(rownames(reference), featNames)
  reference <- reference[featNames,]
  query <- lapply(query, function(x) x[featNames, ])

  if(PhiSpaceAssay == "rank"){
    reference <- RankTransf(reference, PhiSpaceAssay)
    query <- lapply(query, RankTransf, assayname = PhiSpaceAssay)
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


  atlas_re <- SuperPC(
    reference = reference,
    YY = YY,
    ncomp = ncomp,
    selectedFeat = selectedFeat,
    assayName = PhiSpaceAssay,
    regMethod = regMethod,
    center = center,
    scale = scale
  )
  YrefHat <- phenotype(
    phenoAssay = t(assay(reference, PhiSpaceAssay)),
    atlas_re = atlas_re,
    assayName = PhiSpaceAssay
  )
  ## Project query
  PhiSpaceScore_l <- lapply(
    query,
    function(x){
      phenotype(
        phenoAssay = t(assay(x, PhiSpaceAssay)),
        atlas_re = atlas_re,
        assayName = PhiSpaceAssay
      )
    }
  )
  if(length(PhiSpaceScore_l) == 1) PhiSpaceScore_l <- PhiSpaceScore_l[[1]]


  return(
    list(
      ncomp = ncomp,
      impScores = impScores,
      phenoDict = phenoDict,
      selectedFeat = selectedFeat,
      YrefHat = YrefHat,
      PhiSpaceScore = PhiSpaceScore_l,
      center = center,
      scale = scale,
      atlas_re = atlas_re
    )
  )
}
