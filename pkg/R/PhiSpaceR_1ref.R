#' PhiSpace using a single reference
#'
#' @param reference SCE. The reference.
#' @param query SCE or a list of SCE. The query.
#' @param phenotypes Charater. Which types of phenotypes to predict. If `NULL`, then have to specify `response`. Currently only support categorical phenotypes. For continuous ones, specify `response` directly.
#' @param response Named matrix. Rows correspond to cells (columns) in reference; columns correspond to phenotypes. If not `NULL`, then will override `phenotypes`. Can be continuous values such as age and BMI.
#' @param refAssay Character. Which assay in reference to use to train PhiSpace
#' @param queryAssay Character. Which assay in query to use for prediction; by default same as refAssay
#' @param regMethod Character. Regression method: one of PLS and PCA
#' @param ncomp Integer.
#' @param nfeat Integer.
#' @param selectedFeat Character.
#' @param center Logic.
#' @param scale Logic.
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.
#'
#' @return A list
#'
#' @export
PhiSpaceR_1ref <- function(
    reference,
    query,
    phenotypes = NULL,
    response = NULL,
    refAssay = "log1p",
    queryAssay = NULL,
    regMethod = c("PLS", "PCA"),
    ncomp = NULL,
    nfeat = NULL,
    selectedFeat = NULL,
    center = TRUE,
    scale = FALSE,
    DRinfo = FALSE
){

  if(!inherits(query, "list")) query <- list(query)

  # Check if refAssay is in reference
  if(is.null(queryAssay)) queryAssay <- refAssay
  if(!(refAssay %in% assayNames(reference))) stop("refAssay is not present in reference.")

  # Check if queryAssay is in all queries
  allAssayNames <- lapply(query, assayNames)
  allAssayNames <- Reduce(intersect, allAssayNames)
  if(!(queryAssay %in% allAssayNames)) stop("queryAssay needs to be present in every query.")


  regMethod <- match.arg(regMethod)

  ## Build Y matrix
  if(!is.null(response)){ # If response is provided

    YY <- as.matrix(response)
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

  ## Build atlas
  if(is.null(ncomp)) ncomp <- ncol(YY)

  # Define selectedFeat
  if(!is.null(selectedFeat)){ # if selectedFeat provided

    selectedFeat <- intersect(selectedFeat, rownames(reference))
    if(length(selectedFeat) == 0) stop("Reference and query don't share any selected genes.")

    impScores <- NULL

  } else {

    if(!is.null(nfeat)){ # if nfeat has been specified
      impScores <- mvr(
        t(assay(reference, refAssay)),
        YY,
        ncomp,
        method = regMethod,
        center = center, scale = scale
      )$coefficients[,,ncomp]
      selectedFeat <- selectFeat(impScores, nfeat)$selectedFeat
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
    assayName = refAssay,
    regMethod = regMethod,
    center = center,
    scale = scale,
    DRinfo = DRinfo
  )

  if(is.null(impScores)){

    impScores <- atlas_re$reg_re$coefficients[,,ncomp]
  }

  YrefHat <- phenotype(
    phenoAssay = t(assay(reference, refAssay)),
    atlas_re = atlas_re,
    assayName = refAssay
  )
  ## Project query
  PhiSpaceScore_l <- lapply(
    query,
    function(x){
      phenotype(
        phenoAssay = t(assay(x, queryAssay)),
        atlas_re = atlas_re,
        assayName = queryAssay
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
