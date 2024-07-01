#' PhiSpace continuous phenotyping
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
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. Disable to save memory.
#' @param storeUnNorm Store unnormalised raw PhiSpace scores or not. Default is `FALSE`.
#'
#' @return An SCE object with annotation results stored as a reducedDims slot.
#'
#' @export
PhiSpace <- function(
    reference,
    query,
    phenotypes = NULL,
    response = NULL,
    PhiSpaceAssay = "rank",
    regMethod = c("PLS", "PCA"),
    ncomp = NULL,
    nfeat = NULL,
    selectedFeat = NULL,
    center = TRUE,
    scale = FALSE,
    DRinfo = FALSE,
    storeUnNorm = FALSE
){

  PhiRes <- PhiSpaceR_1ref(
    reference = reference,
    query = query,
    phenotypes = phenotypes,
    response = response,
    PhiSpaceAssay = PhiSpaceAssay,
    regMethod = regMethod,
    ncomp = ncomp,
    nfeat = nfeat,
    selectedFeat = selectedFeat,
    center = center,
    scale = scale,
    DRinfo = DRinfo
  )

  if(storeUnNorm) reducedDim(query, "PhiSpaceNonNorm") <- PhiRes$PhiSpaceScore
  reducedDim(query, "PhiSpace") <- normPhiScores(PhiRes$PhiSpaceScore)

  return(
    list(
      annotatedQuery = query,
      OtherResults = PhiRes
    )
  )

}
