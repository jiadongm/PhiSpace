#' PhiSpace annotation using a single reference
#'
#' Core single-reference implementation for PhiSpace continuous phenotyping. Given an annotated reference
#' dataset and one or more query datasets, this function builds a PLS or PCA regression model from the
#' reference and projects the query onto the learned phenotype space. Both raw and normalised PhiSpace
#' scores are returned. For a higher-level wrapper that handles multiple references, see [PhiSpace()].
#'
#' @param reference A `SingleCellExperiment` (SCE) object. The annotated reference dataset. Must contain
#'   an assay named by `refAssay`.
#' @param query An SCE object or a list of SCE objects. The query dataset(s) to annotate.
#'   Each must contain an assay named by `queryAssay`.
#' @param phenotypes Character vector. Column name(s) in `colData(reference)` specifying which
#'   phenotypes (e.g. cell type, disease condition) to predict. If `NULL`, `response` must be provided.
#'   Currently only supports categorical phenotypes; for continuous phenotypes, specify `response` directly.
#' @param response Named matrix. Rows correspond to cells (columns) in `reference`; columns correspond
#'   to phenotypes. If not `NULL`, overrides `phenotypes`. Can contain continuous values such as age or BMI.
#' @param refAssay Character. Name of the assay in `reference` to use for model training (default: `"log1p"`).
#' @param queryAssay Character. Name of the assay in `query` to use for prediction. If `NULL` (default),
#'   uses the same value as `refAssay`.
#' @param regMethod Character. Regression method: one of `"PLS"` (partial least squares, default) or
#'   `"PCA"` (principal component regression).
#' @param ncomp Integer. Number of PLS/PCA components. If `NULL` (default), set to the total number
#'   of phenotype levels.
#' @param nfeat Integer. Number of top-ranked features to select per phenotype for model fitting. If `NULL`
#'   (default), all features are used. See Details.
#' @param selectedFeat Character vector. A set of pre-selected feature names to use for model fitting. If
#'   provided, overrides `nfeat`.
#' @param center Logical. Whether to center the feature matrix before regression (default: `TRUE`). See Details.
#' @param scale Logical. Whether to scale the feature matrix before regression (default: `FALSE`). See Details.
#' @param DRinfo Logical. Whether to return dimension reduction information from PCA or PLS.
#'   Disable to save memory (default: `FALSE`).
#' @param cellTypeThreshold Integer or NULL. If a positive integer, cell types with fewer than this many
#'   cells in the reference will be removed before model fitting. Only used when `phenotypes` is provided.
#'   Default is `NULL` (no filtering).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{ncomp}{Integer. Number of components used.}
#'   \item{impScores}{Matrix. Feature importance scores (regression coefficients at the final component).}
#'   \item{phenoDict}{Data frame mapping phenotype labels to their categories, or `NULL` if `response`
#'     was provided directly.}
#'   \item{selectedFeat}{Character vector. Features used in the final model.}
#'   \item{YrefHat}{Matrix. Raw (unnormalised) predicted scores for the reference cells.}
#'   \item{YrefHatNorm}{Matrix. Normalised predicted scores for the reference cells (via [normPhiScores()]).}
#'   \item{PhiSpaceScore}{Matrix or list of matrices. Raw (unnormalised) PhiSpace scores for the query
#'     dataset(s). A single matrix if one query was provided; a list of matrices if multiple queries.}
#'   \item{PhiSpaceNorm}{Matrix or list of matrices. Normalised PhiSpace scores for the query dataset(s)
#'     (via [normPhiScores()]). Same structure as `PhiSpaceScore`.}
#'   \item{center}{Logical. Whether centering was applied.}
#'   \item{scale}{Logical. Whether scaling was applied.}
#'   \item{atlas_re}{List. Internal model object from `SuperPC()`, containing regression results and
#'     preprocessing parameters needed for prediction.}
#' }
#'
#' @details
#' - **Parameter tuning** By default, PhiSpace takes `ncomp` equal to the total number of phenotypes.
#' For example, if there are 20 cell types and 5 sample sources (e.g. *in vivo* and *in vitro*) defined
#' in the reference, `ncomp` will be set to be 25. By default, PhiSpace doesn't do any feature selection.
#' However, feature selection is very straightforward in PhiSpace. The partial least squares model
#' automatically rank all features according to their contribution to predicting each phenotype. For
#' example, if there are 25 phenotypes then PLS will return 25 rankings of all features. Based on this,
#' customised feature selection can be done in two ways. The user either specify `nfeat`, i.e. the number
#' of features (or markers) used to predict each phenotype; or specify `selectedFeat`, i.e. a subset of
#' features used to predict all phenotypes (e.g. highly variable genes). If `nfeat` is provided, then
#' `selectedFeat` will be automatically defined as the union of the `nfeat` markers of all phenotypes.
#'
#' - **Center and scale** By default, PhiSpace only center the data but do not scale. This is because we
#' often use more features than standard single-cell pipelines (e.g. Seurat selects 2,000 highly variable
#' genes) for a better prediction. This means that we might have included genes with small variance and
#' scaling would greatly inflate the expression levels of these genes.
#'
#' @references
#' Mao J., Deng, Y. and Lê Cao, K.-A. (2024). Phi-Space: Continuous phenotyping of single-cell multi-omics data. bioRxiv.
#'
#' @seealso [PhiSpace()] for the main user-facing wrapper that handles multiple references,
#'   [normPhiScores()] for score normalisation, [tunePhiSpace()] for parameter tuning.
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
    DRinfo = FALSE,
    cellTypeThreshold = NULL
){

  if(!inherits(query, "list")) query <- list(query)

  # Validate cellTypeThreshold
  if(!is.null(cellTypeThreshold)){
    if(!(length(cellTypeThreshold) == 1 && is.numeric(cellTypeThreshold) &&
         cellTypeThreshold == as.integer(cellTypeThreshold) && cellTypeThreshold > 0)){
      stop("cellTypeThreshold must be a positive integer or NULL.")
    }
  }

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
    if(sum(is.na(colData(reference)[,phenotypes]))) stop("phenotypes cannot contain NAs.")

    # Filter rare cell types if cellTypeThreshold is set
    if(!is.null(cellTypeThreshold)){
      for(ph in phenotypes){
        cellTypeCounts <- table(colData(reference)[, ph])
        rareCellTypes <- names(cellTypeCounts[cellTypeCounts < cellTypeThreshold])
        if(length(rareCellTypes) > 0){
          message(
            "Removing cell types with fewer than ", cellTypeThreshold,
            " cells in '", ph, "': ",
            paste(rareCellTypes, " (n=", cellTypeCounts[rareCellTypes], ")",
                  sep = "", collapse = ", ")
          )
          keepCells <- !(as.character(colData(reference)[, ph]) %in% rareCellTypes)
          reference <- reference[, keepCells]
        }
      }
      if(ncol(reference) == 0) stop("No cells remain after filtering rare cell types.")
    }

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

  ## Normalize scores
  YrefHatNorm <- normPhiScores(YrefHat)
  if(is.list(PhiSpaceScore_l)){
    PhiSpaceNorm_l <- lapply(PhiSpaceScore_l, normPhiScores)
  } else {
    PhiSpaceNorm_l <- normPhiScores(PhiSpaceScore_l)
  }

  return(
    list(
      ncomp = ncomp,
      impScores = impScores,
      phenoDict = phenoDict,
      selectedFeat = selectedFeat,
      YrefHat = YrefHat,
      YrefHatNorm = YrefHatNorm,
      PhiSpaceScore = PhiSpaceScore_l,
      PhiSpaceNorm = PhiSpaceNorm_l,
      center = center,
      scale = scale,
      atlas_re = atlas_re
    )
  )
}
