#' Train a PhiSpace model from a reference
#'
#' Fits the PhiSpace regression model on an annotated reference and returns it
#' as a `PhiSpaceModel` object. The model can be saved (for example with
#' [saveRDS()]) and used later to score queries with [predict()][predict.PhiSpaceModel]
#' or [PhiSpace()], without the reference cells.
#'
#' [PhiSpaceR_1ref()] and [PhiSpace()] fit the same model on the genes that the
#' reference shares with the query. `trainPhiSpace()` uses all reference genes,
#' or the genes given in `genes`. To score a query that lacks some genes, train
#' with `genes` set to the genes of that query (for example a targeted spatial
#' panel); the query must contain every gene in the model.
#'
#' @param reference A `SingleCellExperiment` (or `SummarizedExperiment`)
#'   object. The annotated reference dataset.
#' @param phenotypes Character vector. Column name(s) in `colData(reference)`
#'   with the phenotypes to predict. If `NULL`, `response` must be provided.
#' @param response Named matrix. Rows correspond to cells (columns) in
#'   `reference`; columns correspond to phenotypes. If not `NULL`, overrides
#'   `phenotypes`.
#' @param refAssay Character. Name of the assay in `reference` used for
#'   training. Queries must be preprocessed in the same way, so there is no
#'   default.
#' @param regMethod Character. `"PLS"` (default) or `"PCA"`.
#' @param ncomp Integer. Number of components. If `NULL`, the total number of
#'   phenotype levels.
#' @param nfeat Integer. Number of top features per phenotype. If `NULL`, all
#'   genes are used. See [PhiSpaceR_1ref()].
#' @param selectedFeat Character vector. Pre-selected features; overrides
#'   `nfeat`.
#' @param center,scale Logical. Centring and scaling of the features.
#' @param DRinfo Logical. Whether to keep the PLS or PCA scores and loadings.
#' @param cellTypeThreshold Integer or `NULL`. Cell types with fewer cells than
#'   this are removed before training. Only used with `phenotypes`.
#' @param genes Character vector or `NULL`. If given, training uses only these
#'   genes (those present in the reference).
#' @param referenceName,species,geneIdType Character or `NULL`. Optional
#'   descriptions of the reference, stored in the model.
#'
#' @return A `PhiSpaceModel` object: a list with
#' \describe{
#'   \item{atlas_re}{The fitted model, as returned in `PhiSpaceR_1ref()$atlas_re`.}
#'   \item{selectedFeat}{Genes used by the model, in model order.}
#'   \item{impScores}{Feature importance scores (see [PhiSpaceR_1ref()]).}
#'   \item{phenoDict}{Phenotype labels and their categories, or `NULL` if
#'     `response` was given.}
#'   \item{responseNames}{Names of the score columns.}
#'   \item{refAssay}{The training assay.}
#'   \item{nCells}{Number of reference cells used for training.}
#'   \item{referenceName, species, geneIdType}{As supplied, or `NULL`.}
#'   \item{format_version, phispace_version, created}{Model format version,
#'     the PhiSpace version that trained the model, and the date.}
#' }
#'
#' @seealso [predict.PhiSpaceModel()], [PhiSpace()], [PhiSpaceR_1ref()].
#'
#' @export
trainPhiSpace <- function(
    reference,
    phenotypes = NULL,
    response = NULL,
    refAssay,
    regMethod = c("PLS", "PCA"),
    ncomp = NULL,
    nfeat = NULL,
    selectedFeat = NULL,
    center = TRUE,
    scale = FALSE,
    DRinfo = FALSE,
    cellTypeThreshold = NULL,
    genes = NULL,
    referenceName = NULL,
    species = NULL,
    geneIdType = NULL
){

  regMethod <- match.arg(regMethod)
  .validate_cellTypeThreshold(cellTypeThreshold)
  if(!(refAssay %in% assayNames(reference))) stop("refAssay is not present in reference.")

  if(!is.null(response)){
    YY <- as.matrix(response)
    phenoDict <- NULL
  } else {
    if(is.null(phenotypes)) stop("phenotypes and response cannot both be NULL.")
    reference <- .filterRareTypes(reference, phenotypes, cellTypeThreshold)
    coded <- .codePhenotypes(reference, phenotypes)
    YY <- coded$YY
    phenoDict <- coded$phenoDict
  }

  refX <- assay(reference, refAssay)
  featNames <- rownames(refX)
  if(!is.null(genes)){
    featNames <- intersect(featNames, genes)
    if(length(featNames) == 0) stop("None of genes is present in reference.")
  }

  .trainPhiSpace(
    refX = refX, YY = YY, phenoDict = phenoDict, featNames = featNames,
    refAssay = refAssay, regMethod = regMethod, ncomp = ncomp, nfeat = nfeat,
    selectedFeat = selectedFeat, center = center, scale = scale,
    DRinfo = DRinfo, referenceName = referenceName, species = species,
    geneIdType = geneIdType
  )$model
}


#' Score a query with a stored PhiSpace model
#'
#' @param object A `PhiSpaceModel` from [trainPhiSpace()].
#' @param newdata The query: a `SummarizedExperiment` (for example a
#'   `SingleCellExperiment`) or a gene by cell matrix. It must contain every
#'   gene in `object$selectedFeat`.
#' @param assay Character. Assay of `newdata` to use. Defaults to the model's
#'   `refAssay`. With `"rank"`, the genes are re-ranked within each cell over
#'   the model genes, as in [PhiSpaceR_1ref()].
#' @param ... Not used.
#'
#' @return A cell by phenotype matrix of raw PhiSpace scores. Use
#'   [normPhiScores()] to normalise them.
#'
#' @export
predict.PhiSpaceModel <- function(object, newdata, assay = NULL, ...){

  .checkPhiSpaceModel(object)
  if(is.null(assay)) assay <- object$refAssay

  if(methods::is(newdata, "SummarizedExperiment")){
    if(!(assay %in% assayNames(newdata))) stop("assay '", assay, "' is not present in newdata.")
    X <- SummarizedExperiment::assay(newdata, assay)
  } else {
    X <- newdata
  }
  if(is.null(rownames(X))) stop("newdata must have gene names as row names.")

  .predictPhiSpace(object, X, assay)
}


#' @exportS3Method print PhiSpaceModel
print.PhiSpaceModel <- function(x, ...){

  ar <- x$atlas_re
  cat("PhiSpaceModel",
      if(!is.null(x$referenceName)) paste0("(", x$referenceName, ")"), "\n")
  cat(" ", length(x$responseNames), "phenotypes,", length(x$selectedFeat),
      "genes,", x$nCells, "reference cells\n")
  cat("  Assay:", x$refAssay, " Method:", ar$regMethod, " ncomp:", ar$ncomp,
      " center:", ar$center, " scale:", ar$scale, "\n")
  if(!is.null(x$species) || !is.null(x$geneIdType)){
    cat("  Species:", if(is.null(x$species)) "-" else x$species,
        " Gene IDs:", if(is.null(x$geneIdType)) "-" else x$geneIdType, "\n")
  }
  cat("  Trained with PhiSpace", x$phispace_version, "on", x$created, "\n")
  invisible(x)
}


## Model format written by this version of PhiSpace. Increase it when the
## structure of a PhiSpaceModel changes, and keep reading older formats.
.PhiSpaceModel_format <- 1L


## Fit a model on the genes featNames of the gene by cell matrix refX.
## Returns list(model, YrefHat); YrefHat (raw reference scores) only when
## scoreReference = TRUE, while the selected reference matrix is in memory.
.trainPhiSpace <- function(refX, YY, phenoDict, featNames, refAssay, regMethod,
                           ncomp, nfeat, selectedFeat, center, scale, DRinfo,
                           referenceName = NULL, species = NULL,
                           geneIdType = NULL, scoreReference = FALSE){

  if(is.null(ncomp)) ncomp <- ncol(YY)

  if(!is.null(selectedFeat)){ # if selectedFeat provided
    selectedFeat <- intersect(selectedFeat, featNames)
    if(length(selectedFeat) == 0) stop("Reference and query don't share any selected genes.")
    impScores <- NULL
  } else {
    if(!is.null(nfeat)){ # if nfeat has been specified
      impScores <- .coefSlice(
        mvr(
          .cellsByGenes(refX, featNames),
          YY,
          ncomp,
          method = regMethod,
          center = center, scale = scale,
          keepComps = ncomp
        )$coefficients,
        ncomp
      )
      selectedFeat <- selectFeat(impScores, nfeat)$selectedFeat
    } else {
      impScores <- NULL
      selectedFeat <- featNames
    }
  }

  refXsel <- .cellsByGenes(refX, selectedFeat)
  atlas_re <- SuperPC(
    reference = refXsel,
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
    impScores <- .coefSlice(atlas_re$reg_re$coefficients, ncomp)
  }
  YrefHat <- NULL
  if(scoreReference){
    YrefHat <- phenotype(
      phenoAssay = refXsel,
      atlas_re = atlas_re,
      assayName = refAssay
    )
  }

  model <- structure(
    list(
      atlas_re = atlas_re,
      selectedFeat = selectedFeat,
      impScores = impScores,
      phenoDict = phenoDict,
      responseNames = colnames(YY),
      refAssay = refAssay,
      nCells = nrow(YY),
      referenceName = referenceName,
      species = species,
      geneIdType = geneIdType,
      format_version = .PhiSpaceModel_format,
      phispace_version = as.character(utils::packageVersion("PhiSpace")),
      created = format(Sys.Date())
    ),
    class = "PhiSpaceModel"
  )

  list(model = model, YrefHat = YrefHat)
}


## Raw scores of the gene by cell matrix X. assayName decides re-ranking,
## as the query assay name does in PhiSpaceR_1ref().
.predictPhiSpace <- function(model, X, assayName){

  genes <- model$selectedFeat
  missing <- setdiff(genes, rownames(X))
  if(length(missing) > 0){
    stop("The query lacks ", length(missing), " of the ", length(genes),
         " model genes, for example: ",
         paste(utils::head(missing, 5), collapse = ", "),
         ". Train the model with genes = the query genes.")
  }

  phenotype(
    phenoAssay = .cellsByGenes(X, genes),
    atlas_re = model$atlas_re,
    assayName = assayName
  )
}


.checkPhiSpaceModel <- function(model){

  if(!inherits(model, "PhiSpaceModel")) stop("object is not a PhiSpaceModel.")
  if(is.null(model$format_version) || model$format_version > .PhiSpaceModel_format){
    stop("This PhiSpaceModel has format version ", model$format_version,
         "; this PhiSpace version reads format ", .PhiSpaceModel_format,
         " or older. Update PhiSpace.")
  }
  invisible(model)
}


.validate_cellTypeThreshold <- function(cellTypeThreshold){

  if(!is.null(cellTypeThreshold)){
    if(!(length(cellTypeThreshold) == 1 && is.numeric(cellTypeThreshold) &&
         cellTypeThreshold == as.integer(cellTypeThreshold) && cellTypeThreshold > 0)){
      stop("cellTypeThreshold must be a positive integer or NULL.")
    }
  }
}


## Check the phenotype columns for NAs and remove cell types with fewer than
## cellTypeThreshold cells.
.filterRareTypes <- function(reference, phenotypes, cellTypeThreshold){

  if(sum(is.na(colData(reference)[,phenotypes]))) stop("phenotypes cannot contain NAs.")

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

  reference
}


## Coded response matrix and phenotype dictionary.
.codePhenotypes <- function(reference, phenotypes){

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
  list(YY = YY, phenoDict = phenoDict)
}
