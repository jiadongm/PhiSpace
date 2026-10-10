#' Train a PhiSpace model from a reference
#'
#' Fits the PhiSpace regression model on an annotated reference and returns it
#' as a `PhiSpaceModel` object. The model can be saved (for example with
#' [saveRDS()]) and used later to score queries with [predict()][predict.PhiSpaceModel]
#' or [PhiSpace()], without the reference cells.
#'
#' [PhiSpaceR_1ref()] and [PhiSpace()] fit the same model on the genes that the
#' reference shares with the query. `trainPhiSpace()` uses all reference genes,
#' or the genes given in `genes`.
#'
#' A query may lack some of the model genes, for example a targeted spatial
#' panel. With `keepStats = TRUE` (the default), the model stores the
#' cross-products \eqn{X'X} and \eqn{X'Y} of the centred (and scaled)
#' reference over the model genes. [predict()][predict.PhiSpaceModel] then
#' refits the model on the genes that the query has, without the reference
#' cells. The refitted model equals a model trained on the reference with
#' `selectedFeat` set to those genes. Genes are not selected again (`nfeat`):
#' the refit uses the model genes that the query has.
#'
#' The statistics are not kept for the `"rank"` assay. Ranks are computed
#' within each cell over the model genes, so removing genes changes the
#' reference values and an exact refit needs the reference cells. A `"rank"`
#' model therefore scores a query that lacks genes approximately, with a
#' warning (see [predict.PhiSpaceModel()]).
#'
#' \eqn{X'X} has one row and one column per model gene: it takes
#' \eqn{8 G^2} bytes for \eqn{G} genes, for example 72 MB for 3,000 genes and
#' 3.2 GB for 20,000 genes. Use `nfeat`, `selectedFeat` or `genes` to train on
#' fewer genes, or `keepStats = FALSE`.
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
#' @param keepStats Logical. Whether to store the cross-product statistics
#'   that refit the model on a subset of its genes. Ignored for the `"rank"`
#'   assay. See Details.
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
#'   \item{stats}{`NULL`, or a list with `XtX` (genes x genes) and `XtY`
#'     (genes x phenotypes): the cross-products of the centred (and scaled)
#'     reference over the model genes, with the uncoded response.}
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
    keepStats = TRUE,
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
    DRinfo = DRinfo, keepStats = keepStats, referenceName = referenceName,
    species = species, geneIdType = geneIdType
  )$model
}


#' Score a query with a stored PhiSpace model
#'
#' @param object A `PhiSpaceModel` from [trainPhiSpace()].
#' @param newdata The query: a `SummarizedExperiment` (for example a
#'   `SingleCellExperiment`) or a gene by cell matrix, with gene names as row
#'   names. It must share at least one gene with the model.
#' @param assay Character. Assay of `newdata` to use. Defaults to the model's
#'   `refAssay`. With `"rank"`, the genes are re-ranked within each cell over
#'   the model genes, as in [PhiSpaceR_1ref()].
#' @param ... Not used.
#'
#' @details
#' If the query lacks some of the model genes, the model is first restricted
#' to the genes that the query has:
#' * A model with stored statistics (see [trainPhiSpace()]) is refitted
#'   exactly on these genes, with a message. The scores equal those of a model
#'   trained on the reference with `selectedFeat` set to these genes.
#' * A model without statistics (a `"rank"` model, a model trained with
#'   `keepStats = FALSE`, or `PhiSpaceR_1ref()$model`) keeps the coefficients
#'   of these genes and drops the others, with a warning. These scores are
#'   approximate: the error grows with the share of the coefficients that the
#'   missing genes hold, which the warning reports.
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
  cat("  Refits on fewer genes:",
      if(is.null(x$stats)) "approximate (no stored statistics)" else "exact (statistics stored)",
      "\n")
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
## keepStats stores the cross-products for exact refits (not for "rank").
.trainPhiSpace <- function(refX, YY, phenoDict, featNames, refAssay, regMethod,
                           ncomp, nfeat, selectedFeat, center, scale, DRinfo,
                           keepStats = FALSE, referenceName = NULL,
                           species = NULL, geneIdType = NULL,
                           scoreReference = FALSE){

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
  stats <- NULL
  if(keepStats && refAssay != "rank"){
    statBytes <- 8 * length(selectedFeat)^2
    if(statBytes > 1e9){
      message("Storing X'X for ", length(selectedFeat), " genes (",
              format(statBytes/1e9, digits = 2), " GB). Use nfeat, ",
              "selectedFeat or genes for fewer genes, or keepStats = FALSE.")
    }
    stats <- .crossStats(refXsel, YY, atlas_re$reg_re)
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
      stats = stats,
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
## as the query assay name does in PhiSpaceR_1ref(). A model with genes that
## X lacks is first restricted to the genes of X.
.predictPhiSpace <- function(model, X, assayName){

  model <- .modelForGenes(model, rownames(X))
  genes <- model$selectedFeat

  phenotype(
    phenoAssay = .cellsByGenes(X, genes),
    atlas_re = model$atlas_re,
    assayName = assayName
  )
}


## The model restricted to the model genes present in queryGenes, with a
## message (exact refit) or a warning (coefficients of the shared genes only).
.modelForGenes <- function(model, queryGenes){

  genes <- model$selectedFeat
  shared <- genes[genes %in% queryGenes]
  nMissing <- length(genes) - length(shared)
  if(nMissing == 0) return(model)
  if(length(shared) == 0) stop("The query shares none of the ", length(genes), " model genes.")

  lacks <- paste0("The query lacks ", nMissing, " of the ", length(genes),
                  " model genes (for example ",
                  paste(utils::head(setdiff(genes, shared), 3), collapse = ", "),
                  ")")
  if(!is.null(model$stats)){
    message(lacks, "; refitting the model on the ", length(shared), " shared genes.")
  } else {
    B <- .coefSlice(model$atlas_re$reg_re$coefficients, model$atlas_re$ncomp)
    lost <- 1 - sum(B[shared, ]^2)/sum(B^2)
    reason <- if(model$refAssay == "rank"){
      "a \"rank\" model cannot be refitted without the reference cells"
    } else {
      "the model has no stored statistics (see keepStats in trainPhiSpace())"
    }
    warning(lacks, ", and ", reason, ". The scores use the coefficients of the ",
            length(shared), " shared genes and are approximate; the missing ",
            "genes hold ", format(100 * lost, digits = 2),
            "% of the sum of squared coefficients.", call. = FALSE)
  }

  .subsetPhiSpaceModel(model, shared)
}


## The model on genes, a subset of model$selectedFeat in model order. With
## stored statistics the model is refitted on these genes, as if trained with
## selectedFeat = genes. Without them, the coefficients of these genes are
## kept unchanged.
.subsetPhiSpaceModel <- function(model, genes){

  ar <- model$atlas_re
  reg <- ar$reg_re
  ncomp <- ar$ncomp
  idx <- match(genes, model$selectedFeat)

  if(!is.null(model$stats)){
    if(length(genes) < ncomp){
      stop("Only ", length(genes), " model genes are shared with the query, ",
           "fewer than the ", ncomp, " components of the model.")
    }
    S <- model$stats$XtX[idx, idx, drop = FALSE]
    XtY <- model$stats$XtY[idx, , drop = FALSE]
    if(ar$regMethod == "PLS"){
      # pls.fit() centres the response and, with scale = TRUE, scales it
      XtYs <- XtY
      if(ar$scale) XtYs <- sweep(XtYs, 2, reg$Yscals, "/")
      project <- function(r){
        Sr <- as.numeric(S %*% r)
        list(t = NULL, tsq = sum(r * Sr), Xtt = Sr)
      }
      B <- .kernelpls(XtYs, project, ncomp, ncomp)$B
    } else {
      # svdspc.fit(): B = V diag(1/d^2) V'X'Y with the uncentred response
      e <- eigen(S, symmetric = TRUE)
      V <- e$vectors[, seq_len(ncomp), drop = FALSE]
      B <- V %*% (crossprod(V, XtY)/e$values[seq_len(ncomp)])
    }
    B <- matrix(B, length(genes), ncol(XtY), dimnames = list(genes, colnames(XtY)))
    stats <- list(XtX = S, XtY = XtY)
  } else {
    B <- .coefSlice(reg$coefficients, ncomp)[idx, , drop = FALSE]
    stats <- NULL
  }

  reg$coefficients <- array(B, c(dim(B), 1),
                            dimnames = c(dimnames(B), list(paste(ncomp, "comps"))))
  if(!is.null(reg$Xmeans)) reg$Xmeans <- reg$Xmeans[idx]
  if(!is.null(reg$Xscals)) reg$Xscals <- reg$Xscals[idx]
  reg$scores <- NULL
  reg$loadings <- NULL
  ar$reg_re <- reg
  ar$selectedFeat <- genes

  model$atlas_re <- ar
  model$selectedFeat <- genes
  model$impScores <- B
  model$stats <- stats
  model
}


## Cross-products of the centred and scaled cell by gene matrix X, with the
## centring and scaling of the fit reg_re: X'X and X'Y.
.crossStats <- function(X, YY, reg_re){

  X <- .fit_matrix(X)
  YY <- as.matrix(YY)
  mu <- if(is.null(reg_re$Xmeans)) rep(0, ncol(X)) else reg_re$Xmeans
  s <- if(is.null(reg_re$Xscals)) rep(1, ncol(X)) else reg_re$Xscals

  XtX <- (.crossprodBlocks(X) - nrow(X) * tcrossprod(mu))/tcrossprod(s)
  XtY <- (as.matrix(crossprod(X, YY)) - outer(mu, colSums(YY)))/s
  dimnames(XtX) <- list(colnames(X), colnames(X))
  dimnames(XtY) <- list(colnames(X), colnames(YY))

  list(XtX = XtX, XtY = XtY)
}


## X'X for a dense or sparse matrix X. A sparse X is multiplied in blocks of
## rows made dense (about blockBytes each): a dense crossprod() uses BLAS and
## was 13 times faster than the sparse product for a 58,759 x 4,613 matrix
## with 21% non-zero entries.
.crossprodBlocks <- function(X, blockBytes = 2^28){

  if(!inherits(X, "sparseMatrix")) return(crossprod(X))
  size <- max(1L, floor(blockBytes/(8 * ncol(X))))
  out <- NULL
  for(start in seq(1, nrow(X), by = size)){
    rows <- start:min(nrow(X), start + size - 1)
    block <- crossprod(as.matrix(X[rows, , drop = FALSE]))
    out <- if(is.null(out)) block else out + block
  }
  out
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
