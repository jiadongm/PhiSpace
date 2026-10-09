#' Converting single cells to pseudo-bulks.
#'
#' @param sce `SingleCellExpeirment` object.
#' @param clusterid Column name of colData for specifying groups for pseudobulking.
#' @param phenotypes Phenotypes to predict. Can be multiple.
#' @param response Response matrix. Can be continuous.
#' @param assayName Assay used for pseudobulking.
#' @param resampSizes How many pseudo-bulk resamples to generate from each cluster.
#' @param proportion Used to determine unequal resampSizes.
#' @param seed Random seeds
#' @param nPool How many cells to use for pseudobulking.
#' @param calcMean Logical. Calculate mean of sum of expression for pseudobulking.
#'
#' @return An updated `SingleCellExpeirment` object with a pseudobulk assay.
#'
#' @export
pseudoBulk <- function(
    sce,
    phenotypes = NULL,
    response = NULL,
    clusterid = NULL,
    assayName = 'counts',
    resampSizes = 100,
    proportion = NULL,
    calcMean = FALSE,
    nPool = 15,
    seed = 904800
){

  ## Response matrix
  if(is.null(response)){

    if(is.null(phenotypes)) stop("Phenotypes and response cannot both be NULL.")

    YY <- codeY(sce, phenotypes)

    if(is.null(clusterid)) clusterid <- phenotypes[1]
  } else {

    YY <- response
  }

  if(is.null(clusterid)) stop("Need to specify clusterid.")

  ## Resampling indices
  nGenes <- nrow(sce)
  nCells <- ncol(sce)
  # Index list
  cluster <- as.character(colData(sce)[, clusterid])
  cluname <- sort(unique(cluster))
  idxList <- split(1:nCells, cluster)
  names(idxList) <- cluname
  Ncluster <- length(cluname)
  clustSizes <- sapply(idxList, length)
  # How many to resample
  if(is.null(proportion)){

    if(length(resampSizes) == 1){

      resampSizes <- rep(resampSizes, Ncluster)
    } else {

      if(length(resampSizes) != Ncluster) stop("Length of resampSizes has to be either 1 or number of clusters.")
    }

  } else {

    resampSizes <- ceiling(clustSizes * proportion)
  }
  # Resample
  set.seed(seed)
  resampIdx <- lapply(
    1:Ncluster,
    function(x){

      resampSize <- resampSizes[x]
      # Index into the cluster: sample() on a single index i would draw
      # from 1:i
      clustIdx <- idxList[[x]]
      out <- split(
        clustIdx[sample.int(length(clustIdx), nPool*resampSize, replace = TRUE)],
        rep(1:resampSize, rep(nPool, resampSize))
      )
      names(out) <- NULL
      return(out)
    }
  )
  resampIdx <- do.call(c, resampIdx)
  # Indicator matrix A (cells x pseudo-bulks): A[c, b] counts how often cell c
  # was drawn for pseudo-bulk b, so X %*% A gives the sums
  poolSizes <- lengths(resampIdx)
  A <- Matrix::sparseMatrix(
    i = unlist(resampIdx, use.names = FALSE),
    j = rep(seq_along(resampIdx), poolSizes),
    x = 1,
    dims = c(nCells, length(resampIdx))
  )
  # Aggregate X
  XX <- assay(sce, assayName)
  if (!is.matrix(XX) && !methods::is(XX, "Matrix")) XX <- as.matrix(XX)
  XXagg <- as.matrix(XX %*% A)
  if(calcMean) XXagg <- sweep(XXagg, 2, poolSizes, "/")
  dimnames(XXagg) <- list(rownames(XX), NULL)
  # Aggregate Y
  YYagg <- as.matrix(Matrix::crossprod(A, as.matrix(YY))) / poolSizes
  dimnames(YYagg) <- list(NULL, colnames(YY))
  colnames(XXagg) <- rownames(YYagg) <- paste0("PB", 1:ncol(XXagg))

  # Output
  sce <- SingleCellExperiment(list(counts = XXagg), reducedDims = list(response = YYagg))
  if(length(phenotypes)==1) colData(sce)[,phenotypes] <- getClass(YYagg)
  return(sce)
}
