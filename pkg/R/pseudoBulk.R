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
      out <- split(
        sample(idxList[[x]], nPool*resampSize, replace = TRUE),
        rep(1:resampSize, rep(nPool, resampSize))
      )
      names(out) <- NULL
      return(out)
    }
  )
  resampIdx <- do.call(c, resampIdx)
  # Aggregate X
  XX <- assay(sce, assayName)
  XXagg <- sapply(
    resampIdx,
    function(x){

      rowMeans(XX[,x])
    }
  )
  # Aggregate Y
  YYagg <- t(
    sapply(
      resampIdx,
      function(x){

        colMeans(YY[x,])
      }
    )
  )
  colnames(XXagg) <- rownames(YYagg) <- paste0("PB", 1:ncol(XXagg))

  # Output
  SingleCellExperiment(
    list(data = XXagg),
    reducedDims = list(
      response = YYagg
    )
  )


}
