#' Converting single cells to pseudo-bulks.
#'
#' @param query
#' @param clusterid
#' @param assay
#' @param size.factor
#' @param pool.factor
#' @param nPool
#' @param YY
#'
#' @return An updated `SingleCellExpeirment` object with a pseudobulk assay
#'
#' @export
pseudoBulk <- function(query,
                       clusterid,
                       assay = 'counts',
                       size.factor = 1,
                       pool.factor = 1,
                       nPool = 15,
                       YY = NULL)
{
  cluster <- colData(query)[, clusterid]
  cluname <- sort(unique(cluster))
  Ncluster <- length(cluname)
  G <- nrow(query)

  if(length(size.factor)==1){
    sizeFacVec <- rep(size.factor, Ncluster)
  } else if (length(size.factor)==Ncluster) {
    sizeFacVec <- size.factor[cluname]
  } else {
    stop("Length of size.factor has to be either 1 or same as number of clusters/celltypes.")
  }

  outList <- outYYlist <- vector("list", length(cluname))
  outCluster <- c()

  if(!is.null(YY)){
    if(nrow(YY)!=ncol(query)) stop("nrow of YY has to be equal to #cells in query.")
  }

  #Aggregate cells
  for(i in cluname){
    message(paste("\r Now aggregate", i), appendLF = F)
    subQuery <- Matrix::Matrix(assay(query, assay)[, cluster == i], sparse = T)
    if(!is.null(YY)){
      subYY <- YY[cluster == i, , drop=F]
    }
    nCell <- ncol(subQuery)
    nOut <- ceiling(nCell * sizeFacVec[i])
    if(is.null(nPool)) nPool <- ceiling(nCell * pool.factor)

    subOut <- matrix(NA, nrow = G, ncol = nOut)
    colnames(subOut) <- paste(i, '.', 1:nOut, sep = '')
    if(!is.null(YY)){
      subYYout <- matrix(NA, nrow = nOut, ncol = ncol(subYY))
      colnames(subYYout) <- colnames(YY)
    }



    sampleIndex <- sample(1:nCell, nOut * nPool, replace = T) %>%
      matrix(ncol = nPool, nrow = nOut)

    if(nPool > 1){
      for(j in 1:nOut){
        subOut[,j] <- rowSums(subQuery[,sampleIndex[j,]])/nPool
        if(!is.null(YY)){
          subYYout[j,] <- colSums(subYY[sampleIndex[j,], ])/nPool
        }
      }
    }else{
      for(j in 1:nOut){
        subOut[,j] <- subQuery[,sampleIndex[j,]]/nPool
        if(!is.null(YY)){
          subYYout[j,] <- subYY[sampleIndex[j,], ]/nPool
        }
      }
    }

    outList[[i]] <- subOut
    if(!is.null(YY)){
      colnames(subYYout) <- colnames(YY)
      outYYlist[[i]] <- subYYout
    }
    outCluster <- c(outCluster, rep(i, nOut))
  }

  #Collapse the list
  outList <- unname(outList)
  outData <- do.call(cbind, outList)
  rownames(outData) <- rownames(query)

  #remove duplicated pseudo-bulk samples.
  # dup <- duplicated(t(outData)) | duplicated(t(outData)[ncol(outData):1,])[ncol(outData):1]
  # outData <- outData[,!dup]
  # outCluster <- outCluster[!dup]


  sce <- SingleCellExperiment(list(data = outData))
  colData(sce)[,clusterid] <- outCluster
  if(!is.null(YY)){
    outYYlist <- unname(outYYlist)
    outYY <- do.call(rbind, outYYlist)
    rownames(outYY) <- colnames(outData)
    colData(sce) <- cbind(colData(sce), outYY)
  }



  message(paste('\n Sparsity after imputation is', round(mean(outData == 0),3)))
  return(sce)

}
