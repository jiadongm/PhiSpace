#' Subsample from a sce object
#'
#' @param sce SCE object to be sampled from
#' @param key Required for stratefied subsampling
#' @param proportion Proportion of cells to sample
#' @param seed Random seed
#'
#' @return A downsized SCE object
#' @export
subsample <- function(sce, key = NULL, proportion = 0.1, seed = 5202056){

  if(is.null(key)){

    sceSize <- ncol(sce)
    outSize <- max(1, ceiling(sceSize * proportion) )

    set.seed(seed)
    idx <- sample(
      1:ncol(sce),
      outSize
    )

  } else {


    # Stratefied subsampling
    keys <- as.character(colData(sce)[,key])
    keyNames <- unique(keys)

    set.seed(seed)
    idx_list <- lapply(
      keyNames,
      function(keyName){

        fullIdx <- which(keys == keyName)
        fullIdxSize <- length(fullIdx)
        subIdxSize <- max(1, ceiling(fullIdxSize * proportion) )
        subIdx <- sample(fullIdx, subIdxSize)

        return(subIdx)
      }
    )
    idx <- unlist(idx_list)

  }


  return(sce[,idx])
}
