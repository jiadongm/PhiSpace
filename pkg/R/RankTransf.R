#' Rank transform of a `SingleCellExperiment` object.
#'
#' @param sce `SingleCellExperiment` object.
#' @param assayname Character.
#' @param targetAssay Character.
#' @param sparse Logic.
#'
#' @return An updated `SingleCellExperiment` object with a rank transformed assay.
#'
#' @export
RankTransf <- function(sce, assayname = 'counts', targetAssay = 'rank', sparse = TRUE){

  temp <- Matrix::t(assay(sce, assayname))
  # RTassay takes cell by gene matrix as input
  temp <- Matrix::t(RTassay(temp))

  if(sparse){
    temp <- as.sparse.matrix(temp)
  }
  assay(sce, targetAssay) <- temp
  return(sce)
}

