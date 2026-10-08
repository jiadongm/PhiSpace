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

  temp <- assay(sce, assayname)
  # RTassay ranks within columns, i.e. the genes within each cell
  temp <- RTassay(temp)

  if(sparse){
    temp <- as.sparse.matrix(temp)
  }
  assay(sce, targetAssay) <- temp
  return(sce)
}

