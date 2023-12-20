#' Relative counts normalisation.
#'
#' @param sce
#' @param assayName
#' @param targetAssay
#'
#' @return A SingleCellExperiment object containing a new assay after relative counts normalisation.
#'
#' @export
rcTransf <- function(sce, assayName = "counts", targetAssay = "data"){

  libSizes <- colSums(assay(sce, assayName))
  sce <- sce[,libSizes>0]
  libSizes <- libSizes[libSizes > 0]
  assay(sce, targetAssay) <- t(
    t(assay(sce, assayName))/libSizes
  )

  return(sce)
}
