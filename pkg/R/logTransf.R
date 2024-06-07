#' Relative counts normalisation.
#'
#' @param sce SingleCellExperiment object.
#' @param assayName Character.
#' @param targetAssay Character.
#' @param use_log1p Logical. Use shifted log (log1p) transform, using raw counts instead of relative counts.
#' @param scalFactor Scale factor to inflate relative counts.
#' @param sparse Store the transformed assay as sparse matrix or not.
#'
#' @return A SingleCellExperiment object containing a new assay after relative counts normalisation.
#'
#' @export
logTransf <- function(
    sce,
    scalFactor = 10000,
    assayName = "counts",
    targetAssay = "data",
    use_log1p = FALSE,
    sparse = TRUE
){

  libSizes <- colSums(assay(sce, assayName))
  sce <- sce[,libSizes>0]
  libSizes <- libSizes[libSizes > 0]

  if(use_log1p){

    assay(sce, targetAssay) <-
      Matrix::Matrix(
        log1p(
          assay(sce, assayName)
        ),
        sparse = sparse
      )
  } else {

    assay(sce, targetAssay) <-
      Matrix::Matrix(
        log(
          scale(
            assay(sce, assayName),
            center = F,
            scale = libSizes
          ) * scalFactor + 1
        ),
        sparse = sparse
      )
  }

  return(sce)
}
