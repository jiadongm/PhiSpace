#' Minimal quality control: remove features with all zero values.
#'
#' @param sce
#' @param assayName
#'
#' @return
#' @export
zeroFeatQC <- function(sce, assayName = "counts"){

  geneSpars <- rowMeans(assay(sce, assayName) == 0)
  if(max(geneSpars) == 0){
    sce <- sce[geneSpars < 1, ]
    cat("Deleted features with all zeros.")
  } else {
    cat("All features have at least 1 nonzero value.")
  }

  return(sce)
}
