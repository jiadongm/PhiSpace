#' Center log ratio normalisation
#'
#' @param sce SCE object
#' @param assayName Assay to normalise
#' @param targetAssay Name of normalised assay
#' @param offset Offset value added to raw counts
#'
#' @return SCE object
#' @export
CLRnorm <- function(sce, assayName = "counts", targetAssay = "data", offset = 1){

  rawCounts <- assay(sce, assayName)
  logCounts <- log(rawCounts + offset)
  geoMeans <- rowMeans(logCounts)
  assay(sce, targetAssay) <- logCounts - geoMeans

  return(sce)
}
