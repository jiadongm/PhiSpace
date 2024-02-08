#' Run scran normalisation
#'
#' @param sce SCE object.
#' @param assayName Assay name.
#' @param clusts Cluster labels (optional).
#'
#' @return SCE object.
#' @export
scranTransf <- function(
    sce,
    assayName = "counts",
    clusts = NULL
){

  if(is.null(clusts)) clusts <- scran::quickCluster(sce, assay.type = assayName)
  sce <- scran::computeSumFactors(sce, assay.type = assayName, cluster = clusts)
  # If having non-positive size factors
  if(min(sce$sizeFactor) <= 0){
    smallestPositive <- sce$sizeFactor[which(sort(sce$sizeFactor) > 0)[1]]
    sce$sizeFactor[sce$sizeFactor <= 0] <- smallestPositive
  }
  sce <- scuttle::logNormCounts(sce, assay.type = assayName)

  return(sce)
}
