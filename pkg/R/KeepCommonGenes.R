#' Keep common genes of two SingleCellExperiment objects.
#'
#' @param reference SingleCellExperiment object.
#' @param query SingleCellExperiment object.
#' @param dots2hyphen Logic.
#'
#' @return A list containing updated SingleCellExperiment obejcts with matching features.
#'
#' @export
KeepCommonGenes <- function(reference, query, dots2hyphen = F){
  if(dots2hyphen){
    rownames(reference) <- gsub("\\.", "-", rownames(reference))
    rownames(query) <- gsub("\\.", "-", rownames(query))
  }
  commonGenes <- intersect(rownames(query), rownames(reference))
  reference <- reference[commonGenes, ]
  query <- query[commonGenes, ]
  return(list(reference=reference, query=query))
}
