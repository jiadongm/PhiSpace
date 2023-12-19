#' Keep common genes of two SingleCellExperiment objects.
#'
#' @param reference
#' @param query
#' @param dots2hyphen
#'
#' @return
#' @export
keepCommonGenes <- function(reference, query, dots2hyphen = F){
  if(dots2hyphen){
    rownames(reference) <- gsub("\\.", "-", rownames(reference))
    rownames(query) <- gsub("\\.", "-", rownames(query))
  }
  commonGenes <- intersect(rownames(query), rownames(reference))
  reference <- reference[commonGenes, ]
  query <- query[commonGenes, ]
  return(list(reference=reference, query=query))
}
