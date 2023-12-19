#' Convert a matrix to sparse metrix.
#'
#' @param X
#'
#' @return
#' @importFrom methods as
as.sparse.matrix <- function(X){
  X <- as(
    object = as(object = as(object = X, Class = "dMatrix"),
                Class = "generalMatrix"),
    Class = "CsparseMatrix"
  )
  return(X)
}
