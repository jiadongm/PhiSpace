#' Convert a matrix to sparse metrix.
#'
#' @param X Matrix.
#'
#' @return Sparse matrix version of `X`
as.sparse.matrix <- function(X){
  X <- Matrix::Matrix(
    X, sparse = TRUE
  )
  return(X)
}
