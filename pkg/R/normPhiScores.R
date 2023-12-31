#' Normalise PhiSpace phenotypic scores by row or by column.
#'
#' The by-row normalisation is the same as `SingleR:::.trim_byLabel_and_normalize_scores`.
#'
#' @param X Matrix.
#' @param method Character.
#'
#' @return A matrix normalised by row or by column.
#'
#' @export
normPhiScores <- function(X, method = c("col", "row")){
  method <- match.arg(method)

  X <- as.matrix(X)

  if(method == "col"){
    Xout <- apply(X, 2,
                  function(x){
                    x_cent <- x - stats::median(x)
                    x_cent/max(abs(x_cent))
                  })
    dimnames(Xout) <- dimnames(X)
  } else {
    # This is how SingleR normalise their scores for better visual effects
    # See SingleR:::.trim_byLabel_and_normalize_scores
    row_mins <- apply(X, 1, min)
    row_maxs <- apply(X, 1, max)
    Xout <- (X - row_mins)/pmax(row_maxs - row_mins, 1e-08)
    Xout <- Xout^3
  }

  return(Xout)
}
