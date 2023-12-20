#' Dicretise soft annotation to labels.
#'
#' @param X
#' @param labPerSample
#'
#' @return A vector contaning labels.
#'
#' @export
getClass <- function(X, labPerSample = NULL){

  if(is.null(labPerSample)) labPerSample <- 1

  if(labPerSample <= 0) stop("labPerSample has to be postive integer.")

  if(labPerSample == 1){

    out <- colnames(X)[apply(X, 1, which.max)]

  } else {

    out <- sapply(1:nrow(X),
                  function(x){
                    vec <- X[x,,drop=F]
                    ct_sort <- colnames(vec)[order(-vec)]
                    ct_sort[1:labPerSample]
                    # colnames(X)[rank(-vec, ties.method = "last") <= labPerSample]
                  })
    out <- t(out)

  }


  return(out)
}
