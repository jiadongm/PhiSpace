#' Apply rank transform to a gene by cell matrix.
#'
#' @param X A gene by cell matrix.
#'
#' @return Rank transformed cell by gene matrix.
#'
RTassay <- function(X){

  temp <- Matrix::Matrix(
    apply(X, 2, rank, ties.method = "min") - 1,
    sparse = TRUE
  )
  temp <- temp/(nrow(temp) - 1)
  return(temp)
}

#' Calcualte classification errors.
#'
#' @param classQuery Character vector.
#' @param classOriginal Character vector.
#' @param labPerSample Integer.
#'
#' @return A list containing:
#' \item{err}{}
#' \item{byClassErrs}{}
#'
classErr <- function(classQuery, classOriginal, labPerSample = NULL){

  if(is.null(labPerSample)) labPerSample <- 1

  if(labPerSample == 1){

    out_errs <- rep(NA, 2)

    # Overall error
    out_errs[1] <- sum(classQuery != classOriginal)/length(classOriginal)

    # Balanced error
    classL <- unique(classOriginal)
    errs <- sapply(1:length(classL),
                   function(x){
                     cl <- classL[x]
                     idx <- (classOriginal == cl)
                     # Per class classification error
                     sum(classQuery[idx] != classOriginal[idx])/sum(idx)
                   })
    names(errs) <- classL
    out_errs[2] <- mean(errs)

    return(list(err = out_errs, byClassErrs = sort(errs)))

  } else {

    if(labPerSample <= 0) stop("labPerSample has to be positive integer.")

    # Overall errors
    ove_errs <-
      sapply(1:nrow(classQuery),
             function(x){
               out <- sum(2 - sum(classQuery[x,] %in% classOriginal[x,]))/labPerSample
             })
    ove_errs <- sum(ove_errs)/length(ove_errs)

    # Balanced error
    out_list <- vector("list", labPerSample)
    for(erIdx in 1:ncol(classOriginal)){

      classL <- unique(classOriginal[,erIdx])
      errs <- sapply(1:length(classL),
                     function(x){
                       cl <- classL[x]
                       idx <- (classOriginal[,erIdx] == cl)
                       sum(classQuery[idx, erIdx] != classOriginal[idx, erIdx])/sum(idx)
                     })
      names(errs) <- classL
      out_list[[erIdx]] <- errs

    }

    byClassErrs <- do.call(`c`, out_list)
    BER <- mean(byClassErrs)

    return(
      list(
        err = c(ove_errs, BER),
        byClassErrs = byClassErrs
      )
    )

  }
}


#' Create partitions of data for cross-validation.
#'
#' @param x Vector.
#' @param n Integer.
#'
#' @return A partition of index vector `x` to `n` folds.
#'
split2 <- function(x, n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}





#' Turn a matrix to a data.frame withe newly defined colnames.
#'
#' @param mat Matrix.
#' @param key String.
#'
#' @return Data frame.
#' @export
reNameCols <- function(
  mat,
  key = "comp"
){

  mat <- as.data.frame(mat)
  colnames(mat) <- paste0(key, 1:ncol(mat))

  return(mat)
}


## Center and scale by matrix operations
scal <- function(X, center = NULL, scale = NULL){

  ones <- rep(1, nrow(X))
  if(!is.null(center)){

    X <- X - outer(ones, center)
  }

  if(!is.null(scale)){

    X <- X / outer(ones, scale)
  }

  return(X)
}



## Double centring
#' Double center a matrix by column and row means, resulting in a new matrix with zero row and column means.
#'
#' @param X A matrix or an object convertible to one.
#'
#' @return A double-centred matrix.
#' @export
doubleCent <- function(X){

  X <- as.matrix(X)
  X <- sweep(X, 1, rowMeans(X))
  X <- sweep(X, 2, colMeans(X))
  return(X)
}


