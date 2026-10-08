MATLAB_cols <- c(rgb(54, 70, 157, maxColorValue = 255),
                 rgb(61, 146, 185, maxColorValue = 255),
                 rgb(126, 203, 166, maxColorValue = 255),
                 rgb(204, 234, 156, maxColorValue = 255),
                 rgb(249, 252, 181, maxColorValue = 255),
                 rgb(255, 226, 144, maxColorValue = 255),
                 rgb(253, 164, 93, maxColorValue = 255),
                 rgb(234, 95, 70, maxColorValue = 255),
                 rgb(185, 30, 72, maxColorValue = 255))

censor <- function(vec, quant = 0.05){

  vec <- as.vector(vec)

  cutoff <- stats::quantile(vec, 1-quant)
  vec[vec <= cutoff] <- cutoff
  vec
}


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



## Predictor matrix for model fitting: sparse input stays sparse (as a
## dgCMatrix); any other input becomes a base matrix.
.fit_matrix <- function(X){

  if(inherits(X, "sparseMatrix")){
    X <- methods::as(methods::as(methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  } else {
    X <- as.matrix(X)
  }

  return(X)
}


## Column standard deviations. For a dgCMatrix they are computed from the
## non-zero entries and the number of zeros, without a dense copy.
.colSds <- function(X, means = colMeans(X)){

  if(inherits(X, "CsparseMatrix")){

    nnz <- diff(X@p)
    dev2 <- X
    dev2@x <- (X@x - rep(means, nnz))^2
    ss <- colSums(dev2) + (nrow(X) - nnz) * means^2
    out <- sqrt(ss/(nrow(X) - 1))
    names(out) <- colnames(X)
  } else {

    out <- apply(X, 2, stats::sd)
  }

  return(out)
}


## Coefficient matrix for ncomp components from a coefficient array
## (features x responses x stored components). Works for arrays that store
## every component (older objects) and for arrays that store only some.
.coefSlice <- function(coefs, ncomp){

  if(length(dim(coefs)) != 3) return(as.matrix(coefs))

  sliceNames <- dimnames(coefs)[[3]]
  sliceName <- paste(ncomp, "comps")
  if(!is.null(sliceNames) && sliceName %in% sliceNames){
    k <- match(sliceName, sliceNames)
  } else if(is.null(sliceNames) && dim(coefs)[3] >= ncomp){
    k <- ncomp
  } else {
    stop("Coefficients for ", ncomp, " components are not stored.")
  }

  out <- coefs[, , k]
  dim(out) <- dim(coefs)[1:2]
  dimnames(out) <- dimnames(coefs)[1:2]

  return(out)
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
  X <- sweep(X, 1, rowMeans(X, na.rm = T))
  X <- sweep(X, 2, colMeans(X, na.rm = T))
  return(X)
}


