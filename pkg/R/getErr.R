#' Compute different errors for cross-validation.
#'
#' @param ncompGrid Vector.
#' @param XXtrain Matrix.
#' @param YYtrain Matrix.
#' @param XXtest Matrix.
#' @param YYtest Matrix.
#' @param selectedFeat Vector.
#' @param regMethod Character.
#' @param assayName Character.
#' @param center Logic.
#' @param scale Logic.
#'
#' @return A list
getErr <- function(ncompGrid,
                   XXtrain, YYtrain, XXtest, YYtest,
                   center = TRUE,
                   scale = FALSE,
                   selectedFeat = NULL,
                   regMethod = 'PLS',
                   assayName = 'rank')
{

  if(is.null(selectedFeat)) selectedFeat <- colnames(XXtrain)


  if(assayName == 'rank'){
    XX_RT <- RTassay(XXtrain[, selectedFeat])
    XXtest_RT <- RTassay(XXtest[, selectedFeat])
  } else {
    XX_RT <- XXtrain[, selectedFeat]
    XXtest_RT <- XXtest[, selectedFeat]
  }

  nfeat <- ncol(XX_RT)

  XX_RT_cent <- scale(XX_RT, center = center, scale = scale)

  ## Regression using largest ncomp in grid (for selection)
  ncomp <- max(ncompGrid)
  BhatAll <- mvr(XX_RT_cent, YYtrain,
                 ncomp = min(ncomp, nfeat),
                 method = regMethod,
                 center = FALSE, scale = FALSE)$coefficients

  ## Calculate errors
  if(ncomp <= nfeat){
    ncompV <- ncompGrid
  } else {
    ncompV <- seq(1, nfeat, length.out = length(ncompGrid))
  }
  YYtest <- as.matrix(YYtest)



  out <-
    sapply(ncompV,
           function(ncomp){

             Bhat <- BhatAll[,,ncomp]

             if(scale){
               scal <- attr(XX_RT_cent, "scaled:scale")
             } else {
               scal <- FALSE
             }

             XXtestCent <- scale(XXtest_RT,
                                 scale = scal,
                                 center = attr(XX_RT_cent, "scaled:center") )
             Yhat <- scale(XXtestCent %*% Bhat,
                           scale = F,
                           center = -colMeans(YYtrain))

             Ers <- norm(abs(YYtest-Yhat), type = 'F')/norm(abs(YYtest), type = 'F')

             return(Ers)
           }
    )


  return(out)
}
