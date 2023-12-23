#' Compute different errors for cross-validation.
#'
#' @param ncompGrid Vector.
#' @param XX Matrix.
#' @param YY Matrix.
#' @param XXtest Matrix.
#' @param YYtest Matrix.
#' @param labPerSample Integer.
#' @param selectedFeat Vector.
#' @param regMethod Character.
#' @param assayName Character.
#' @param mode Character.
#' @param normYY Logic.
#'
#' @return A list
getErr <- function(ncompGrid,
                   XX, YY, XXtest, YYtest,
                   labPerSample = NULL,
                   selectedFeat = NULL,
                   regMethod = 'PLS',
                   assayName = 'rank',
                   mode = "supervised",
                   normYY = F)
{

  ## Rank transform again after feature selection
  if(!is.null(selectedFeat)){
    selectedFeat <- selectedFeat
  } else {
    selectedFeat <- colnames(XX)
  }

  if(assayName == 'rank'){
    XX_RT <- RTassay(XX[, selectedFeat])
    XXtest_RT <- RTassay(XXtest[, selectedFeat])
  } else {
    XX_RT <- XX[, selectedFeat]
    XXtest_RT <- XXtest[, selectedFeat]
  }

  Nselect <- ncol(XX_RT)

  ## Regression using largest ncomp in grid (for selection)
  ncomp <- max(ncompGrid)
  precond_re2 <- mvr(XX_RT, YY,
                     ncomp = min(ncomp, Nselect),
                     method = regMethod)
  BhatAll2 <- precond_re2$coefficients

  ## Calculate errors
  if(ncomp <= Nselect){
    ncompV <- ncompGrid
  } else {
    ncompV <- seq(1, Nselect, length.out = length(ncompGrid))
  }
  YYtest <- as.matrix(YYtest)


  if(normYY){
    classOrigin <- getClass(normPhiScores(YYtest), labPerSample = labPerSample)
  } else {
    classOrigin <- getClass(YYtest, labPerSample = labPerSample)
  }

  if(mode == "supervised"){
    out <-
      sapply(ncompV,
             function(ncomp){
               Bhat2 <- BhatAll2[,,ncomp]
               XXtestCent <- scale(XXtest_RT,
                                   scale = F,
                                   center = colMeans(XX_RT))
               Yhat <- scale(XXtestCent %*% Bhat2,
                             scale = F,
                             center = -colMeans(as.matrix(YY)))

               classQuery <- getClass(Yhat, labPerSample)
               classQueryNorm <- getClass(normPhiScores(Yhat), labPerSample)


               Ers <- c(norm(abs(YYtest-Yhat), type = 'F')/norm(abs(YYtest), type = 'F'),
                        norm(abs(YYtest-Yhat), type = 'I'),
                        classErr(classQueryNorm, classOrigin, labPerSample)$err,
                        classErr(classQuery, classOrigin, labPerSample)$err
               )
               names(Ers) <- c("Frob ratio",
                               "inf norm",
                               'overall_norm', 'balanced_norm',
                               'overall', 'balanced')

               return(Ers)
             }
      )
  } else {
    out <-
      sapply(ncompV,
             function(ncomp){
               Bhat2 <- BhatAll2[,,ncomp]
               XXtestCent <- scale(XXtest_RT,
                                   scale = F,
                                   center = colMeans(XX_RT))
               Yhat <- scale(XXtestCent %*% Bhat2,
                             scale = F,
                             center = -colMeans(YY))

               Ers <- c(norm(abs(YYtest-Yhat), type = 'F')/norm(abs(YYtest), type = 'F'),
                        norm(abs(YYtest-Yhat), type = 'I'))
               names(Ers) <- c("Frob ratio",
                               "inf norm")

               return(Ers)
             }
      )
  }


  return(out)
}
