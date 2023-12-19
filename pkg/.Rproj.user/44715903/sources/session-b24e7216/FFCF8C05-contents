#' Calculate errors for the cross-validation for selecting Nselect.
#'
#' @param NselectV
#' @param orderByCol
#' @param ncomp
#' @param XX
#' @param YY
#' @param XXtest
#' @param YYtest
#' @param labPerSample
#' @param regMethod
#' @param assayName
#' @param mode
#' @param normYY
#'
#' @return
getErrNselect <- function(NselectV,
                          orderByCol,
                          ncomp,
                          XX, YY, XXtest, YYtest,
                          labPerSample = NULL,
                          regMethod = 'PLS',
                          assayName = 'logcounts',
                          mode = "supervised",
                          normYY = F)
{

  ## By default each sample has one label, eg cell type
  if(is.null(labPerSample)){
    labPerSample <- 1
  }

  ## Calculate errors
  YY <- as.matrix(YY)
  YYtest <- as.matrix(YYtest)
  if(normYY){
    classOrigin <- getClass(normPhiScores(YYtest), labPerSample)
  } else {
    classOrigin <- getClass(YYtest, labPerSample)
  }

  out <-
    sapply(NselectV,
           function(Nselect){

             selectFeat <- unique(
               as.vector(
                 orderByCol[1:Nselect, ]
               )
             )

             ## Prepare predictor matrix
             if(assayName == 'rank'){
               XX_RT <- RTassay(XX[, selectFeat])
               XXtest_RT <- RTassay(XXtest[, selectFeat])
             } else {
               XX_RT <- XX[, selectFeat]
               XXtest_RT <- XXtest[, selectFeat]
             }

             ## Regression using Nselect
             precond_re2 <- mvr(XX_RT, YY,
                                ncomp = ncomp,
                                method = regMethod)
             Bhat2 <- precond_re2$coefficients[,,ncomp]

             ## Prediction
             XXtestCent <- scale(XXtest_RT,
                                 scale = F,
                                 center = colMeans(XX_RT))
             Yhat <- scale(XXtestCent %*% Bhat2,
                           scale = F,
                           center = -colMeans(as.matrix(YY)))

             if(mode == "supervised"){

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
             } else {

               Ers <- c(norm(abs(YYtest-Yhat), type = 'F')/norm(abs(YYtest), type = 'F'),
                        norm(abs(YYtest-Yhat), type = 'I'))
               names(Ers) <- c("Frob ratio",
                               "inf norm")
             }


             return(Ers)
           }
    )


  return(out)
}
