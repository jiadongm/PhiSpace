#' Calculate errors for the cross-validation for selecting nfeat.
#'
#' @param nfeatV Vector.
#' @param impScores Matrix.
#' @param ncomp Integer.
#' @param XXtrain Matrix.
#' @param YYtrain Matrix.
#' @param XXtest Matrix.
#' @param YYtest Matrix.
#' @param regMethod Character.
#' @param center Logic.
#' @param scale Logic.
#' @param assayName Character.
#' @param loss Type of loss function.
#'
#' @return A list.
getErr_nfeat <- function(
    nfeatV,
    impScores,
    ncomp,
    XXtrain, YYtrain, XXtest, YYtest,
    center = TRUE,
    scale = FALSE,
    regMethod = 'PLS',
    assayName = 'logcounts',
    loss = c("exponential", "Frobenius")
  )
{

  loss <- match.arg(loss)

  ## Calculate errors
  YYtrain <- as.matrix(YYtrain)
  YYtest <- as.matrix(YYtest)

  out <-
    sapply(nfeatV,
           function(nfeat){

             selectedFeat <- selectFeat(impScores, nfeat)

             ## Prepare predictor matrix
             if(assayName == 'rank'){
               XX_RT <- RTassay(XXtrain[, selectedFeat])
               XXtest_RT <- RTassay(XXtest[, selectedFeat])
             } else {
               XX_RT <- XXtrain[, selectedFeat]
               XXtest_RT <- XXtest[, selectedFeat]
             }

             XX_RT_cent <- scale(XX_RT, center = center, scale = scale)

             ## Regression using nfeat
             Bhat <- mvr(XX_RT_cent, YYtrain,
                         ncomp = ncomp,
                         method = regMethod,
                         center = FALSE, scale = FALSE)$coefficients[,,ncomp]

             ## Prediction
             if(scale){
               scal <- attr(XX_RT_cent, "scaled:scale")
             } else {
               scal <- FALSE
             }
             XXtestCent <- scale(XXtest_RT,
                                 scale = scal,
                                 center = attr(XX_RT_cent, "scaled:center"))
             Yhat <- scale(XXtestCent %*% Bhat,
                           scale = F,
                           center = -colMeans(as.matrix(YYtrain)))



             if(loss == "Frobenius"){

               # Frobenius norm
               Ers <- norm(abs(YYtest-Yhat), type = 'F')/norm(abs(YYtest), type = 'F')
             } else {

               # Exponential loss
               Ers <- 1/nrow(YYtest) * sum(
                 exp(
                   - YYtest * Yhat
                 )
               )
             }


             return(Ers)
           }
    )


  return(out)
}
