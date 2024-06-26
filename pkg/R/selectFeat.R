#' Select nfeat top features from each column (label) of impScores
#'
#' @param impScores Matrix. Each column corresponds to a phenotype (eg a cell type) and contains importance scores of features for predicting that phenotype.
#' @param nfeat Number of top features to select from each column of impScores.
#' @param absVal Take absolute value or not.
#'
#' @return Vector of selected features.
#' @export
selectFeat <- function(
    impScores,
    nfeat = 5,
    absVal = TRUE
  ){

  impScores <- as.matrix(impScores)

  if(is.null(rownames(impScores))) stop("impScores has to have row names corresponding to feature names.")

  orderByCol <- apply(
    impScores,
    2,
    function(x){
      names(x) <- rownames(impScores)

      if(absVal){
        toSort <- abs(x)
      } else {
        toSort <- x
      }

      names(sort(toSort, decreasing = T))
    }
  )

  selectFeat <- unique(
    as.vector(
      orderByCol[1:nfeat, ]
    )
  )

  list(
    orderedFeatMat = orderByCol,
    selectedFeat = selectFeat
  )
}
