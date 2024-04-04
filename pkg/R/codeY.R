#' Convert discrete annotations to dummy variable matrices.
#'
#' @param sce SingleCellExperiment object.
#' @param phenotypes Which phenotype annotations (e.g. cell types) to convert.
#' @param method Which coding to use, either (-1,1) or (0,1).
#'
#' @return
#' \item{YY}{A dummy variable matrix encoding discrete annotation}
#'
#' @export
codeY <- function(sce, phenotypes, method = c("-1,1", "0,1")){

  method <- match.arg(method)

  YY <- vector("list", length(phenotypes))

  for(phIdx in 1:length(phenotypes)){
    ph <- phenotypes[phIdx]
    labs <- unique(as.character(colData(sce)[,ph]))
    YY_temp <- sapply(1:length(labs),
                      function(x){
                        lab <- labs[x]
                        out <- as.numeric(colData(sce)[,ph] == lab)
                        if(method == "-1,1") out[out == 0] <- -1
                        out
                      })
    colnames(YY_temp) <- labs
    YY[[phIdx]] <- YY_temp
  }
  YY <- do.call("cbind", YY)
  rownames(YY) <- colnames(sce)

  return(YY)
}


#' Convert discrete annotations to dummy variable matrices.
#'
#' @param vec Factor contanining class labels.
#' @param rowNames Optional. Row names for output dummy matrix (eg cell barcode).
#' @param method Which coding to use, either (-1,1) or (0,1).
#'
#' @return
#' \item{YY}{A dummy variable matrix encoding discrete annotation}
#'
#' @export
codeY_vec <- function(vec, rowNames = NULL, method = c("-1,1", "0,1")){

  vec <- as.factor(vec)

  method <- match.arg(method)

  labs <- levels(vec)
  YY <- sapply(
    1:length(labs),
    function(x){
      lab <- labs[x]
      out <- as.numeric(vec == lab)
      if(method == "-1,1") out[out == 0] <- -1
      out
    }
  )
  colnames(YY) <- labs

  if(!is.null(rowNames)) rownames(YY) <- rowNames

  return(YY)
}
