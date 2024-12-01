#' Create cell type correlation matrices
#'
#' A generic function to compute a correlation matrix for different input types, including `matrix`,
#' `data.frame` and `SingleCellExperiment` objects.
#'
#' @param input The input data. Can be a `matrix`, `data.frame`, or `SingleCellExperiment`.
#' @param mask A zero-one matrix with the same number of rows as the input. Optional. If provided, each column of `mask` defines a subset of rows used to compute a separate correlation matrix.
#' @param reducedDimName A character string specifying the name of the reduced dimension to use for `SingleCellExperiment` objects. Defaults to `"PhiSpace"`.
#' @param groupBy A column of `colData(input)` when input is an SCE object. If provided, mask will be defined according to `colData(input)[,groupBy]`.
#'
#' @return If `mask` is provided, returns a list of correlation matrices (one per column of `mask`).
#' Otherwise, returns a single correlation matrix.
#'
#' @export
setGeneric(
  "cellTypeCorMat",
  function(input,
           mask = NULL,
           reducedDimName = "PhiSpace",
           groupBy = NULL){
    standardGeneric("cellTypeCorMat")
  }
)



#' @rdname cellTypeCorMat
#' @aliases cellTypeCorMat,matrix-method
#' @param input A `matrix` containing numeric data for correlation computation.
#' @param mask A zero-one matrix with the same number of rows as the input. Optional. If provided,
#' each column of `mask` defines a subset of rows used to compute a separate correlation matrix.
#' @details This method calculates the correlation matrix directly for the input `matrix`.
#'
#' @export
setMethod(
  "cellTypeCorMat", "matrix",
  function(input,
           mask = NULL,
           reducedDimName = "PhiSpace",
           groupBy = NULL) {

    if (!is.null(mask)) {

      mask <- as.matrix(mask)

      if (!all(mask %in% c(0, 1)) || nrow(mask) != nrow(input)) {
        stop("mask must be a zero-one matrix with the same number of rows as input.")
      }
      # Compute correlation matrices for each column of the mask
      result <- lapply(
        seq_len(ncol(mask)),
        function(i){
          subset_rows <- (mask[, i] == 1)
          if(sum(subset_rows) == 0) stop(paste0("The ", i,"th column of mask doesn't contain any 1s."))
          if(sum(subset_rows) < 10) warning(paste0("The ", i,"th column of mask contains fewer than 10 1s."))
          stats::cor(input[subset_rows, , drop = FALSE], use = "pairwise.complete.obs")
        })
      return(result)
    } else {
      # Compute single correlation matrix
      return(stats::cor(input, use = "pairwise.complete.obs"))
    }
  }
)



#' @rdname cellTypeCorMat
#' @aliases cellTypeCorMat,data.frame-method
#' @details This method first converts the `data.frame` to a `matrix` and then calculates the correlation matrix.
#'
#' @export
setMethod(
  "cellTypeCorMat", "data.frame",
  function(input,
           mask = NULL,
           reducedDimName = "PhiSpace",
           groupBy = NULL) {
    cellTypeCorMat(as.matrix(input), mask)
  }
)





#' @rdname cellTypeCorMat
#' @aliases cellTypeCorMat,SingleCellExperiment-method
#' @details This method uses the reduced dimension `"PhiSpace"` of the `SingleCellExperiment` object to compute the correlation matrix.
#' If `"PhiSpace"` is not available, an error is raised.
#'
#' @export
setMethod(
  "cellTypeCorMat", "SingleCellExperiment",
  function(input,
           mask = NULL,
           reducedDimName = "PhiSpace",
           groupBy = NULL) {

    if (!reducedDimName %in% reducedDimNames(input)) {
      stop(paste(reducedDimName, "is not a valid reduced dimension in the SCE object."))
    }

    if(!is.null(groupBy)){

      if(!is.null(mask)) warning("Both mask and groupBy were provided, will use groupBy.")
      if(!(groupBy %in% colnames(colData(input)))) stop("groupBy is not a column of colData(input).")

      mask <- codeY_vec(colData(input)[,groupBy], method = "0,1")
    }

    data <- reducedDim(input, reducedDimName)
    cellTypeCorMat(data, mask)
  }
)
