#' Rank Features by Discriminative Power
#'
#' Uses supervised learning methods (PLS-DA, PLS, or DWD) to rank features by
#' their ability to discriminate between groups or predict a continuous response.
#' Based on PhiSpace's correlatePhiSpace but generalized to work with any
#' numeric data (not just PhiSpace scores).
#'
#' @param data A matrix or data.frame where rows are observations and columns
#'   are features. Can also be a SpatialExperiment, SingleCellExperiment, or
#'   SummarizedExperiment object.
#' @param response Either a character string specifying a column name in colData
#'   (when \code{data} is an experiment object), or a vector of response values.
#'   For classification (PLSDA, DWD), should be factor or character.
#'   For regression (PLS), should be numeric.
#' @param method Character string specifying the method: "PLSDA" for classification,
#'   "PLS" for regression, or "DWD" for binary classification. Default is "PLSDA".
#' @param source Character string specifying data source when \code{data} is an
#'   experiment object: "assay" or "reducedDim". Default is "reducedDim".
#' @param assay_name Character string specifying which assay to use when
#'   \code{source = "assay"}. Default is "logcounts".
#' @param reducedDim_name Character string specifying which reduced dimension
#'   to use when \code{source = "reducedDim"}. Default is "PhiSpace".
#' @param ncomp Integer specifying number of components. If NULL (default),
#'   set automatically: number of classes for PLSDA, min(10, ncol(data)/2) for PLS.
#'   Not used for DWD.
#' @param center Logical indicating whether to center features. Default is TRUE.
#' @param scale Logical indicating whether to scale features. Default is FALSE.
#' @param dwd_params List of parameters for DWD method (see Details).
#' @param seed Integer seed for reproducibility. Default is NULL.
#'
#' @return A list with class "FeatureRanking" containing:
#'   \item{method}{The method used}
#'   \item{importance_scores}{Matrix of feature importance scores}
#'   \item{scores}{Component or discriminant scores for observations}
#'   \item{model}{The fitted model object}
#'   \item{feature_ranking}{Data frame with features ranked by importance}
#'   \item{response_summary}{Summary of the response variable}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' \strong{Methods:}
#' \itemize{
#'   \item \strong{PLSDA}: Partial Least Squares Discriminant Analysis for
#'     multi-class classification. Features ranked by coefficient magnitude.
#'   \item \strong{PLS}: Partial Least Squares regression for continuous response.
#'     Features ranked by coefficient magnitude.
#'   \item \strong{DWD}: Distance Weighted Discrimination for binary classification.
#'     Features ranked by discriminant weights.
#' }
#'
#' \strong{DWD Parameters} (in \code{dwd_params} list):
#' \itemize{
#'   \item \code{kernel}: Kernel function (default: vanilladot())
#'   \item \code{qval}: q-value parameter (default: 1)
#'   \item \code{lambda}: Regularization parameter or sequence (default: auto-tuned)
#'   \item \code{cv_folds}: Cross-validation folds (default: 5)
#' }
#'
#' @examples
#' \dontrun{
#' # With matrix input - PLSDA
#' expr_mat <- t(assay(spe, "logcounts"))
#' result <- rankFeatures(expr_mat, response = spe$cluster, method = "PLSDA")
#' head(result$feature_ranking, 10)
#'
#' # With experiment object - from assay
#' result <- rankFeatures(spe, response = "cluster", method = "PLSDA",
#'                        source = "assay", assay_name = "logcounts")
#'
#' # From reduced dimensions (e.g., PhiSpace scores)
#' result <- rankFeatures(spe, response = "cluster", method = "PLSDA",
#'                        source = "reducedDim", reducedDim_name = "PhiSpace")
#'
#' # PLS for continuous response
#' result <- rankFeatures(expr_mat, response = spe$spatial_score, method = "PLS")
#'
#' # DWD for binary classification
#' result <- rankFeatures(expr_mat, response = spe$treatment, method = "DWD")
#' }
#'
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom stats sd
#' @importFrom utils modifyList
#' @export
rankFeatures <- function(data,
                         response,
                         method = c("PLSDA", "PLS", "DWD"),
                         source = c("reducedDim", "assay"),
                         assay_name = "logcounts",
                         reducedDim_name = "PhiSpace",
                         ncomp = NULL,
                         center = TRUE,
                         scale = FALSE,
                         dwd_params = list(),
                         seed = NULL) {

  method <- match.arg(method)
  source <- match.arg(source)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract feature matrix
  if (inherits(data, "SpatialExperiment") ||
      inherits(data, "SingleCellExperiment") ||
      inherits(data, "SummarizedExperiment")) {

    # Extract response from colData if character
    if (is.character(response) && length(response) == 1) {
      if (!response %in% colnames(SummarizedExperiment::colData(data))) {
        stop("Response '", response, "' not found in colData. ",
             "Available: ", paste(colnames(SummarizedExperiment::colData(data)), collapse = ", "))
      }
      YY_raw <- SummarizedExperiment::colData(data)[[response]]
      response_name <- response
    } else {
      if (length(response) != ncol(data)) {
        stop("Length of response (", length(response),
             ") does not match number of observations (", ncol(data), ")")
      }
      YY_raw <- response
      response_name <- "user_provided"
    }

    # Extract feature data
    if (source == "assay") {
      if (!assay_name %in% SummarizedExperiment::assayNames(data)) {
        stop("Assay '", assay_name, "' not found. Available: ",
             paste(SummarizedExperiment::assayNames(data), collapse = ", "))
      }
      # Transpose to get observations x features
      XX <- t(as.matrix(SummarizedExperiment::assay(data, assay_name)))
    } else {
      if (!reducedDim_name %in% SingleCellExperiment::reducedDimNames(data)) {
        stop("reducedDim '", reducedDim_name, "' not found. Available: ",
             paste(SingleCellExperiment::reducedDimNames(data), collapse = ", "))
      }
      XX <- SingleCellExperiment::reducedDim(data, reducedDim_name)
    }

  } else if (is.matrix(data) || is.data.frame(data)) {
    XX <- as.matrix(data)
    YY_raw <- response
    response_name <- "user_provided"

    if (length(YY_raw) != nrow(XX)) {
      stop("Length of response (", length(YY_raw),
           ") does not match number of observations (", nrow(XX), ")")
    }
  } else {
    stop("'data' must be a matrix, data.frame, or experiment object")
  }

  # Remove NAs
  valid_idx <- !is.na(YY_raw)
  if (sum(!valid_idx) > 0) {
    warning("Removing ", sum(!valid_idx), " observations with NA response")
    XX <- XX[valid_idx, , drop = FALSE]
    YY_raw <- YY_raw[valid_idx]
  }

  # Prepare response based on method
  if (method == "DWD") {
    YY_factor <- as.factor(YY_raw)
    if (nlevels(YY_factor) != 2) {
      stop("DWD requires binary response. Found ", nlevels(YY_factor), " levels")
    }
    YY <- as.numeric(YY_factor)
    YY <- 2 * (YY - 1.5)  # Convert to -1/1
    response_summary <- list(
      type = "binary",
      levels = levels(YY_factor),
      n_per_class = table(YY_factor)
    )
  } else if (method == "PLSDA") {
    YY_factor <- as.factor(YY_raw)
    YY <- codeY_vec(YY_factor)
    response_summary <- list(
      type = "categorical",
      levels = levels(YY_factor),
      n_classes = nlevels(YY_factor),
      n_per_class = table(YY_factor)
    )
  } else {  # PLS
    if (is.factor(YY_raw) || is.character(YY_raw)) {
      warning("Response appears categorical but method='PLS'. ",
              "Consider method='PLSDA' for classification")
    }
    YY <- as.numeric(YY_raw)
    response_summary <- list(
      type = "continuous",
      range = range(YY),
      mean = mean(YY),
      sd = stats::sd(YY)
    )
  }

  # Set default ncomp
  if (is.null(ncomp)) {
    if (method == "PLSDA") {
      ncomp <- ncol(YY)
    } else if (method == "PLS") {
      ncomp <- min(10, floor(ncol(XX) / 2))
    }
  }

  # Fit model
  if (method == "DWD") {
    result <- .fit_dwd(XX, YY, center, scale, dwd_params)
  } else {
    result <- .fit_pls(XX, YY, ncomp, center, scale, method)
  }

  # Add metadata
  result$method <- method
  result$response_name <- response_name
  result$response_summary <- response_summary
  result$parameters <- list(
    ncomp = if (method != "DWD") ncomp else NA,
    center = center,
    scale = scale,
    source = if (inherits(data, "SummarizedExperiment")) source else "matrix",
    assay_name = if (source == "assay") assay_name else NA,
    reducedDim_name = if (source == "reducedDim") reducedDim_name else NA,
    dwd_params = if (method == "DWD") dwd_params else NA
  )

  class(result) <- c("FeatureRanking", "list")
  return(result)
}


#' Internal function to fit PLS or PLSDA
#' @keywords internal
.fit_pls <- function(XX, YY, ncomp, center, scale, method) {

  # Fit PLS model using PhiSpace's internal mvr function
  pls_model <- mvr(
    X = XX,
    Y = YY,
    ncomp = ncomp,
    method = "PLS",
    center = center,
    scale = scale
  )

  # Extract coefficients from final component
  # The coefficients array is [features x responses x components]
  # We want [features x responses] for the final component
  if (length(dim(pls_model$coefficients)) == 3) {
    importance_scores <- pls_model$coefficients[, , ncomp, drop = FALSE]
    # Convert 3D array to 2D matrix
    dim(importance_scores) <- dim(importance_scores)[1:2]
    rownames(importance_scores) <- rownames(pls_model$coefficients)
    colnames(importance_scores) <- if (is.matrix(YY)) colnames(YY) else "Response"
  } else {
    # Handle case where coefficients might already be 2D
    importance_scores <- as.matrix(pls_model$coefficients)
    colnames(importance_scores) <- if (is.matrix(YY)) colnames(YY) else "Response"
  }

  # Get component scores
  scores <- pls_model$scores

  # Rank features by absolute importance
  if (ncol(importance_scores) == 1) {
    # Single response
    feature_ranking <- data.frame(
      feature = rownames(importance_scores),
      importance = importance_scores[, 1],
      abs_importance = abs(importance_scores[, 1]),
      stringsAsFactors = FALSE
    )
  } else {
    # Multiple responses (PLSDA case) - use Euclidean norm across responses
    feature_ranking <- data.frame(
      feature = rownames(importance_scores),
      importance_norm = sqrt(rowSums(importance_scores^2)),
      stringsAsFactors = FALSE
    )
    # Add individual response importances
    for (i in 1:ncol(importance_scores)) {
      feature_ranking[[paste0("importance_", colnames(importance_scores)[i])]] <-
        importance_scores[, i]
    }
  }

  # Sort by importance
  if ("abs_importance" %in% colnames(feature_ranking)) {
    feature_ranking <- feature_ranking[order(-feature_ranking$abs_importance), ]
  } else {
    feature_ranking <- feature_ranking[order(-feature_ranking$importance_norm), ]
  }
  rownames(feature_ranking) <- NULL

  return(list(
    importance_scores = importance_scores,
    scores = scores,
    model = pls_model,
    feature_ranking = feature_ranking
  ))
}


#' Internal function to fit DWD
#'
#' @importFrom stats predict
#' @keywords internal
.fit_dwd <- function(XX, YY, center, scale, dwd_params) {

  # Check if kerndwd is available
  if (!requireNamespace("kerndwd", quietly = TRUE)) {
    stop("Package 'kerndwd' is required for DWD method. Please install it with:\n",
         "  install.packages('kerndwd')")
  }

  # Set default DWD parameters
  default_params <- list(
    kernel = kerndwd::vanilladot(),
    qval = 1,
    lambda = 10^(seq(3, -3, length.out = 50)),
    cv_folds = 5,
    eps = 1e-5,
    maxit = 1e5
  )

  # Override with user parameters
  dwd_params <- utils::modifyList(default_params, dwd_params)

  # Center and/or scale data
  if (center || scale) {
    XX <- scale(XX, center = center, scale = scale)
  }

  # Perform cross-validation if lambda is a sequence
  if (length(dwd_params$lambda) > 1) {
    cv_result <- kerndwd::cv.kerndwd(
      x = XX,
      y = YY,
      kern = dwd_params$kernel,
      qval = dwd_params$qval,
      lambda = dwd_params$lambda,
      eps = dwd_params$eps,
      maxit = dwd_params$maxit,
      nfolds = dwd_params$cv_folds
    )
    optimal_lambda <- cv_result$lambda.min
  } else {
    optimal_lambda <- dwd_params$lambda
    cv_result <- NULL
  }

  # Fit final DWD model
  dwd_model <- kerndwd::kerndwd(
    x = XX,
    y = YY,
    kern = dwd_params$kernel,
    qval = dwd_params$qval,
    lambda = optimal_lambda,
    eps = dwd_params$eps,
    maxit = dwd_params$maxit
  )

  # Extract discriminant weights (loadings)
  A <- dwd_model$alpha[-1, , drop = FALSE]
  importance_scores <- crossprod(XX, A)
  colnames(importance_scores) <- "DWD_weight"

  # Compute DWD scores
  scores <- predict(dwd_model, dwd_params$kernel, XX, XX, type = "link")
  colnames(scores) <- "DWD_score"

  # Rank features by absolute importance
  feature_ranking <- data.frame(
    feature = rownames(importance_scores),
    importance = importance_scores[, 1],
    abs_importance = abs(importance_scores[, 1]),
    stringsAsFactors = FALSE
  )
  feature_ranking <- feature_ranking[order(-feature_ranking$abs_importance), ]
  rownames(feature_ranking) <- NULL

  return(list(
    importance_scores = importance_scores,
    scores = scores,
    model = dwd_model,
    feature_ranking = feature_ranking,
    cv_result = cv_result,
    optimal_lambda = optimal_lambda
  ))
}
