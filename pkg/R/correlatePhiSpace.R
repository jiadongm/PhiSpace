#' Correlate PhiSpace scores with a response variable
#'
#' This function performs regression or classification analysis to identify
#' which PhiSpace cell state scores are most associated with a response variable.
#' It supports both continuous (regression) and discrete (classification) responses
#' using various methods including PLS, PLSDA, and DWD.
#'
#' @param spe A SpatialExperiment or SingleCellExperiment object containing
#'   PhiSpace scores in reducedDim
#' @param response Either a character string specifying a column name in colData(spe),
#'   or a vector of response values. For classification methods (PLSDA, DWD),
#'   this should be a factor or will be coerced to a factor.
#' @param method Character string specifying the analysis method. One of:
#'   \itemize{
#'     \item "PLS" - Partial Least Squares regression (for continuous response)
#'     \item "PLSDA" - Partial Least Squares Discriminant Analysis (for discrete response)
#'     \item "DWD" - Distance Weighted Discrimination (for binary classification)
#'   }
#' @param ncomp Integer specifying the number of components to use.
#'   If NULL (default), will be set automatically based on the method:
#'   \itemize{
#'     \item For PLSDA: number of classes - 1
#'     \item For PLS: min(10, ncol(PhiSpace)/2)
#'     \item For DWD: not applicable (DWD finds a single discriminant direction)
#'   }
#' @param reducedDimName Character string specifying which reducedDim slot
#'   contains the PhiSpace scores. Default is "PhiSpace".
#' @param center Logical indicating whether to center the PhiSpace scores
#'   before analysis. Default is TRUE.
#' @param scale Logical indicating whether to scale the PhiSpace scores
#'   before analysis. Default is FALSE.
#' @param dwd_params List of additional parameters for DWD (only used if method="DWD"):
#'   \itemize{
#'     \item kernel: kernel function (default: vanilladot())
#'     \item qval: q-value parameter for DWD (default: 1)
#'     \item lambda: regularization parameter or sequence (default: auto-tuned via CV)
#'     \item cv_folds: number of cross-validation folds (default: 5)
#'   }
#' @param seed Integer seed for reproducibility. Default is NULL (no seed set).
#'
#' @return A list with class "PhiSpaceCorrelation" containing:
#'   \item{method}{The method used}
#'   \item{importance_scores}{Matrix of feature importance scores (coefficients or weights)}
#'   \item{scores}{Predicted scores or component scores for each observation}
#'   \item{model}{The fitted model object}
#'   \item{feature_ranking}{Data frame with features ranked by importance}
#'   \item{response_summary}{Summary information about the response variable}
#'   \item{parameters}{List of parameters used in the analysis}
#'
#' @details
#' \strong{Method Selection:}
#' \itemize{
#'   \item Use \strong{PLS} for continuous responses (e.g., survival time, gene expression)
#'   \item Use \strong{PLSDA} for multi-class classification (e.g., disease subtypes, clusters)
#'   \item Use \strong{DWD} for binary classification when interpretability is important
#' }
#'
#' \strong{Feature Importance:}
#' \itemize{
#'   \item For PLS/PLSDA: regression coefficients from the final component
#'   \item For DWD: the discriminant weights
#' }
#'
#' The importance scores indicate which cell state features (PhiSpace dimensions)
#' are most strongly associated with the response. Positive values indicate
#' positive association, negative values indicate negative association.
#'
#' @examples
#' \dontrun{
#' # PLSDA example: associate PhiSpace scores with clusters
#' result <- correlatePhiSpace(
#'   spe = query_spe,
#'   response = "PhiClust",
#'   method = "PLSDA"
#' )
#'
#' # View top important features
#' head(result$feature_ranking, 10)
#'
#' # Plot importance scores
#' plot(result, nfeat = 20)
#'
#' # PLS example: associate with a continuous variable
#' result <- correlatePhiSpace(
#'   spe = query_spe,
#'   response = "spatial_distance",
#'   method = "PLS",
#'   ncomp = 5
#' )
#'
#' # DWD example: binary classification
#' result <- correlatePhiSpace(
#'   spe = query_spe,
#'   response = "cancer_vs_normal",
#'   method = "DWD",
#'   dwd_params = list(qval = 1, cv_folds = 10)
#' )
#' }
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @export
correlatePhiSpace <- function(
    spe,
    response,
    method = c("PLSDA", "PLS", "DWD"),
    ncomp = NULL,
    reducedDimName = "PhiSpace",
    center = TRUE,
    scale = FALSE,
    dwd_params = list(),
    seed = NULL
) {

  # Match and validate method
  method <- match.arg(method)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract PhiSpace scores
  if (!reducedDimName %in% reducedDimNames(spe)) {
    stop("reducedDim '", reducedDimName, "' not found in the object. ",
         "Available reducedDims: ", paste(reducedDimNames(spe), collapse = ", "))
  }
  XX <- reducedDim(spe, reducedDimName)

  # Extract response variable
  if (is.character(response) && length(response) == 1) {
    if (!response %in% colnames(colData(spe))) {
      stop("Response variable '", response, "' not found in colData. ",
           "Available columns: ", paste(colnames(colData(spe)), collapse = ", "))
    }
    YY_raw <- colData(spe)[[response]]
    response_name <- response
  } else {
    if (length(response) != nrow(XX)) {
      stop("Length of response vector (", length(response),
           ") does not match number of observations (", nrow(XX), ")")
    }
    YY_raw <- response
    response_name <- "user_provided"
  }

  # Remove any NA values
  valid_idx <- !is.na(YY_raw)
  if (sum(!valid_idx) > 0) {
    warning("Removing ", sum(!valid_idx), " observations with NA response values")
    XX <- XX[valid_idx, , drop = FALSE]
    YY_raw <- YY_raw[valid_idx]
  }

  # Validate and prepare response based on method
  if (method == "DWD") {
    # DWD requires binary response
    YY_factor <- as.factor(YY_raw)
    if (nlevels(YY_factor) != 2) {
      stop("DWD requires a binary response variable. Found ",
           nlevels(YY_factor), " levels: ",
           paste(levels(YY_factor), collapse = ", "))
    }
    # Convert to -1/1 encoding for DWD
    YY <- as.numeric(YY_factor)
    YY <- 2 * (YY - 1.5)  # Convert 1,2 to -1,1
    response_summary <- list(
      type = "binary",
      levels = levels(YY_factor),
      n_per_class = table(YY_factor)
    )

  } else if (method == "PLSDA") {
    # PLSDA requires factor response
    YY_factor <- as.factor(YY_raw)
    YY <- codeY_vec(YY_factor)
    response_summary <- list(
      type = "categorical",
      levels = levels(YY_factor),
      n_classes = nlevels(YY_factor),
      n_per_class = table(YY_factor)
    )

  } else {  # PLS
    # PLS can handle continuous response
    if (is.factor(YY_raw) || is.character(YY_raw)) {
      warning("Response appears to be categorical but method='PLS' was specified. ",
              "Consider using method='PLSDA' for classification.")
    }
    YY <- as.numeric(YY_raw)
    response_summary <- list(
      type = "continuous",
      range = range(YY),
      mean = mean(YY),
      sd = stats::sd(YY)
    )
  }

  # Set default ncomp if not provided
  if (is.null(ncomp)) {
    if (method == "PLSDA") {
      ncomp <- ncol(YY)  # Number of classes for PLSDA
    } else if (method == "PLS") {
      ncomp <- min(10, floor(ncol(XX) / 2))
    }
    # DWD doesn't use ncomp
  }

  # Fit model based on method
  if (method == "DWD") {
    result <- .fit_dwd(XX, YY, center, scale, dwd_params)

  } else {  # PLS or PLSDA
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
    reducedDimName = reducedDimName,
    dwd_params = if (method == "DWD") dwd_params else NA
  )

  class(result) <- c("PhiSpaceCorrelation", "list")
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
