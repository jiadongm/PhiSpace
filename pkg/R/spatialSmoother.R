#' Spatial smoother wrapper for SpatialExperiment and SingleCellExperiment objects
#'
#' @param object SpatialExperiment or SingleCellExperiment object
#' @param assay2smooth Name of assay to smooth (for gene expression smoothing)
#' @param smoothedAssay Name for the new smoothed assay (default: paste0(assay2smooth, "_smoothed"))
#' @param smoothReducedDim Logical, whether to smooth reduced dimensions instead of gene expression
#' @param reducedDim2smooth Name of reduced dimension to smooth (default: "PhiSpace")
#' @param smoothedReducedDim Name for the new smoothed reduced dimension (default: paste0(reducedDim2smooth, "_smoothed"))
#' @param x_coord Column name in colData for x coordinates (required for SCE objects)
#' @param y_coord Column name in colData for y coordinates (required for SCE objects)
#' @param k Number of nearest neighbors (default: 10)
#' @param kernel Kernel function: "gaussian", "uniform", or "linear" (default: "gaussian")
#' @param sigma Bandwidth parameter for Gaussian kernel (default: auto-computed)
#' @param include_self Whether to include the cell itself in smoothing (default: TRUE)
#' @return Modified object with new smoothed assay or reduced dimension
spatialSmoother <- function(object,
                            assay2smooth = "logcounts",
                            smoothedAssay = NULL,
                            smoothReducedDim = FALSE,
                            reducedDim2smooth = "PhiSpace",
                            smoothedReducedDim = NULL,
                            x_coord = NULL,
                            y_coord = NULL,
                            k = 10,
                            kernel = "linear",
                            sigma = NULL,
                            include_self = TRUE) {


  # Check object type and extract coordinates
  is_spe <- inherits(object, "SpatialExperiment")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_spe && !is_sce) {
    stop("Object must be a SpatialExperiment or SingleCellExperiment")
  }

  # Extract spatial coordinates
  if (is_spe) {
    # For SpatialExperiment objects, use spatialCoords()
    coords_matrix <- spatialCoords(object)
    if (ncol(coords_matrix) < 2) {
      stop("SpatialExperiment object must have at least 2 spatial coordinates")
    }
    coordinates <- data.frame(
      x = coords_matrix[, 1],
      y = coords_matrix[, 2]
    )
    cat("Using spatial coordinates from SpatialExperiment object\n")

  } else if (is_sce) {
    # For SingleCellExperiment objects, use specified columns from colData
    if (is.null(x_coord) || is.null(y_coord)) {
      stop("For SingleCellExperiment objects, x_coord and y_coord must be specified")
    }

    if (!x_coord %in% colnames(colData(object))) {
      stop(paste("Column", x_coord, "not found in colData"))
    }
    if (!y_coord %in% colnames(colData(object))) {
      stop(paste("Column", y_coord, "not found in colData"))
    }

    coordinates <- data.frame(
      x = colData(object)[[x_coord]],
      y = colData(object)[[y_coord]]
    )
    cat("Using coordinates from colData columns:", x_coord, "and", y_coord, "\n")
  }

  # Check for missing coordinates
  if (any(is.na(coordinates$x)) || any(is.na(coordinates$y))) {
    stop("Missing values detected in spatial coordinates")
  }

  # Decide what to smooth
  if (smoothReducedDim) {
    # Set default name for smoothed reduced dimension if not provided
    if (is.null(smoothedReducedDim)) {
      smoothedReducedDim <- paste0(reducedDim2smooth, "_smoothed")
    }

    # Smooth reduced dimensions
    cat("Smoothing reduced dimensions:", reducedDim2smooth, "\n")

    if (!reducedDim2smooth %in% reducedDimNames(object)) {
      stop(paste("Reduced dimension", reducedDim2smooth, "not found in object"))
    }

    # Get reduced dimension matrix (cells x dimensions)
    redDim_matrix <- reducedDim(object, reducedDim2smooth)

    # Transpose for smoothing function (dimensions x cells)
    redDim_t <- t(redDim_matrix)

    # Smooth the reduced dimensions
    smoothed_redDim_t <- spatial_smooth_expression(
      expression_matrix = redDim_t,
      coordinates = coordinates,
      k = k,
      kernel = kernel,
      sigma = sigma,
      include_self = include_self
    )

    # Transpose back (cells x dimensions)
    smoothed_redDim <- t(smoothed_redDim_t)

    # Add smoothed reduced dimension to object
    reducedDim(object, smoothedReducedDim) <- smoothed_redDim

    cat("Added smoothed reduced dimension:", smoothedReducedDim, "\n")

  } else {
    # Set default name for smoothed assay if not provided
    if (is.null(smoothedAssay)) {
      smoothedAssay <- paste0(assay2smooth, "_smoothed")
    }

    # Smooth gene expression
    cat("Smoothing gene expression from assay:", assay2smooth, "\n")

    if (!assay2smooth %in% assayNames(object)) {
      stop(paste("Assay", assay2smooth, "not found in object"))
    }

    # Get expression matrix
    expr_matrix <- assay(object, assay2smooth)

    # Smooth the expression matrix
    smoothed_expr <- spatial_smooth_expression(
      expression_matrix = expr_matrix,
      coordinates = coordinates,
      k = k,
      kernel = kernel,
      sigma = sigma,
      include_self = include_self
    )

    # Add smoothed assay to object
    assay(object, smoothedAssay) <- smoothed_expr

    cat("Added smoothed assay:", smoothedAssay, "\n")
  }

  # Add metadata about smoothing parameters
  metadata_entry <- list(
    function_used = "spatialSmoother",
    smoothReducedDim = smoothReducedDim,
    k = k,
    kernel = kernel,
    sigma = sigma,
    include_self = include_self,
    timestamp = Sys.time()
  )

  if (smoothReducedDim) {
    metadata_entry$reducedDim2smooth <- reducedDim2smooth
    metadata_entry$output_name <- smoothedReducedDim
  } else {
    metadata_entry$assay2smooth <- assay2smooth
    metadata_entry$output_name <- smoothedAssay
  }

  if (is_sce && !is_spe) {
    metadata_entry$x_coord <- x_coord
    metadata_entry$y_coord <- y_coord
  }

  # Add to metadata
  S4Vectors::metadata(object)$spatialSmoothing <- c(
    S4Vectors::metadata(object)$spatialSmoothing, list(metadata_entry)
  )

  return(object)
}





#' Spatial smoothing of gene expression using KNN with kernel weighting
#'
#' @param expression_matrix Gene x Cell matrix of expression values
#' @param coordinates Data frame with columns 'x' and 'y' for cell coordinates
#' @param k Number of nearest neighbors (default: 10)
#' @param kernel Kernel function: "gaussian", "uniform", or "linear" (default: "gaussian")
#' @param sigma Bandwidth parameter for Gaussian kernel (default: auto-computed)
#' @param include_self Whether to include the cell itself in smoothing (default: TRUE)
#' @return Smoothed gene x cell expression matrix
spatial_smooth_expression <- function(
    expression_matrix,
    coordinates,
    k = 10,
    kernel = "linear",
    sigma = NULL,
    include_self = TRUE
) {

  # Input validation
  if (ncol(expression_matrix) != nrow(coordinates)) {
    stop("Number of cells in expression matrix must match number of coordinate rows")
  }

  if (!all(c("x", "y") %in% colnames(coordinates))) {
    stop("Coordinates must have 'x' and 'y' columns")
  }

  n_cells <- ncol(expression_matrix)
  n_genes <- nrow(expression_matrix)

  # Adjust k if including self
  search_k <- if(include_self) k else k + 1

  # Find k-nearest neighbors for each cell
  cat("Finding k-nearest neighbors...\n")
  knn_result <- FNN::get.knnx(data = coordinates[, c("x", "y")],
                              query = coordinates[, c("x", "y")],
                              k = search_k)

  knn_indices <- knn_result$nn.index
  knn_distances <- knn_result$nn.dist

  # Remove self if not including it
  if (!include_self) {
    # Remove the first column (self) from indices and distances
    knn_indices <- knn_indices[, -1, drop = FALSE]
    knn_distances <- knn_distances[, -1, drop = FALSE]
  }

  # Compute kernel weights
  cat("Computing kernel weights...\n")
  weights <- compute_kernel_weights(knn_distances, kernel, sigma)

  # Initialize smoothed expression matrix
  smoothed_expr <- matrix(0, nrow = n_genes, ncol = n_cells)
  rownames(smoothed_expr) <- rownames(expression_matrix)
  colnames(smoothed_expr) <- colnames(expression_matrix)

  # Perform smoothing for each cell
  cat("Performing spatial smoothing...\n")
  pb <- utils::txtProgressBar(min = 0, max = n_cells, style = 3)

  for (i in 1:n_cells) {
    # Get neighbors and weights for current cell
    neighbors <- knn_indices[i, ]
    cell_weights <- weights[i, ]

    # Normalize weights to sum to 1
    cell_weights <- cell_weights / sum(cell_weights)

    # Compute weighted average for all genes
    smoothed_expr[, i] <- as.vector(expression_matrix[, neighbors] %*% cell_weights)

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  cat("Spatial smoothing completed!\n")
  return(smoothed_expr)
}



#' Compute kernel weights based on distances
#'
#' @param distances Matrix of distances (cells x neighbors)
#' @param kernel Kernel type: "gaussian", "uniform", or "linear"
#' @param sigma Bandwidth parameter (auto-computed if NULL)
#' @return Matrix of weights
compute_kernel_weights <- function(distances, kernel, sigma = NULL) {

  if (kernel == "gaussian") {
    # Auto-compute sigma if not provided (median of all distances)
    if (is.null(sigma)) {
      sigma <- stats::median(distances[distances > 0])
    }
    # Gaussian kernel: exp(-d^2 / (2*sigma^2))
    weights <- exp(-distances^2 / (2 * sigma^2))

  } else if (kernel == "uniform") {
    # Uniform kernel: all neighbors get equal weight
    weights <- matrix(1, nrow = nrow(distances), ncol = ncol(distances))

  } else if (kernel == "linear") {
    # Linear kernel: weight = 1 - d/max_d (inverse distance)
    max_dist <- apply(distances, 1, max)
    weights <- 1 - distances / max_dist
    weights[weights < 0] <- 0  # Ensure non-negative weights

  } else {
    stop("Supported kernels: 'gaussian', 'uniform', 'linear'")
  }

  return(weights)
}


