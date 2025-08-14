#' Find spatial niches by clustering PhiSpace scores
#'
#' Identifies spatial niches by performing k-means clustering on PhiSpace cell type scores.
#' Optionally performs PCA before clustering to reduce dimensionality and improve clustering performance.
#' Can test multiple numbers of niches and return all results for comparison.
#'
#' @param spe SpatialExperiment or SingleCellExperiment object containing PhiSpace scores.
#' @param reducedDimName Name of the reducedDim containing PhiSpace scores (default: "PhiSpace").
#' @param n_niches Integer or vector of integers. Number(s) of niches (clusters) to identify.
#'   If a vector is provided, multiple clustering results will be generated.
#' @param use_pca Logical. Whether to perform PCA before clustering (default: TRUE).
#' @param ncomp Integer. Number of principal components to use if use_pca = TRUE.
#'   If NULL, defaults to half the number of cell types in PhiSpace scores.
#' @param center Logical. Whether to center data for PCA (default: TRUE).
#' @param scale Logical. Whether to scale data for PCA (default: FALSE).
#' @param kmeans_iter Integer. Maximum number of iterations for k-means (default: 500).
#' @param kmeans_nstart Integer. Number of random starts for k-means (default: 20).
#' @param kmeans_algorithm Character. Algorithm for k-means clustering (default: "Lloyd").
#' @param seed Integer. Random seed for reproducibility (default: 94863).
#' @param store_pca Logical. Whether to store PCA results in the object (default: FALSE).
#' @param pca_name Character. Name for stored PCA results if store_pca = TRUE (default: "PhiSpace_PCA").
#'
#' @return The input object with niche assignments added to colData. If a single n_niches is provided,
#'   results are stored as "spatial_niches". If multiple n_niches are provided, results are stored as
#'   "spatial_niches_k7", "spatial_niches_k8", etc. If store_pca = TRUE, PCA results are also stored in reducedDims.
#'
#' @details
#' This function identifies spatial niches by clustering cells/spots based on their PhiSpace
#' cell type composition scores. The workflow consists of:
#'
#' 1. Extract PhiSpace scores from the specified reducedDim
#' 2. Optionally perform PCA to reduce dimensionality and focus on major patterns
#' 3. Apply k-means clustering to identify niches (for each k in n_niches)
#' 4. Store niche assignments in colData
#'
#' PCA is recommended when you have many cell types (>20) as it can improve clustering
#' by focusing on the main axes of variation in cell type composition.
#'
#' When multiple values of n_niches are provided, the function will test each one and store
#' all results, allowing you to compare different numbers of niches and choose the most
#' appropriate one for your analysis.
#'
#' @examples
#' \dontrun{
#' # Basic usage - find 9 niches using PCA
#' spe <- findNiches(spe, n_niches = 9)
#'
#' # Test multiple numbers of niches
#' spe <- findNiches(spe, n_niches = c(7, 8, 9, 10))
#'
#' # Find niches without PCA
#' spe <- findNiches(spe, n_niches = 6, use_pca = FALSE)
#'
#' # Custom PCA settings
#' spe <- findNiches(spe, n_niches = c(6, 8, 10), ncomp = 15, store_pca = TRUE)
#'
#' # Access niche assignments
#' table(spe$spatial_niches)  # for single n_niches
#' table(spe$spatial_niches_k8)  # for multiple n_niches
#'
#' # Visualize niches spatially
#' VizSpatial(spe, groupBy = "spatial_niches")
#' VizSpatial(spe, groupBy = "spatial_niches_k8")
#' }
#'
#' @export
findNiches <- function(spe,
                       reducedDimName = "PhiSpace",
                       n_niches,
                       use_pca = TRUE,
                       ncomp = NULL,
                       center = TRUE,
                       scale = FALSE,
                       kmeans_iter = 500L,
                       kmeans_nstart = 20L,
                       kmeans_algorithm = "Lloyd",
                       seed = 94863,
                       store_pca = FALSE,
                       pca_name = "PhiSpace_PCA") {

  # Input validation
  if (!inherits(spe, c("SpatialExperiment", "SingleCellExperiment"))) {
    stop("Input must be a SpatialExperiment or SingleCellExperiment object")
  }

  if (!reducedDimName %in% reducedDimNames(spe)) {
    stop(paste("reducedDim", reducedDimName, "not found in object"))
  }

  if (missing(n_niches) || !is.numeric(n_niches) || any(n_niches < 1)) {
    stop("n_niches must be a positive integer or vector of positive integers")
  }

  # Sort n_niches for consistent output
  n_niches <- sort(unique(n_niches))
  multiple_k <- length(n_niches) > 1

  if (multiple_k) {
    cat("Testing", length(n_niches), "different numbers of niches:", paste(n_niches, collapse = ", "), "\n")
  }

  # Extract PhiSpace scores
  phi_scores <- reducedDim(spe, reducedDimName)

  if (nrow(phi_scores) != ncol(spe)) {
    stop("Number of rows in PhiSpace scores must match number of cells/spots")
  }

  cat("Found PhiSpace scores with", ncol(phi_scores), "cell types for", nrow(phi_scores), "cells/spots\n")

  # Determine matrix for clustering
  if (use_pca) {
    # Set default ncomp if not provided
    if (is.null(ncomp)) {
      ncomp <- max(2, min(30, floor(ncol(phi_scores) / 2)))
    }

    # Validate ncomp
    if (ncomp >= ncol(phi_scores)) {
      warning("ncomp >= number of cell types, setting ncomp to ncol(phi_scores) - 1")
      ncomp <- ncol(phi_scores) - 1
    }

    cat("Performing PCA with", ncomp, "components\n")

    # Perform PCA
    pca_result <- getPC(X = phi_scores,
                        ncomp = ncomp,
                        center = center,
                        scale = scale)

    mat2clust <- pca_result$scores

    # Store PCA results if requested
    if (store_pca) {
      reducedDim(spe, pca_name) <- mat2clust
      cat("PCA results stored as", pca_name, "\n")
    }

  } else {
    cat("Using PhiSpace scores directly for clustering\n")
    mat2clust <- phi_scores
  }

  # Validate n_niches
  if (any(n_niches >= nrow(mat2clust))) {
    invalid_k <- n_niches[n_niches >= nrow(mat2clust)]
    stop("All values in n_niches must be less than the number of cells/spots. Invalid values: ",
         paste(invalid_k, collapse = ", "))
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Perform clustering for each k
  cluster_results <- list()
  metadata_entries <- list()

  for (k in n_niches) {
    cat("Performing k-means clustering with k =", k, "\n")

    # Perform k-means clustering
    cluster_result <- stats::kmeans(
      x = mat2clust,
      centers = k,
      iter.max = kmeans_iter,
      nstart = kmeans_nstart,
      algorithm = kmeans_algorithm
    )

    # Create niche labels
    niche_labels <- factor(cluster_result$cluster,
                           levels = 1:k,
                           labels = paste0("Niche_", 1:k))

    # Store results
    cluster_results[[as.character(k)]] <- niche_labels

    # Create column name
    if (multiple_k) {
      col_name <- paste0("spatial_niches_k", k)
    } else {
      col_name <- "spatial_niches"
    }

    # Store in colData
    colData(spe)[[col_name]] <- niche_labels

    # Print distribution
    cat("Niche distribution for k =", k, ":\n")
    print(table(niche_labels))
    cat("\n")

    # Store metadata
    metadata_entry <- list(
      k = k,
      column_name = col_name,
      within_ss = cluster_result$withinss,
      total_ss = cluster_result$totss,
      between_ss = cluster_result$betweenss,
      tot_within_ss = cluster_result$tot.withinss
    )

    metadata_entries[[as.character(k)]] <- metadata_entry
  }

  # Add metadata about the analysis
  metadata_entry <- list(
    function_used = "findNiches",
    n_niches_tested = n_niches,
    use_pca = use_pca,
    reducedDimName = reducedDimName,
    seed = seed,
    clustering_results = metadata_entries,
    timestamp = Sys.time()
  )

  if (use_pca) {
    metadata_entry$ncomp = ncomp
    metadata_entry$center = center
    metadata_entry$scale = scale
    metadata_entry$variance_explained = sum(pca_result$props)
  }

  S4Vectors::metadata(spe)$nicheAnalysis <- c(
    S4Vectors::metadata(spe)$nicheAnalysis,
    list(metadata_entry)
  )

  if (multiple_k) {
    cat("Spatial niches identified for", length(n_niches), "different k values\n")
    cat("Results stored in colData as:", paste0("spatial_niches_k", n_niches, collapse = ", "), "\n")
  } else {
    cat("Spatial niches identified and stored in colData$spatial_niches\n")
  }

  return(spe)
}
