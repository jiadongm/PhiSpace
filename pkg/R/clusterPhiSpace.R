#' K-means clustering on PhiSpace principal components
#'
#' Performs k-means clustering on principal components derived from PhiSpace
#' scores or any other matrix-like data. This is useful for identifying spatial
#' niches or clusters based on cell state compositions or gene expression patterns.
#'
#' @param x Either a SpatialExperiment/SingleCellExperiment object containing
#'   PhiSpace scores in reducedDim, OR a matrix-like object (matrix, sparse matrix,
#'   data frame) with features in columns and observations in rows.
#' @param k Integer specifying the number of clusters. Either provide k or
#'   k_range but not both.
#' @param k_range Integer vector of length 2 specifying range of k values to
#'   test (e.g., c(5, 15)). Will use elbow method or silhouette to select
#'   optimal k. Either provide k or k_range but not both.
#' @param select_k_method Character string specifying method to select optimal k
#'   when k_range is provided. Options: "elbow" (total within-cluster SS) or
#'   "silhouette" (average silhouette width). Default is "silhouette".
#' @param ncomp Integer specifying number of principal components to use for
#'   clustering. If NULL (default), uses min(30, nfeatures - 1). If "all",
#'   uses nfeatures - 1 components.
#' @param reducedDimName Character string specifying which reducedDim slot
#'   contains the PhiSpace scores. Only used when x is a SingleCellExperiment
#'   or SpatialExperiment. Default is "PhiSpace".
#' @param use_assay Character string specifying which assay to use if extracting
#'   data from a SingleCellExperiment/SpatialExperiment object instead of using
#'   reducedDim. If NULL (default), uses reducedDim specified by reducedDimName.
#'   Common options: "logcounts", "counts", "normcounts".
#' @param nstart Integer specifying number of random starts for k-means.
#'   Default is 20.
#' @param iter.max Integer specifying maximum number of iterations.
#'   Default is 500.
#' @param algorithm Character string specifying k-means algorithm. Options are
#'   "Hartigan-Wong", "Lloyd", "Forgy", "MacQueen". Default is "Lloyd".
#' @param center Logical indicating whether to center data before PCA.
#'   Default is TRUE.
#' @param scale Logical indicating whether to scale data before PCA.
#'   Default is FALSE.
#' @param seed Integer seed for reproducibility. Default is NULL (no seed set).
#' @param return_pca Logical indicating whether to return PCA results.
#'   Default is TRUE.
#' @param store_in_colData Logical indicating whether to store cluster assignments
#'   in the colData of the input object (only applicable when x is an SCE/SPE).
#'   Default is FALSE.
#' @param cluster_name Character string specifying the column name for cluster
#'   assignments if store_in_colData is TRUE. Default is "PhiClust".
#'
#' @return A list with class "PhiSpaceClustering" containing:
#'   \item{clusters}{Factor vector of cluster assignments}
#'   \item{cluster_centers}{Matrix of cluster centers in PC space}
#'   \item{kmeans_result}{Full kmeans object from stats::kmeans}
#'   \item{pca_result}{PCA results (if return_pca = TRUE)}
#'   \item{pc_scores}{Matrix of PC scores used for clustering}
#'   \item{optimal_k}{Selected k value (relevant when k_range is used)}
#'   \item{k_selection}{List with selection metrics (if k_range was used)}
#'   \item{parameters}{List of parameters used}
#'   \item{spe}{Updated object (if store_in_colData = TRUE and x is SCE/SPE)}
#'
#' @details
#' This function implements a common workflow for spatial clustering:
#' 1. Extract data (from reducedDim, assay, or use directly if matrix)
#' 2. Perform PCA to reduce dimensionality
#' 3. Select top principal components
#' 4. Apply k-means clustering
#'
#' The function can either:
#' - Use a fixed k value (specify k parameter)
#' - Automatically select k from a range (specify k_range parameter)
#'
#' When k_range is provided, the function tests all k values in the range and
#' selects the optimal k using either:
#' - **Silhouette method** (default): Maximizes average silhouette width
#' - **Elbow method**: Identifies elbow point in total within-cluster sum of squares
#'
#' @section Input Types:
#' The function accepts multiple input types:
#' - **SingleCellExperiment/SpatialExperiment**: Uses reducedDim (default) or assay
#' - **Matrix**: Standard R matrix with observations in rows
#' - **Sparse matrix**: dgCMatrix or similar sparse formats
#' - **Data frame**: Coerced to matrix
#'
#' @section PCA Parameters:
#' The number of PCs to use affects clustering resolution:
#' - Fewer PCs (e.g., 10-15): Broader, more general clusters
#' - More PCs (e.g., 30-50): Finer, more specific clusters
#' - Default (30): Good balance for most applications
#'
#' @section K-means Algorithm:
#' - **Lloyd** (default): Standard algorithm, good balance of speed and quality
#' - **Hartigan-Wong**: Often better results but slower
#' - **MacQueen**: Faster but may converge to local optima
#' - **Forgy**: Similar to Lloyd
#'
#' @examples
#' \dontrun{
#' # Example 1: Using SingleCellExperiment with PhiSpace scores
#' result <- clusterPhiSpace(
#'   x = lung_data,
#'   k = 9,
#'   ncomp = 30,
#'   seed = 123
#' )
#'
#' # Example 2: Using matrix directly
#' phi_matrix <- reducedDim(lung_data, "PhiSpace")
#' result <- clusterPhiSpace(
#'   x = phi_matrix,
#'   k = 9,
#'   ncomp = 30
#' )
#'
#' # Example 3: Cluster on normalized counts instead
#' result <- clusterPhiSpace(
#'   x = lung_data,
#'   use_assay = "logcounts",
#'   k = 9,
#'   ncomp = 50
#' )
#'
#' # Example 4: Using data frame
#' df <- as.data.frame(reducedDim(lung_data, "PhiSpace"))
#' result <- clusterPhiSpace(x = df, k = 9)
#'
#' # Example 5: Automatic k selection
#' result <- clusterPhiSpace(
#'   x = phi_matrix,
#'   k_range = c(5, 15),
#'   select_k_method = "silhouette"
#' )
#' }
#'
#' @importFrom stats kmeans
#' @export
clusterPhiSpace <- function(
    x,
    k = NULL,
    k_range = NULL,
    select_k_method = c("silhouette", "elbow"),
    ncomp = NULL,
    reducedDimName = "PhiSpace",
    use_assay = NULL,
    nstart = 20,
    iter.max = 500,
    algorithm = c("Lloyd", "Hartigan-Wong", "MacQueen", "Forgy"),
    center = TRUE,
    scale = FALSE,
    seed = NULL,
    return_pca = TRUE,
    store_in_colData = FALSE,
    cluster_name = "PhiClust"
) {

  # Match arguments
  algorithm <- match.arg(algorithm)
  select_k_method <- match.arg(select_k_method)

  # Validate k and k_range
  if (is.null(k) && is.null(k_range)) {
    stop("Either 'k' or 'k_range' must be specified")
  }
  if (!is.null(k) && !is.null(k_range)) {
    stop("Specify either 'k' or 'k_range', not both")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine input type and extract data
  is_sce <- inherits(x, "SingleCellExperiment") || inherits(x, "SpatialExperiment")

  if (is_sce) {
    # Input is SingleCellExperiment or SpatialExperiment
    if (!is.null(use_assay)) {
      # Extract from assay
      if (!use_assay %in% assayNames(x)) {
        stop("Assay '", use_assay, "' not found in the object. ",
             "Available assays: ", paste(assayNames(x), collapse = ", "))
      }
      data_matrix <- t(as.matrix(assay(x, use_assay)))
      data_source <- paste0("assay:", use_assay)
    } else {
      # Extract from reducedDim
      if (!reducedDimName %in% reducedDimNames(x)) {
        stop("reducedDim '", reducedDimName, "' not found in the object. ",
             "Available reducedDims: ", paste(reducedDimNames(x), collapse = ", "))
      }
      data_matrix <- reducedDim(x, reducedDimName)
      data_source <- paste0("reducedDim:", reducedDimName)
    }
    original_object <- x
  } else {
    # Input is matrix-like object
    if (is.data.frame(x)) {
      data_matrix <- as.matrix(x)
      data_source <- "data.frame"
    } else if (is.matrix(x)) {
      data_matrix <- x
      data_source <- "matrix"
    } else if (inherits(x, "sparseMatrix") || inherits(x, "dgCMatrix")) {
      data_matrix <- as.matrix(x)
      data_source <- "sparse matrix"
    } else {
      # Try to coerce to matrix
      tryCatch({
        data_matrix <- as.matrix(x)
        data_source <- "coerced to matrix"
      }, error = function(e) {
        stop("Input 'x' must be a SingleCellExperiment, SpatialExperiment, ",
             "matrix, sparse matrix, or data.frame. Cannot coerce object of class '",
             paste(class(x), collapse = ", "), "' to matrix.")
      })
    }
    original_object <- NULL

    # Validate matrix dimensions
    if (nrow(data_matrix) < 2) {
      stop("Input matrix must have at least 2 observations (rows)")
    }
    if (ncol(data_matrix) < 2) {
      stop("Input matrix must have at least 2 features (columns)")
    }
  }

  # Ensure matrix has row names
  if (is.null(rownames(data_matrix))) {
    rownames(data_matrix) <- paste0("obs_", seq_len(nrow(data_matrix)))
  }

  # Determine number of PCs
  max_ncomp <- ncol(data_matrix) - 1
  if (is.null(ncomp)) {
    ncomp <- min(30, max_ncomp)
  } else if (is.character(ncomp) && ncomp == "all") {
    ncomp <- max_ncomp
  } else if (ncomp > max_ncomp) {
    warning("ncomp (", ncomp, ") exceeds maximum possible (", max_ncomp,
            "). Using ", max_ncomp, " components.")
    ncomp <- max_ncomp
  }

  # Perform PCA
  pca_result <- getPC(
    data_matrix,
    ncomp = max_ncomp,
    center = center,
    scale = scale
  )

  # Extract PC scores for clustering
  pc_scores <- pca_result$scores[, 1:ncomp, drop = FALSE]

  # Determine k value(s)
  if (!is.null(k_range)) {
    # Automatic k selection
    k_select_result <- .select_optimal_k(
      pc_scores,
      k_range,
      select_k_method,
      nstart = nstart,
      iter.max = iter.max,
      algorithm = algorithm
    )
    k <- k_select_result$optimal_k
    k_selection <- k_select_result
  } else {
    k_selection <- NULL
  }

  # Perform k-means clustering
  kmeans_result <- kmeans(
    x = pc_scores,
    centers = k,
    iter.max = iter.max,
    nstart = nstart,
    algorithm = algorithm
  )

  # Create cluster factor
  clusters <- factor(kmeans_result$cluster)
  names(clusters) <- rownames(pc_scores)

  # Store in colData if requested and possible
  if (store_in_colData) {
    if (!is_sce) {
      warning("store_in_colData = TRUE but input is not a SingleCellExperiment/SpatialExperiment. ",
              "Cluster assignments will not be stored in the object.")
    } else {
      colData(original_object)[[cluster_name]] <- clusters
    }
  }

  # Prepare output
  result <- list(
    clusters = clusters,
    cluster_centers = kmeans_result$centers,
    kmeans_result = kmeans_result,
    pca_result = if (return_pca) pca_result else NULL,
    pc_scores = pc_scores,
    optimal_k = k,
    k_selection = k_selection,
    parameters = list(
      k = k,
      k_range = k_range,
      select_k_method = if (!is.null(k_range)) select_k_method else NA,
      ncomp = ncomp,
      nstart = nstart,
      iter.max = iter.max,
      algorithm = algorithm,
      center = center,
      scale = scale,
      data_source = data_source,
      input_type = if (is_sce) class(x)[1] else class(x)[1]
    ),
    spe = if (store_in_colData && is_sce) original_object else NULL
  )

  class(result) <- c("PhiSpaceClustering", "list")
  return(result)
}


#' Internal function to select optimal k
#' @keywords internal
.select_optimal_k <- function(pc_scores, k_range, method, nstart, iter.max, algorithm) {

  k_values <- seq(k_range[1], k_range[2])

  if (method == "silhouette") {
    # Silhouette method
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package 'cluster' is required for silhouette method. ",
           "Please install it with: install.packages('cluster')")
    }

    sil_widths <- sapply(k_values, function(k) {
      km <- kmeans(pc_scores, centers = k, nstart = nstart,
                   iter.max = iter.max, algorithm = algorithm)
      sil <- cluster::silhouette(km$cluster, dist(pc_scores))
      mean(sil[, 3])
    })

    optimal_k <- k_values[which.max(sil_widths)]

    return(list(
      optimal_k = optimal_k,
      method = "silhouette",
      k_values = k_values,
      silhouette_widths = sil_widths
    ))

  } else {  # elbow method

    wss <- sapply(k_values, function(k) {
      km <- kmeans(pc_scores, centers = k, nstart = nstart,
                   iter.max = iter.max, algorithm = algorithm)
      km$tot.withinss
    })

    # Find elbow using second derivative
    if (length(k_values) >= 3) {
      second_deriv <- diff(diff(wss))
      # Elbow is where second derivative is maximum (most negative)
      elbow_idx <- which.max(second_deriv) + 1
      optimal_k <- k_values[elbow_idx]
    } else {
      optimal_k <- k_values[1]
    }

    return(list(
      optimal_k = optimal_k,
      method = "elbow",
      k_values = k_values,
      wss = wss
    ))
  }
}


#' Print method for PhiSpaceClustering objects
#' @param x A PhiSpaceClustering object
#' @param ... Additional arguments (not used)
#' @method print PhiSpaceClustering
#' @export
print.PhiSpaceClustering <- function(x, ...) {
  cat("PhiSpace K-means Clustering\n")
  cat("===========================\n\n")

  cat("Input type:", x$parameters$input_type, "\n")
  cat("Data source:", x$parameters$data_source, "\n")
  cat("Number of clusters (k):", x$optimal_k, "\n")

  if (!is.null(x$k_selection)) {
    cat("Selection method:", x$k_selection$method, "\n")
    cat("k range tested:", paste(range(x$k_selection$k_values), collapse = " to "), "\n")
  }

  cat("Number of PCs used:", ncol(x$pc_scores), "\n")
  cat("Algorithm:", x$parameters$algorithm, "\n")
  cat("Number of starts:", x$parameters$nstart, "\n")

  cat("\nCluster sizes:\n")
  print(table(x$clusters))

  cat("\nCluster quality:\n")
  cat("  Total within-cluster SS:",
      round(x$kmeans_result$tot.withinss, 2), "\n")
  cat("  Between-cluster SS:",
      round(x$kmeans_result$betweenss, 2), "\n")
  cat("  Total SS:",
      round(x$kmeans_result$totss, 2), "\n")
  cat("  Between/Total ratio:",
      round(x$kmeans_result$betweenss / x$kmeans_result$totss, 3), "\n")

  if (!is.null(x$pca_result)) {
    var_explained <- sum(x$pca_result$var.exp[1:ncol(x$pc_scores)])
    cat("\nVariance explained by selected PCs:",
        round(var_explained * 100, 1), "%\n")
  }

  cat("\nUse plot() to visualize results\n")
  cat("Use summary() for more details\n")
}


#' Summary method for PhiSpaceClustering objects
#' @param object A PhiSpaceClustering object
#' @param ... Additional arguments (not used)
#' @method summary PhiSpaceClustering
#' @export
summary.PhiSpaceClustering <- function(object, ...) {
  cat("PhiSpace K-means Clustering - Detailed Summary\n")
  cat("==============================================\n\n")

  print(object)

  if (!is.null(object$k_selection)) {
    cat("\n\nK selection details:\n")
    if (object$k_selection$method == "silhouette") {
      cat("Silhouette widths for each k:\n")
      sil_df <- data.frame(
        k = object$k_selection$k_values,
        avg_silhouette = round(object$k_selection$silhouette_widths, 4)
      )
      print(sil_df)
    } else {
      cat("Within-cluster SS for each k:\n")
      wss_df <- data.frame(
        k = object$k_selection$k_values,
        total_withinss = round(object$k_selection$wss, 2)
      )
      print(wss_df)
    }
  }

  cat("\n\nWithin-cluster SS by cluster:\n")
  print(round(object$kmeans_result$withinss, 2))

  if (!is.null(object$pca_result)) {
    cat("\n\nPCA variance explained:\n")
    n_show <- min(10, length(object$pca_result$var.exp))
    pca_df <- data.frame(
      PC = paste0("PC", 1:n_show),
      Variance = round(object$pca_result$var.exp[1:n_show], 4),
      Cumulative = round(cumsum(object$pca_result$var.exp)[1:n_show], 4)
    )
    print(pca_df)
  }

  invisible(object)
}


#' Plot method for PhiSpaceClustering objects
#' @param x A PhiSpaceClustering object
#' @param type Character string specifying plot type. Options: "pca" (PC1 vs PC2),
#'   "elbow" (elbow plot if k_range was used), "silhouette" (silhouette plot if
#'   k_range was used), "variance" (PCA variance explained). Default is "pca".
#' @param ... Additional arguments passed to plotting functions
#' @method plot PhiSpaceClustering
#' @export
plot.PhiSpaceClustering <- function(x, type = c("pca", "elbow", "silhouette", "variance"), ...) {

  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  if (type == "pca") {
    # PC1 vs PC2 colored by cluster
    plot_df <- data.frame(
      PC1 = x$pc_scores[, 1],
      PC2 = x$pc_scores[, 2],
      Cluster = x$clusters
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = Cluster)) +
      ggplot2::geom_point(alpha = 0.6, size = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "K-means Clustering on PhiSpace PCs",
        subtitle = paste0("k = ", x$optimal_k, ", ", ncol(x$pc_scores), " PCs used"),
        x = paste0("PC1 (", round(x$pca_result$var.exp[1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(x$pca_result$var.exp[2] * 100, 1), "%)")
      )

    return(p)

  } else if (type == "elbow") {

    if (is.null(x$k_selection) || x$k_selection$method != "elbow") {
      stop("Elbow plot only available when k_range was used with method='elbow'")
    }

    plot_df <- data.frame(
      k = x$k_selection$k_values,
      wss = x$k_selection$wss
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = wss)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_point(data = subset(plot_df, k == x$optimal_k),
                          color = "red", size = 5) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Elbow Method for Optimal k",
        subtitle = paste("Selected k =", x$optimal_k),
        x = "Number of Clusters (k)",
        y = "Total Within-Cluster Sum of Squares"
      )

    return(p)

  } else if (type == "silhouette") {

    if (is.null(x$k_selection) || x$k_selection$method != "silhouette") {
      stop("Silhouette plot only available when k_range was used with method='silhouette'")
    }

    plot_df <- data.frame(
      k = x$k_selection$k_values,
      silhouette = x$k_selection$silhouette_widths
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = silhouette)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_point(data = subset(plot_df, k == x$optimal_k),
                          color = "red", size = 5) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Silhouette Method for Optimal k",
        subtitle = paste("Selected k =", x$optimal_k),
        x = "Number of Clusters (k)",
        y = "Average Silhouette Width"
      )

    return(p)

  } else {  # variance

    if (is.null(x$pca_result)) {
      stop("PCA results not available. Set return_pca=TRUE when calling clusterPhiSpace()")
    }

    n_show <- min(30, length(x$pca_result$var.exp))
    plot_df <- data.frame(
      PC = 1:n_show,
      Variance = x$pca_result$var.exp[1:n_show]
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC, y = Variance)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::geom_vline(xintercept = ncol(x$pc_scores),
                          linetype = "dashed", color = "red") +
      ggplot2::annotate("text", x = ncol(x$pc_scores), y = max(plot_df$Variance),
                        label = paste("Used:", ncol(x$pc_scores), "PCs"),
                        hjust = -0.1, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "PCA Variance Explained",
        subtitle = paste("First", n_show, "principal components"),
        x = "Principal Component",
        y = "Proportion of Variance Explained"
      )

    return(p)
  }
}
