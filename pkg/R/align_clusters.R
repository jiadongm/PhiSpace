# Older version before 14 Aug 2025
# align_clusters <- function(clust, clust_ref){
#
#   clust <- as.factor(clust)
#   clust_ref <- as.factor(clust_ref)
#   lvls <- levels(clust)
#   lvls_old <- levels(clust_ref)
#
#   sol <- clue::solve_LSAP(
#     table(clust, clust_ref), T
#   )
#   kclust <- length(levels(clust_ref))
#   adj <- (1:kclust)[sol]
#   temp <- rep(NA, length(clust))
#   for(jj in 1:length(adj)){
#
#     temp[clust == lvls[jj]] <- lvls_old[adj[jj]]
#   }
#   clust <- factor(temp, levels = lvls_old)
#
#   return(clust)
# }


#' Align clustering results to a reference.
#'
#' Given a reference clustering result, align one or more query clustering results
#' to the reference using the Hungarian algorithm. Can handle both single clustering
#' vectors and multiple clusterings in a matrix/data.frame format.
#'
#' @param clust Query clustering result(s). Can be:
#'   - A vector (numeric, character, or factor) for single clustering
#'   - A matrix/data.frame where each column is a clustering result
#'   - A list of clustering vectors
#' @param clust_ref Reference clustering result (vector).
#' @param clust_names Optional names for the clustering results (for matrix/list input).
#' @param handle_unmatched How to handle unmatched clusters:
#'   - "label": create new labels like "Unmatched_X" (default)
#'   - "NA": assign NA to unmatched clusters
#'   - "nearest": assign to nearest matching cluster
#' @param verbose Logical. Print alignment progress for multiple clusterings.
#'
#' @return
#'   - If input is a single vector: aligned clustering (factor)
#'   - If input is matrix/data.frame: matrix with aligned clusterings
#'   - If input is a list: list with aligned clusterings
#' @export
#'
#' @examples
#' # Single clustering alignment
#' ref <- factor(c("A", "A", "B", "B", "C", "C"))
#' query <- factor(c(2, 2, 1, 1, 3, 3))
#' aligned <- align_clusters(query, ref)
#'
#' # Multiple clusterings alignment
#' query_matrix <- cbind(
#'   method1 = c(1, 1, 2, 2, 3, 3),
#'   method2 = c("a", "a", "b", "b", "c", "c"),
#'   method3 = c(2, 2, 2, 1, 1, 1)
#' )
#' aligned_all <- align_clusters(query_matrix, ref)
#'
#' # List of clusterings
#' query_list <- list(
#'   kmeans = c(1, 1, 2, 2, 3, 3),
#'   hierarchical = c("x", "x", "y", "y", "z", "z")
#' )
#' aligned_list <- align_clusters(query_list, ref)
#'
align_clusters <- function(clust,
                           clust_ref,
                           clust_names = NULL,
                           handle_unmatched = c("label", "NA", "nearest"),
                           verbose = FALSE) {

  # Input validation
  handle_unmatched <- match.arg(handle_unmatched)

  if (length(clust_ref) == 0) {
    stop("Reference clustering cannot be empty.")
  }

  clust_ref <- as.factor(clust_ref)

  # Determine input type and process accordingly
  if (is.list(clust) && !is.data.frame(clust)) {
    # List of clusterings
    return(align_clusters_list(clust, clust_ref, clust_names,
                               handle_unmatched, verbose))

  } else if (is.matrix(clust) || is.data.frame(clust)) {
    # Matrix or data.frame of clusterings
    return(align_clusters_matrix(clust, clust_ref, clust_names,
                                 handle_unmatched, verbose))

  } else {
    # Single clustering vector
    return(align_single_clustering(clust, clust_ref, handle_unmatched))
  }
}


#' Align a single clustering to reference
#'
#' Internal function for aligning a single clustering vector
#'
#' @keywords internal
align_single_clustering <- function(clust, clust_ref, handle_unmatched = "label") {

  # Validate inputs
  if (length(clust) == 0) {
    stop("Query clustering cannot be empty.")
  }

  if (length(clust) != length(clust_ref)) {
    warning("Query and reference clusterings have different lengths. ",
            "Make sure they correspond to the same samples.")
  }

  # Convert to factors
  clust <- as.factor(clust)
  clust_ref <- as.factor(clust_ref)

  # Get unique labels
  lvls_query <- levels(clust)
  lvls_ref <- levels(clust_ref)

  # Create contingency table
  cont_table <- table(clust, clust_ref)

  # Handle unequal cluster numbers
  n_query <- length(lvls_query)
  n_ref <- length(lvls_ref)

  if (n_query != n_ref) {
    # Pad the contingency table to make it square
    max_dim <- max(n_query, n_ref)
    padded_table <- matrix(0, nrow = max_dim, ncol = max_dim)
    padded_table[1:n_query, 1:n_ref] <- cont_table
    cont_table <- padded_table
  }

  # Solve the assignment problem (maximize overlap)
  sol <- clue::solve_LSAP(cont_table, maximum = TRUE)

  # Create mapping from query to reference labels
  mapping <- rep(NA_character_, n_query)
  for (i in 1:n_query) {
    ref_idx <- sol[i]
    if (ref_idx <= n_ref) {
      mapping[i] <- lvls_ref[ref_idx]
    } else {
      # Handle unmatched clusters
      if (handle_unmatched == "label") {
        mapping[i] <- paste0("Unmatched_", lvls_query[i])
      } else if (handle_unmatched == "NA") {
        mapping[i] <- NA_character_
      } else if (handle_unmatched == "nearest") {
        # Assign to cluster with highest overlap
        overlaps <- cont_table[i, 1:n_ref]
        if (sum(overlaps) > 0) {
          mapping[i] <- lvls_ref[which.max(overlaps)]
        } else {
          mapping[i] <- lvls_ref[1]  # Default to first cluster
        }
      }
    }
  }

  # Apply the mapping
  if (handle_unmatched == "label") {
    all_levels <- c(lvls_ref,
                    paste0("Unmatched_", lvls_query[!lvls_query %in% mapping]))
  } else {
    all_levels <- lvls_ref
  }

  aligned <- factor(mapping[match(clust, lvls_query)], levels = all_levels)

  # Remove unused levels
  aligned <- droplevels(aligned)

  return(aligned)
}


#' Align multiple clusterings in matrix format
#'
#' @keywords internal
align_clusters_matrix <- function(clust_matrix,
                                  clust_ref,
                                  clust_names = NULL,
                                  handle_unmatched = "label",
                                  verbose = FALSE) {

  # Validate dimensions
  if (nrow(clust_matrix) != length(clust_ref)) {
    warning("Number of rows in clustering matrix doesn't match reference length. ",
            "Make sure rows correspond to the same samples.")
  }

  # Set column names if provided
  if (!is.null(clust_names)) {
    if (length(clust_names) != ncol(clust_matrix)) {
      stop("Length of clust_names must match number of columns in clustering matrix.")
    }
    colnames(clust_matrix) <- clust_names
  }

  # Initialize result matrix
  n_clusterings <- ncol(clust_matrix)
  aligned_matrix <- matrix(NA,
                           nrow = nrow(clust_matrix),
                           ncol = n_clusterings)
  colnames(aligned_matrix) <- colnames(clust_matrix)
  rownames(aligned_matrix) <- rownames(clust_matrix)

  # Align each clustering
  for (i in 1:n_clusterings) {
    if (verbose) {
      col_name <- colnames(clust_matrix)[i]
      if (is.null(col_name)) col_name <- paste0("Clustering_", i)
      cat(sprintf("Aligning clustering %d/%d: %s\n", i, n_clusterings, col_name))
    }

    aligned_clust <- align_single_clustering(clust_matrix[, i],
                                             clust_ref,
                                             handle_unmatched)
    aligned_matrix[, i] <- as.character(aligned_clust)
  }

  # Convert to data.frame with factors
  aligned_df <- as.data.frame(aligned_matrix, stringsAsFactors = TRUE)
  colnames(aligned_df) <- colnames(clust_matrix)

  return(aligned_df)
}


#' Align multiple clusterings in list format
#'
#' @keywords internal
align_clusters_list <- function(clust_list,
                                clust_ref,
                                clust_names = NULL,
                                handle_unmatched = "label",
                                verbose = FALSE) {

  # Validate list
  if (!all(sapply(clust_list, function(x) is.vector(x) || is.factor(x)))) {
    stop("All elements in the list must be vectors or factors.")
  }

  # Check lengths
  lengths <- sapply(clust_list, length)
  if (!all(lengths == length(clust_ref))) {
    warning("Not all clusterings have the same length as the reference. ",
            "Make sure they correspond to the same samples.")
  }

  # Set names if provided
  if (!is.null(clust_names)) {
    if (length(clust_names) != length(clust_list)) {
      stop("Length of clust_names must match length of clustering list.")
    }
    names(clust_list) <- clust_names
  }

  # Align each clustering
  aligned_list <- list()
  n_clusterings <- length(clust_list)

  for (i in 1:n_clusterings) {
    name <- names(clust_list)[i]
    if (is.null(name)) name <- paste0("Clustering_", i)

    if (verbose) {
      cat(sprintf("Aligning clustering %d/%d: %s\n", i, n_clusterings, name))
    }

    aligned_list[[name]] <- align_single_clustering(clust_list[[i]],
                                                    clust_ref,
                                                    handle_unmatched)
  }

  return(aligned_list)
}
