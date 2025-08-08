# spatialSampler: Spatially-aware subsampling for SpatialExperiment objects
#
# Description: Subsample cells from SpatialExperiment objects while maintaining
#              even spatial distribution using grid, k-means, or random methods

#' Spatially-aware subsampling of SpatialExperiment objects
#'
#' @param spe SpatialExperiment object
#' @param prop Proportion of cells to sample (between 0 and 1)
#' @param method Sampling method: "grid", "kmeans", or "random" (default: "grid")
#' @param grid_size Number of grid cells per dimension for grid method (default: auto-calculated)
#' @param seed Random seed for reproducibility (default: 123)
#' @param min_cells_per_region Minimum cells required per spatial region (default: 1)
#' @param balance_regions Whether to balance sampling across spatial regions (default: TRUE)
#'
#' @return Subsampled SpatialExperiment object
#'
#' @export
spatialSampler <- function(spe,
                           prop = 0.1,
                           method = "grid",
                           grid_size = NULL,
                           seed = 123,
                           min_cells_per_region = 1,
                           balance_regions = TRUE) {

  # Input validation
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Object must be a SpatialExperiment")
  }

  if (prop <= 0 || prop > 1) {
    stop("prop must be between 0 and 1")
  }

  if (!method %in% c("grid", "kmeans", "random")) {
    stop("method must be 'grid', 'kmeans', or 'random'")
  }

  set.seed(seed)

  # Get spatial coordinates
  coords_matrix <- spatialCoords(spe)
  if (ncol(coords_matrix) < 2) {
    stop("SpatialExperiment object must have at least 2 spatial coordinates")
  }

  coordinates <- data.frame(
    x = coords_matrix[, 1],
    y = coords_matrix[, 2],
    cell_id = seq_len(ncol(spe))
  )

  n_total <- nrow(coordinates)
  n_target <- round(n_total * prop)

  cat("Spatial sampling:", n_target, "cells from", n_total, "total cells\n")
  cat("Method:", method, "\n")

  # Perform spatial sampling based on method
  if (method == "random") {
    selected_indices <- sample_random(coordinates, n_target)

  } else if (method == "grid") {
    selected_indices <- sample_grid(coordinates, n_target, grid_size,
                                    min_cells_per_region, balance_regions)

  } else if (method == "kmeans") {
    selected_indices <- sample_kmeans(coordinates, n_target,
                                      min_cells_per_region, balance_regions)
  }

  cat("Actually sampled:", length(selected_indices), "cells\n")

  # Subset the SpatialExperiment object
  spe_subset <- spe[, selected_indices]

  # Add metadata about sampling
  metadata_entry <- list(
    function_used = "spatialSampler",
    original_n_cells = n_total,
    sampled_n_cells = length(selected_indices),
    target_prop = prop,
    actual_prop = length(selected_indices) / n_total,
    method = method,
    seed = seed,
    timestamp = Sys.time()
  )

  if (method == "grid") {
    metadata_entry$grid_size <- if (is.null(grid_size)) "auto" else grid_size
  }

  S4Vectors::metadata(spe_subset)$spatialSampling <- c(S4Vectors::metadata(spe_subset)$spatialSampling, list(metadata_entry))

  return(spe_subset)
}

#' Random sampling (baseline method)
#'
#' @param coordinates Spatial coordiantes
#' @param n_target Number of targets
sample_random <- function(coordinates, n_target) {
  return(sample(coordinates$cell_id, n_target))
}



#' Grid-based spatial sampling
#'
#' @param coordinates Spatial coordinates
#' @param n_target Number of neighbours
#' @param grid_size Grid size
#' @param min_cells_per_region Mimimum cells per region
#' @param balance_regions Balance regions
sample_grid <- function(coordinates, n_target, grid_size, min_cells_per_region, balance_regions) {

  # Auto-calculate grid size if not provided
  if (is.null(grid_size)) {
    # Aim for roughly sqrt(n_target) regions
    grid_size <- max(2, round(sqrt(n_target)))
  }

  cat("Using grid size:", grid_size, "x", grid_size, "\n")

  # Create spatial grid
  x_breaks <- seq(min(coordinates$x), max(coordinates$x), length.out = grid_size + 1)
  y_breaks <- seq(min(coordinates$y), max(coordinates$y), length.out = grid_size + 1)

  # Assign cells to grid cells
  coordinates$x_bin <- cut(coordinates$x, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  coordinates$y_bin <- cut(coordinates$y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  coordinates$grid_cell <- paste0(coordinates$x_bin, "_", coordinates$y_bin)

  # Count cells per grid cell
  grid_counts <- table(coordinates$grid_cell)
  occupied_grids <- names(grid_counts)[grid_counts >= min_cells_per_region]

  cat("Grid cells with sufficient cells:", length(occupied_grids), "out of", grid_size^2, "\n")

  if (balance_regions) {
    # Sample equally from each occupied grid cell
    cells_per_grid <- floor(n_target / length(occupied_grids))
    remainder <- n_target - (cells_per_grid * length(occupied_grids))

    selected_indices <- c()

    for (i in seq_along(occupied_grids)) {
      grid_id <- occupied_grids[i]
      grid_cells <- coordinates$cell_id[coordinates$grid_cell == grid_id]

      # Add one extra cell to some grids to handle remainder
      n_sample <- cells_per_grid + (if (i <= remainder) 1 else 0)
      n_sample <- min(n_sample, length(grid_cells))

      if (n_sample > 0) {
        selected_indices <- c(selected_indices, sample(grid_cells, n_sample))
      }
    }

  } else {
    # Sample proportionally from each grid cell
    selected_indices <- c()

    for (grid_id in occupied_grids) {
      grid_cells <- coordinates$cell_id[coordinates$grid_cell == grid_id]
      grid_prop <- length(grid_cells) / sum(grid_counts[occupied_grids])
      n_sample <- round(n_target * grid_prop)
      n_sample <- min(n_sample, length(grid_cells))

      if (n_sample > 0) {
        selected_indices <- c(selected_indices, sample(grid_cells, n_sample))
      }
    }
  }

  return(selected_indices)
}

#' K-means based spatial sampling
#'
#' @param coordinates Spatial coordinates
#' @param n_target Number of targets
#' @param min_cells_per_region Minimum cells per region
#' @param balance_regions Balance regions
sample_kmeans <- function(coordinates, n_target, min_cells_per_region, balance_regions) {

  # Determine number of clusters (spatial regions)
  n_clusters <- max(2, min(round(sqrt(n_target)), round(nrow(coordinates) / min_cells_per_region)))

  cat("Using", n_clusters, "spatial clusters\n")

  # Perform k-means clustering on coordinates
  coord_matrix <- as.matrix(coordinates[, c("x", "y")])

  # Use robust k-means (avoid convergence warnings)
  kmeans_result <- tryCatch({
    stats::kmeans(scale(coord_matrix), centers = n_clusters, algorithm = "Lloyd",
           iter.max = 100, nstart = 25)
  }, error = function(e) {
    # Fallback to simpler clustering
    stats::kmeans(coord_matrix, centers = n_clusters, algorithm = "Lloyd", nstart = 10)
  })

  coordinates$spatial_cluster <- kmeans_result$cluster

  # Count cells per cluster
  cluster_counts <- table(coordinates$spatial_cluster)
  valid_clusters <- names(cluster_counts)[cluster_counts >= min_cells_per_region]

  cat("Valid spatial clusters:", length(valid_clusters), "out of", n_clusters, "\n")

  if (balance_regions) {
    # Sample equally from each cluster
    cells_per_cluster <- floor(n_target / length(valid_clusters))
    remainder <- n_target - (cells_per_cluster * length(valid_clusters))

    selected_indices <- c()

    for (i in seq_along(valid_clusters)) {
      cluster_id <- valid_clusters[i]
      cluster_cells <- coordinates$cell_id[coordinates$spatial_cluster == cluster_id]

      # Add one extra cell to some clusters to handle remainder
      n_sample <- cells_per_cluster + (if (i <= remainder) 1 else 0)
      n_sample <- min(n_sample, length(cluster_cells))

      if (n_sample > 0) {
        selected_indices <- c(selected_indices, sample(cluster_cells, n_sample))
      }
    }

  } else {
    # Sample proportionally from each cluster
    selected_indices <- c()

    for (cluster_id in valid_clusters) {
      cluster_cells <- coordinates$cell_id[coordinates$spatial_cluster == cluster_id]
      cluster_prop <- length(cluster_cells) / sum(cluster_counts[valid_clusters])
      n_sample <- round(n_target * cluster_prop)
      n_sample <- min(n_sample, length(cluster_cells))

      if (n_sample > 0) {
        selected_indices <- c(selected_indices, sample(cluster_cells, n_sample))
      }
    }
  }

  return(selected_indices)
}


