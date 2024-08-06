#' Align clustering results.
#'
#' Given a reference clustering result, represented by a vector of cluster labels, align a query clustering result to the reference one.
#'
#' @param clust Query clustering result, can be a numerical, character or factor vector.
#' @param clust_ref  Reference clustering result.
#'
#' @return Re-aligned clustering result, represented as a factor vector.
#' @export
#'
align_clusters <- function(clust, clust_ref){

  clust <- as.factor(clust)
  clust_ref <- as.factor(clust_ref)
  lvls <- levels(clust)
  lvls_old <- levels(clust_ref)

  sol <- clue::solve_LSAP(
    table(clust, clust_ref), T
  )
  kclust <- length(levels(clust_ref))
  adj <- (1:kclust)[sol]
  temp <- rep(NA, length(clust))
  for(jj in 1:length(adj)){

    temp[clust == lvls[jj]] <- lvls_old[adj[jj]]
  }
  clust <- factor(temp, levels = lvls_old)

  return(clust)
}
