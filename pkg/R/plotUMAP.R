#' Compute UMAP.
#'
#' @param dat Input data.
#' @param computPC Whether to compute PC first.
#' @param ncomp Number of PC components
#' @param center Center or not.
#' @param scale Scale or not.
#'
#' @return UMAP coordinates.
#' @export
computUMAP <- function(
    dat,
    computPC = TRUE,
    ncomp = 30,
    center = TRUE,
    scale = FALSE
){

  if(computPC){

    plot_dat <- getPC(
      dat,
      ncomp = ncomp,
      center = center,
      scale = scale
    )$scores
  } else {

    plot_dat <- dat
  }

  umap::umap(
    query_pc$scores
  )$layout %>%
    reNameCols()
}
