#' Compute UMAP.
#'
#' @param dat Input data.
#' @param computPC Whether to compute PC first.
#' @param ncomp Number of PC components
#' @param center Center or not.
#' @param scale Scale or not.
#' @param config Configuration of UMAP (see `umap` package).
#'
#' @return UMAP coordinates.
#' @export
computeUMAP <- function(
    dat,
    computPC = TRUE,
    ncomp = 30,
    center = TRUE,
    scale = FALSE,
    config = NULL
){

  if(is.null(rownames(dat))) warning("Rownames (sample names) of input data is null.")

  if(computPC){

    dat <- getPC(
      dat,
      ncomp = ncomp,
      center = center,
      scale = scale
    )$scores
  }

  if(is.null(config)) config <- umap::umap.defaults

  out <- umap::umap(dat, config = config)$layout %>% reNameCols()
  rownames(out) <- rownames(dat)

  return(out)
}
