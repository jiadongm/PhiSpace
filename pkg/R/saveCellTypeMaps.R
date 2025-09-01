#' Save spatial cell type annotation heatmaps
#'
#' Creates and optionally saves spatial heatmaps for all cell types from PhiSpace
#' continuous annotations. This function generates publication-ready visualizations
#' showing the spatial distribution of cell type scores across tissue sections.
#'
#' @param sce A `SingleCellExperiment` or `SpatialExperiment` object containing
#'   spatial coordinates and PhiSpace annotations in `reducedDim` slot.
#' @param methodName Character. Name of the reduced dimension containing cell type
#'   scores (default: "PhiSpace").
#' @param tissueName Character. Name identifier for the tissue/sample, used in
#'   output filenames when saving plots (default: "SpatialSample").
#' @param coordNames Character vector of length 2. Column names in `colData(sce)`
#'   containing x and y spatial coordinates. For `SpatialExperiment` objects,
#'   this will be automatically detected if not specified (default: c("x", "y")).
#' @param freeColScale Logical. Whether to use independent color scales for each
#'   cell type (TRUE) or a unified scale across all cell types (FALSE). Unified
#'   scaling facilitates comparison between cell types (default: FALSE).
#' @param quants Numeric vector of length 2. Quantiles used to set color scale
#'   limits when `freeColScale = FALSE` (default: c(0.1, 1)).
#' @param psize Numeric. Point size for spatial coordinates (default: 0.5).
#' @param savePlots Logical. Whether to save plots to disk (default: TRUE).
#' @param returnPlots Logical. Whether to return the list of ggplot objects
#'   (default: FALSE). Note: renamed from `returnPlot` for clarity.
#' @param legendPosition Character. Position of color legend. One of "none",
#'   "left", "right", "bottom", "top" (default: "none").
#' @param censQuant Numeric between 0 and 1. Quantile threshold for censoring
#'   low values to reduce noise in visualization (default: 1, no censoring).
#' @param ctypes Character vector. Specific cell types to plot. If NULL, all
#'   cell types from the annotation will be plotted (default: NULL).
#' @param outputDir Character. Directory path where plots will be saved.
#'   Directory will be created if it doesn't exist (default: "figs/spatialCellTypeHeatmaps").
#' @param width Numeric. Plot width in inches when saving (default: 10).
#' @param height Numeric. Plot height in inches when saving (default: 10).
#' @param fignrow Integer. Number of rows in multi-panel figures (default: 4).
#' @param figncol Integer. Number of columns in multi-panel figures (default: 4).
#' @param fsize_title Numeric. Font size for plot titles (default: 10).
#' @param plot_margin Numeric vector of length 4. Plot margins in points
#'   (top, right, bottom, left) (default: c(0,0,0,0)).
#' @param fileFormat Character. File format for saved plots. One of "png", "pdf",
#'   "jpeg", "tiff" (default: "png").
#' @param dpi Numeric. Resolution for raster formats (default: 300).
#'
#' @return If `returnPlots = TRUE`, returns a named list of ggplot objects,
#'   one for each cell type. If `savePlots = TRUE`, plots are saved to the
#'   specified directory. Invisible NULL otherwise.
#'
#' @details
#' This function creates spatial heatmaps showing the distribution of PhiSpace
#' cell type scores across tissue sections. Each cell type is visualized as a
#' separate heatmap where color intensity represents the confidence score for
#' that cell type at each spatial location.
#'
#' **Key features:**
#' - Automatic handling of both SingleCellExperiment and SpatialExperiment objects
#' - Flexible color scaling options (unified vs. independent scales)
#' - Multi-panel output for efficient visualization of many cell types
#' - Customizable censoring to reduce visualization noise
#' - Multiple output formats supported
#'
#' **Color scaling:**
#' - `freeColScale = FALSE`: Uses the same color scale across all cell types,
#'   facilitating direct comparison of scores between different cell types
#' - `freeColScale = TRUE`: Each cell type uses its own optimal color scale,
#'   maximizing contrast within each individual plot
#'
#' **Output organization:**
#' When saving plots, multiple cell types are arranged in multi-panel figures
#' (e.g., 4x4 grids). If there are more cell types than can fit in one figure,
#' multiple files will be created with sequential numbering.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' saveCellTypeMaps(sce_with_phispace)
#'
#' # Customize visualization parameters
#' saveCellTypeMaps(
#'   sce = spatial_sce,
#'   tissueName = "MouseBrain_Section1",
#'   psize = 0.8,
#'   censQuant = 0.9,
#'   freeColScale = TRUE
#' )
#'
#' # Return plots for further customization
#' plots <- saveCellTypeMaps(
#'   sce = spatial_sce,
#'   savePlots = FALSE,
#'   returnPlots = TRUE
#' )
#'
#' # Plot specific cell types only
#' saveCellTypeMaps(
#'   sce = spatial_sce,
#'   ctypes = c("Neurons", "Astrocytes", "Microglia"),
#'   fignrow = 1, figncol = 3
#' )
#' }
#'
#' @seealso
#' \code{\link{VizSpatial}} for single cell type visualization
#' \code{\link{PhiSpace}} for generating cell type annotations
#'
#' @export
saveCellTypeMaps <- function(
    sce,
    methodName = "PhiSpace",
    tissueName = "SpatialSample",
    coordNames = c("x", "y"),
    freeColScale = FALSE,
    quants = c(0.1, 1),
    psize = 0.5,
    savePlots = TRUE,
    returnPlots = FALSE,  # renamed for clarity
    legendPosition = "none",
    censQuant = 1,
    ctypes = NULL,
    outputDir = "figs/spatialCellTypeHeatmaps",
    width = 10,
    height = 10,
    fignrow = 4,
    figncol = 4,
    fsize_title = 10,
    plot_margin = c(0, 0, 0, 0),
    fileFormat = "png",
    dpi = 300
) {

  # Input validation
  if (!inherits(sce, c("SingleCellExperiment", "SpatialExperiment"))) {
    stop("sce must be a SingleCellExperiment or SpatialExperiment object")
  }

  if (!methodName %in% reducedDimNames(sce)) {
    stop(paste("Method", methodName, "not found in reducedDims. Available:",
               paste(reducedDimNames(sce), collapse = ", ")))
  }

  if (!fileFormat %in% c("png", "pdf", "jpeg", "tiff")) {
    stop("fileFormat must be one of: png, pdf, jpeg, tiff")
  }

  if (length(quants) != 2 || quants[1] >= quants[2]) {
    stop("quants must be a vector of length 2 with quants[1] < quants[2]")
  }

  if (censQuant < 0 || censQuant > 1) {
    stop("censQuant must be between 0 and 1")
  }

  # Handle SpatialExperiment coordinate detection
  if (inherits(sce, "SpatialExperiment") && all(coordNames == c("x", "y"))) {
    spatial_coords <- spatialCoordsNames(sce)
    if (length(spatial_coords) >= 2) {
      coordNames <- spatial_coords[1:2]
      message("Using spatial coordinates: ", paste(coordNames, collapse = ", "))
    }
  }

  # Check coordinate columns
  missing_coords <- coordNames[!coordNames %in% colnames(colData(sce))]
  if (length(missing_coords) > 0) {
    stop("Coordinate columns not found in colData: ", paste(missing_coords, collapse = ", "))
  }

  # Extract scores matrix
  scores <- as.matrix(reducedDim(sce, methodName))

  if (ncol(scores) == 0) {
    stop("No cell type scores found in ", methodName)
  }

  # Set up cell types to plot
  if (is.null(ctypes)) {
    ctypes <- sort(colnames(scores))
  } else {
    missing_ctypes <- ctypes[!ctypes %in% colnames(scores)]
    if (length(missing_ctypes) > 0) {
      warning("Cell types not found in scores matrix: ", paste(missing_ctypes, collapse = ", "))
      ctypes <- ctypes[ctypes %in% colnames(scores)]
    }
  }

  if (length(ctypes) == 0) {
    stop("No valid cell types to plot")
  }

  # Set up color scale limits
  if (freeColScale) {
    lmts <- NULL
  } else {
    lmts <- stats::quantile(scores[, ctypes, drop = FALSE], quants, na.rm = TRUE)
  }

  # Create plots
  message("Creating plots for ", length(ctypes), " cell types...")
  outPlots <- vector("list", length(ctypes))
  names(outPlots) <- ctypes

  # Extract coordinate data once
  coord_data <- as.data.frame(colData(sce)[, coordNames, drop = FALSE])

  for (i in seq_along(ctypes)) {
    ctype <- ctypes[i]

    # Get and process scores for this cell type
    cols <- scores[, ctype]
    if (censQuant < 1) {
      cols <- censor(cols, quant = censQuant)
    }

    # Create plot data
    plot_data <- coord_data
    plot_data$score <- cols
    plot_data <- plot_data[order(plot_data$score), ]  # Order by score for better visualization

    # Create plot
    p <- ggplot(plot_data, aes(x = .data[[coordNames[1]]], y = .data[[coordNames[2]]])) +
      geom_point(aes(colour = .data[["score"]]), size = psize, alpha = 1, stroke = 0) +
      theme_void() +
      scale_colour_gradientn(
        colours = MATLAB_cols,
        limits = lmts,
        name = "Score"
      ) +
      ggtitle(ctype) +
      theme(
        legend.position = legendPosition,
        plot.title = element_text(size = fsize_title, hjust = 0.5, vjust = -1),
        plot.margin = unit(plot_margin, "points")
      )

    outPlots[[i]] <- p
  }

  # Save plots if requested
  if (savePlots) {
    # Create output directory
    full_output_dir <- file.path(outputDir, tissueName)
    if (!dir.exists(full_output_dir)) {
      dir.create(full_output_dir, recursive = TRUE)
      message("Created directory: ", full_output_dir)
    }

    npheno <- length(outPlots)
    nFigPerPlot <- fignrow * figncol
    nFiles <- ceiling(npheno / nFigPerPlot)

    message("Saving ", nFiles, " file(s) to ", full_output_dir)

    for (ii in seq_len(nFiles)) {
      idx_start <- (ii - 1) * nFigPerPlot + 1
      idx_end <- min(ii * nFigPerPlot, npheno)

      filename <- file.path(
        full_output_dir,
        paste0(tissueName, "_", methodName, "_", ii, ".", fileFormat)
      )

      # Create arranged plot
      arranged_plot <- ggpubr::ggarrange(
        plotlist = outPlots[idx_start:idx_end],
        ncol = figncol,
        nrow = fignrow
      )

      # Save with appropriate function based on format
      if (fileFormat == "pdf") {
        ggsave(filename, arranged_plot, width = width, height = height, device = "pdf")
      } else {
        ggsave(filename, arranged_plot, width = width, height = height, dpi = dpi)
      }

      message("Saved: ", basename(filename))
    }
  }

  # Return plots if requested
  if (returnPlots) {
    return(outPlots)
  } else {
    return(invisible(NULL))
  }
}
