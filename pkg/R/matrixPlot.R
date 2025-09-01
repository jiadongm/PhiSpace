#' Score matrix plot
#'
#' Density and pairwise scatter plots for visualising score matrices.
#'
#'
#' @param scores Matrix of scores to be plotted
#' @param max_ncomp Default NULL. Number of first components to plot. If specified, will override comp_idx and var_names.
#' @param comp_idx Default NULL. Indices of components to plot. Will override var_names if specified.
#' @param var_names Default NULL. Character vector of variable names (column names) to plot. Will be overridden by max_ncomp or comp_idx if they are specified.
#' @param pointAlpha Alpha value.
#' @param pointSize Point size.
#' @param manualCol Manual specification of colours.
#' @param manualAlpha Manual specification of alpha colours.
#' @param colBy Numeric or charactor vectors to specify colour of points.
#' @param fsize Figure font size.
#' @param returnPlotList Logical. Whether to return individual plots.
#' @param legendTitle Legend title.
#' @param compName Name of the components, default is "comp", so that the 1st column is named comp1 etc.
#' @param legendPosition Location of the legend when plotting two scores.
#'
#' @export
matrixPlot <- function(
    scores,
    max_ncomp = NULL,
    comp_idx = NULL,
    var_names = NULL,
    colBy = NULL,
    pointAlpha = NULL,
    pointSize = 1,
    manualCol = NULL,
    manualAlpha = NULL,
    fsize = 14,
    returnPlotList = F,
    legendTitle = "",
    compName = "comp",
    legendPosition = c("right", "left",  "bottom", "top", "inside", "none")
){

  legendPosition <- match.arg(legendPosition)

  # Convert scores to matrix and ensure it has proper structure
  scores <- as.matrix(scores)

  # Handle case where scores has no column names
  if(is.null(colnames(scores))) {
    colnames(scores) <- paste0(compName, 1:ncol(scores))
  }

  # Convert to data frame for ggplot
  scores <- as.data.frame(scores)

  # Determine which variables to plot based on priority: max_ncomp > comp_idx > var_names
  if(!is.null(max_ncomp)) {
    # Use first max_ncomp columns
    comp_idx <- 1:max_ncomp
    var_names <- NULL
  } else if(!is.null(comp_idx)) {
    # Use specified column indices
    var_names <- NULL
  } else if(!is.null(var_names)) {
    # Use specified variable names
    # Check if all requested variables exist
    missing_vars <- var_names[!var_names %in% colnames(scores)]
    if(length(missing_vars) > 0) {
      stop(paste0("These variables are not found in the matrix: ", paste(missing_vars, collapse = ", ")))
    }
    # Convert variable names to indices for internal processing
    comp_idx <- match(var_names, colnames(scores))
  } else {
    stop("Need to specify either max_ncomp, comp_idx, or var_names.")
  }

  # If column names don't follow the expected pattern, rename them temporarily
  original_colnames <- colnames(scores)
  temp_colnames <- paste0(compName, 1:ncol(scores))

  # Create a mapping between temporary and original names
  colname_mapping <- stats::setNames(original_colnames, temp_colnames)

  # Temporarily rename columns for internal processing
  colnames(scores) <- temp_colnames

  # Check if requested components exist (now using temporary names)
  requested_cols <- paste0(compName, comp_idx)
  if(!all(requested_cols %in% colnames(scores))){
    stop(paste0("Component indices out of range. Matrix only has ", ncol(scores), " columns."))
  }

  Ngroups <- length(unique(colBy))

  if(length(comp_idx) >= 3){

    # Density plots on diagonal
    out_diag <- vector("list", length(comp_idx))

    for(comp_i in 1:length(comp_idx)){

      var2plot <- paste0(compName, comp_idx[comp_i])
      # Get original column name for axis label
      original_name <- colname_mapping[var2plot]

      p <- scores %>%
        ggplot(aes(x = !!sym(var2plot))) +
        geom_density(bw = "sj") +
        theme_bw(base_size = fsize) +
        theme(
          legend.position = "none",
          axis.title.y = element_blank()
        ) +
        xlab(original_name)  # Use original column name for label

      if(is.null(colBy)){

        p <- p +
          geom_jitter(aes(y = 0), height = diff(layer_scales(p)$y$range$range)/20, size = pointSize, stroke = 0)
        out_diag[[comp_i]] <- p

      } else {

        p <- p +
          geom_jitter(aes(y = 0, colour = colBy),
                      height = diff(layer_scales(p)$y$range$range)/20)

        if(!is.null(manualCol)){
          p <- p + scale_color_manual(values = manualCol)
        } else {

          if(is.numeric(colBy)) p <- p + scale_colour_gradientn(colours = MATLAB_cols)

        }

        out_diag[[comp_i]] <- p

      }
    }

    # Get the legend
    if(!is.null(colBy)){
      suppressWarnings(
        p_legend <- cowplot::get_legend(
          p + theme(
            legend.position = "right"
          ) + labs(colour = legendTitle)

        )
      )
    }

    # Scatter plots on non-diagonal
    combs <- utils::combn(comp_idx, 2) %>% t()
    out_nondiag <- vector("list", nrow(combs))
    for(comb in 1:nrow(combs)){

      var1 <- paste0(compName, combs[comb, 1])
      var2 <- paste0(compName, combs[comb, 2])

      # Get original column names for axis labels
      original_name1 <- colname_mapping[var1]
      original_name2 <- colname_mapping[var2]

      if(is.null(colBy)){

        p <- scores %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(size = pointSize, stroke = 0) +
          theme_bw(base_size = fsize) +
          theme(
            legend.position = "none",
            axis.ticks = element_blank()
          ) +
          xlab(original_name1) + ylab(original_name2)  # Use original names

      } else {

        p <-
          scores %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(aes(colour = colBy), size = pointSize, stroke = 0) +
          theme_bw(base_size = fsize) +
          theme(
            legend.position = "none",
            axis.ticks = element_blank()
          ) +
          xlab(original_name1) + ylab(original_name2)  # Use original names

        if(!is.null(manualCol)){
          p <- p + scale_color_manual(values = manualCol)
        } else {

          if(is.numeric(colBy)) p <- p + scale_colour_gradientn(colours = MATLAB_cols)
        }
      }

      out_nondiag[[comb]] <- p

    }

    # Arrange plots
    out <- c(out_nondiag, out_diag)
    # Arrange diagonal
    diagIdx <- 1
    toAdd <- length(comp_idx)
    for(kk in 1:length(comp_idx)){

      out[[diagIdx]] <- out_diag[[kk]]
      diagIdx <- diagIdx + toAdd
      toAdd <- toAdd - 1
    }
    # Non-diagonal, arrange column by column
    startIdxMat <- 2
    toAdd <- length(comp_idx) - 2
    startIdxNondiag <- 1
    for(kk in 1:(length(comp_idx)-1)){
      endIdxMat <- startIdxMat + toAdd
      endIdxNondiag <- startIdxNondiag + toAdd
      out[startIdxMat:endIdxMat] <- out_nondiag[startIdxNondiag:endIdxNondiag]
      startIdxMat <- endIdxMat + 2
      startIdxNondiag <- endIdxNondiag + 1
      toAdd <- toAdd - 1
    }

    layoutM <- matrix(NA, length(comp_idx), length(comp_idx))
    layoutM[lower.tri(layoutM, diag = T)] <- 1:length(out)
    if(!is.null(colBy)){
      out[[length(out) + 1]] <- p_legend
      layoutM[1, length(comp_idx)] <- length(out)
    }
    p <- gridExtra::grid.arrange(grobs = out, layout_matrix = layoutM)

    ## Return
    if(returnPlotList){
      return(list(matrixPlot = p, plotList = out))
    } else {
      return(p)
    }

  } else if (length(comp_idx) == 1){

    var2plot <- paste0(compName, comp_idx)
    # Get original column name for axis label
    original_name <- colname_mapping[var2plot]

    ## Only density plot is given
    out <-
      scores %>%
      ggplot(aes(x = !! sym(var2plot))) +
      geom_density(bw = "sj") +
      theme_bw(base_size = fsize) +
      ylab("Density") +
      xlab(original_name)  # Use original column name

    if(is.null(colBy)){
      out <-
        out +
        geom_jitter(
          aes(y=0),
          height = diff(layer_scales(out)$y$range$range)/20,
          size = pointSize,
          shape = 16,
          stroke = 0
        )

    } else {
      out <- out +
        geom_jitter(
          aes(y=0, colour = colBy),
          height = diff(layer_scales(out)$y$range$range)/20,
          size = pointSize,
          shape = 16,
          stroke = 0
        ) + labs(colour = legendTitle)
    }

    if(!is.null(manualCol)){
      out <- out + scale_color_manual(values = manualCol)
    } else {
      if(is.numeric(colBy)) out <- out + scale_colour_gradientn(colours = MATLAB_cols)
    }

    # Colour scale
    if(is.null(manualAlpha)){

      if(is.null(pointAlpha)) pointAlpha <- 1

      out <- out +
        scale_alpha_manual(values = rep(pointAlpha, Ngroups))

    } else {
      out <- out +
        scale_alpha_manual(values = manualAlpha)
    }

    return(out)

  } else {

    ## Plot 2 components

    var1 <- paste0(compName, comp_idx[1])
    var2 <- paste0(compName, comp_idx[2])

    # Get original column names for axis labels
    original_name1 <- colname_mapping[var1]
    original_name2 <- colname_mapping[var2]

    if(is.null(colBy)){

      p <- scores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(size = pointSize, stroke = 0) +
        theme_bw(base_size = fsize) +
        theme(legend.position = legendPosition) +
        xlab(original_name1) + ylab(original_name2)  # Use original names

    } else {

      p <-
        scores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(aes(colour = colBy), size = pointSize, stroke = 0) +
        theme_bw(base_size = fsize) + labs(colour = legendTitle) +
        theme(legend.position = legendPosition) +
        xlab(original_name1) + ylab(original_name2)  # Use original names

      if(!is.null(manualCol)){
        p <- p + scale_color_manual(values = manualCol)
      } else {
        if(is.numeric(colBy)) p <- p + scale_colour_gradientn(colours = MATLAB_cols)
      }
    }
  }

  return(p)
}
