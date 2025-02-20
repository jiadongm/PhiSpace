#' Score matrix plot
#'
#' Density and pairwise scatter plots for visualising score matrices.
#'
#'
#' @param scores Matrix of scores to be plotted
#' @param max_ncomp Default NULL. Number of first components to plot. If specified, will override comp_idx.
#' @param comp_idx Default NULL. Indices of components to plot.
#' @param pointAlpha Alpha value.
#' @param pointSize Point size.
#' @param manualCol Manual specification of colours.
#' @param manualAlpha Manual specification of alpha colours.
#' @param colBy Numeric or charactor vectors to specify colour of points.
#' @param fsize Figure font size.
#' @param returnPlotList Logical. Whether to return individual plots.
#' @param legendTitle Legend title.
#' @param compName Name of the components, default is "comp", so that the 1st column is named comp1 etc.
#'
#' @export
matrixPlot <- function(
    scores,
    max_ncomp = NULL,
    comp_idx = NULL,
    colBy = NULL,
    pointAlpha = NULL,
    pointSize = 1,
    manualCol = NULL,
    manualAlpha = NULL,
    fsize = 14,
    returnPlotList = F,
    legendTitle = "",
    compName = "comp"
  ){

  if(is.null(max_ncomp) & is.null(comp_idx)){
    stop("Need to specify either max_ncomp or comp_idx.")
  }

  if(!is.null(max_ncomp)) comp_idx <- 1:max_ncomp


  if(!all(paste0(compName, comp_idx) %in% colnames(scores))){

    missingComps <-
      paste0(compName, comp_idx)[!(paste0(compName, comp_idx) %in% colnames(scores))]
    stop(paste0("These components are missing from scores: ", missingComps))

  }


  Ngroups <- length(unique(colBy))

  if(length(comp_idx) >= 3){

    # Density plots on diagonal
    out_diag <- vector("list", length(comp_idx))

    for(comp_i in 1:length(comp_idx)){

      var2plot <- paste0(compName, comp_idx[comp_i])

      p <- scores %>%
        ggplot(aes(x = !!sym(var2plot))) +
        geom_density(bw = "sj") +
        theme_bw(base_size = fsize) +
        theme(
          legend.position = "none",
          axis.title.y = element_blank()
        )



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
          ) + labs(
            colour = legendTitle
          )

        )
      )
    }


    # Scatter plots on non-diagonal
    combs <- utils::combn(comp_idx, 2) %>% t()
    out_nondiag <- vector("list", nrow(combs))
    for(comb in 1:nrow(combs)){

      var1 <- paste0(compName, combs[comb, 1])
      var2 <- paste0(compName, combs[comb, 2])

      if(is.null(colBy)){

        p <- scores %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(size = pointSize, stroke = 0) +
          theme_bw(base_size = fsize) +
          theme(
            legend.position = "none",
            axis.ticks = element_blank()
          ) +
          scale_x_continuous(name = NULL, labels = NULL) +
          scale_y_continuous(name = NULL, labels = NULL)


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
          scale_x_continuous(name = NULL, labels = NULL) +
          scale_y_continuous(name = NULL, labels = NULL)

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
    # diag(layoutM) <- (nrow(combs)+1):(nrow(combs)+length(comp_idx))
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

    ## Only density plot is given
    out <-
      scores %>%
      ggplot(aes(x = !! sym(var2plot))) +
      geom_density(bw = "sj") +
      theme_bw(base_size = fsize) +
      ylab("Density")
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
        ) +
        labs(
          colour = legendTitle
        )
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

    if(is.null(colBy)){

      p <- cores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(size = pointSize, stroke = 0)


    } else {

      p <-
        scores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(aes(colour = colBy), size = pointSize, stroke = 0) +
        theme_bw(base_size = fsize) +
        labs(
          colour = legendTitle
        )

      if(!is.null(manualCol)){
        p <- p +
          scale_color_manual(values = manualCol)
      } else {
        if(is.numeric(colBy)) p <- p + scale_colour_gradientn(colours = MATLAB_cols)
      }
    }

  }

  return(p)

}

