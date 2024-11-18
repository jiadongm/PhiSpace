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
#' @param groupKey Column name of plot_dat storing group information.
#'
#' @export
matrixPlot <- function(
    scores,
    max_ncomp = NULL,
    comp_idx = NULL,
    groupKey = NULL,
    pointAlpha = NULL,
    pointSize = 1,
    manualCol = NULL,
    manualAlpha = NULL,
    fsize = 14
  ){

  if(is.null(max_ncomp) & is.null(comp_idx)){
    stop("Need to specify either max_ncomp or comp_idx.")
  }

  if(!is.null(max_ncomp)) comp_idx <- 1:max_ncomp


  if(!all(paste0("comp", comp_idx) %in% colnames(scores))){

    missingComps <-
      paste0("comp", comp_idx)[!(paste0("comp", comp_idx) %in% colnames(scores))]
    stop(paste0("These components are missing from scores: ", missingComps))

  }


  if(length(comp_idx) >= 3){

    # Density plots on diagonal
    out_diag <- vector("list", length(comp_idx))

    for(comp_i in 1:length(comp_idx)){

      var2plot <- paste0("comp", comp_idx[comp_i])

      p <-
        scores %>%
        ggplot(aes(x = !! sym(var2plot))) +
        geom_density(bw = "sj") +
        scale_y_continuous(name = NULL)



      if(is.null(groupKey)){

        p <- p +
          geom_jitter(aes(y = 0),
                      height = diff(layer_scales(p)$y$range$range)/20)
        out_diag[[comp_i]] <- p +
          theme_bw(base_size = fsize)

      } else {

        p <- p +
          geom_jitter(aes(y = 0, colour = !! sym(groupKey)),
                      height = diff(layer_scales(p)$y$range$range)/20)

        if(!is.null(manualCol)){
          p <- p +
            scale_color_manual(values = manualCol)
        }

        out_diag[[comp_i]] <- p +
          theme_bw(base_size = fsize) +
          theme(legend.position = "none")



      }



    }

    # Get the legend
    if(!is.null(groupKey)) p_legend <- cowplot::get_legend(p)


    # Scatter plots on non-diagonal
    combs <- utils::combn(comp_idx, 2) %>% t()
    out_nondiag <- vector("list", nrow(combs))
    for(comb in 1:nrow(combs)){

      var1 <- paste0("comp", combs[comb, 1])
      var2 <- paste0("comp", combs[comb, 2])

      if(is.null(groupKey)){

        p <-
          scores %>%
          ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
          geom_point(size = pointSize, stroke = 0) +
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
          geom_point(aes(colour = !! sym(groupKey)), size = pointSize, stroke = 0) +
          theme(
            legend.position = "none",
            axis.ticks = element_blank()
          ) +
          scale_x_continuous(name = NULL, labels = NULL) +
          scale_y_continuous(name = NULL, labels = NULL)

        if(!is.null(manualCol)){
          p <- p + scale_color_manual(values = manualCol)
        }
      }

      p <- p + theme_bw(base_size = fsize)

      out_nondiag[[comb]] <- p

    }


    # Arrange plots
    out <- c(out_nondiag, out_diag)
    if(!is.null(groupKey)){
      out[[length(out) + 1]] <- p_legend
    }
    layoutM <- matrix(NA, length(comp_idx), length(comp_idx))
    layoutM[upper.tri(layoutM, diag = F)] <- 1:nrow(combs)
    diag(layoutM) <- (nrow(combs)+1):(nrow(combs)+length(comp_idx))
    if(!is.null(groupKey)) layoutM[ceiling(length(comp_idx)/2)-1, ceiling(length(comp_idx)/2)+1] <- length(out)
    p <- gridExtra::grid.arrange(grobs = out, layout_matrix = layoutM)

  } else if (length(comp_idx) == 1){

    var2plot <- paste0("comp", comp_idx)

    ## Only density plot is given
    out <-
      scores %>%
      ggplot(aes(x = !! sym(var2plot))) +
      geom_density(bw = "sj")
    if(is.null(groupKey)){
      out <-
        out +
        geom_jitter(aes(y=0),
                    height = diff(layer_scales(out)$y$range$range)/20,
                    size = pointSize,
                    shape = 16)
    } else {
      out <-
        out +
        geom_jitter(aes(y=0, colour = !! sym(groupKey), alpha = !! sym(groupKey)),
                    height = diff(layer_scales(out)$y$range$range)/20,
                    size = pointSize,
                    shape = 16)
    }

    if(!is.null(manualCol)){
      out <- out +
        scale_color_manual(values = manualCol)
    }

    if(is.null(manualAlpha)){

      Ngroups <- length(unique(scores[,groupKey]))

      if(is.null(pointAlpha)) pointAlpha <- 1

      out <- out +
        scale_alpha_manual(values = rep(pointAlpha, Ngroups))

    } else {
      out <- out +
        scale_alpha_manual(values = manualAlpha)
    }


  } else {

    ## Plot 2 components

    var1 <- paste0("comp", comp_idx[1])
    var2 <- paste0("comp", comp_idx[2])

    if(is.null(groupKey)){

      p <-
        scores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(size = pointSize)


    } else {

      p <-
        scores %>%
        ggplot(aes(x = !! sym(var1), y = !! sym(var2))) +
        geom_point(aes(colour = !! sym(groupKey)), size = pointSize)

      if(!is.null(manualCol)){
        p <- p +
          scale_color_manual(values = manualCol)
      }
    }

  }

  return(p)

}

