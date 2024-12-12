MATLAB_cols <- c(rgb(54, 70, 157, maxColorValue = 255),
                 rgb(61, 146, 185, maxColorValue = 255),
                 rgb(126, 203, 166, maxColorValue = 255),
                 rgb(204, 234, 156, maxColorValue = 255),
                 rgb(249, 252, 181, maxColorValue = 255),
                 rgb(255, 226, 144, maxColorValue = 255),
                 rgb(253, 164, 93, maxColorValue = 255),
                 rgb(234, 95, 70, maxColorValue = 255),
                 rgb(185, 30, 72, maxColorValue = 255))

censor <- function(vec, quant = 0.05){

  vec <- as.vector(vec)

  cutoff <- stats::quantile(vec, 1-quant)
  vec[vec <= cutoff] <- cutoff
  vec
}

#' Plot formatting for presentation and publication
#'
#' @param p A ggplot2 object.
#' @param theme_overall Theme of p.
#' @param base_size base_size of theme.
#' @param legend.position Legend position.
#' @param x.text.blank Logical. Set x axis text to blank.
#' @param y.text.blank Logical. Set y axis text to blank.
#' @param axis.title.bold Logical. Set title to be bold.
#' @param legend.title.bold Logical. Set legend title to be bold.
#' @param axis.title.x.blank Logical. Set x axis title to be bold.
#' @param axis.title.y.blank Logical. Set y axis title to be bold.
#' @param legend.title.blank Logical. Set legend title to be blank.
#' @param legend.spacing.x Legend spacing.
#' @param legend.key.spacing Legend key spacing.
#' @param legend.key.spacing.y Legend key spacing on vertical direction.
#' @param legend.box.margin Legend box margin.
#' @param legend.box.spacing Legend box spacing.,
#' @param legend.key.size Legend key size.
#'
#' @return A ggplot2 object.
basicPlotFormat <- function(
    p, theme_overall = "bw",
    base_size = 8, legend.position = "top",
    x.text.blank = F, y.text.blank = F,
    axis.title.bold = T, legend.title.bold = T,
    axis.title.x.blank = F, axis.title.y.blank = F,
    ## Legend format
    legend.title.blank = F,
    legend.spacing.x = 0,
    legend.key.spacing = 0,
    legend.key.spacing.y = 0,
    legend.box.margin = 0,
    legend.box.spacing = 0,
    legend.key.size = 8
){

  theme_fun <- switch (
    theme_overall,
    classic = theme_classic,
    ggplot2 = theme_gray,
    bw = theme_bw,
    linedraw = theme_linedraw,
    light = theme_light,
    dark = theme_dark,
    minimal = theme_minimal,
    void = theme_void
  )

  p <- p +
    theme_fun(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend.position,
      legend.key = element_blank(),
      legend.spacing.x = unit(legend.spacing.x, "pt"),
      legend.key.spacing = unit(legend.key.spacing, "pt"),
      legend.key.spacing.y = unit(legend.key.spacing.y, "pt"),
      legend.box.margin = margin(
        legend.box.margin,
        legend.box.margin,
        legend.box.margin,
        legend.box.margin
      ),
      legend.box.spacing = unit(legend.box.spacing, "pt"),
      legend.background = element_blank(),
      legend.key.size = unit(legend.key.size, "pt")
    )

  if(axis.title.bold) p <- p + theme(axis.title = element_text(face = "bold"))
  if(!legend.title.blank){

    if(legend.title.bold) p <- p + theme(legend.title = element_text(face = "bold"))
  } else {

    p <- p + theme(legend.title = element_blank())
  }

  if(x.text.blank){
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45))
  }

  if(y.text.blank){
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
      )
  }

  if(axis.title.x.blank) p <- p + theme(axis.title.x = element_blank())
  if(axis.title.y.blank) p <- p + theme(axis.title.y = element_blank())

  # Legend
  # p <- p +
  #   theme(
  #     legend.key.size = unit(legend.key.size, "cm")
  #   )

  # p <- p +
  #     theme(
  #       axis.text = element_text(size = axis.text.x.fsize),
  #       axis.title.x = element_text(
  #         size = axis.title.x.fsize,
  #         margin = margin(
  #           axis.title.x.margin,
  #           axis.title.x.margin,
  #           axis.title.x.margin,
  #           axis.title.x.margin
  #         )
  #       ),
  #       axis.title.y = element_text(
  #         size = axis.title.y.fsize,
  #         margin = margin(
  #           axis.title.y.margin,
  #           axis.title.y.margin,
  #           axis.title.y.margin,
  #           axis.title.y.margin
  #         )
  #       ),
  #       legend.text = element_text(size = legend.text.fsize),
  #       legend.title = element_text(size = legend.title.fsize, face = "bold"),
  #       legend.position = legend.position,
  #       legend.title.position = legend.title.position,
  #       legend.key = element_blank(),
  #       legend.spacing.x = unit(legend.spacing.x, "pt"),
  #       legend.key.spacing = unit(legend.key.spacing, "pt"),
  #       legend.key.spacing.y = unit(legend.key.spacing.y, "pt"),
  #       legend.box.margin = margin(
  #         legend.box.margin,
  #         legend.box.margin,
  #         legend.box.margin,
  #         legend.box.margin
  #       ),
  #       legend.box.spacing = unit(legend.box.spacing, "pt"),
  #       legend.background = element_blank()
  #     ) +
  #     guides(
  #       color = guide_legend(
  #         override.aes = list(
  #           size = legend.obj.size
  #         ),
  #         nrow = 1
  #       ),
  #       shape = guide_legend(
  #         override.aes = list(
  #           size = legend.obj.size
  #         )
  #       )
  #     )

  return(p)
}





#' Visualising features (eg gene expression), metadata (eg cell type) or reducedDim component (eg PhiSpace cell type score) in 2D tissue space.
#'
#' The 2D spatial coordinates (x and y) should be stored as two columns in colData.
#'
#' @param obj Input object; either `SingleCellExperiemnt` or `SpatialExperiment`. If `SingleCellExperiemnt`, then `colData(obj)` has to contain columns named by `x_coord` and `y_coord` below.
#' @param x_coord x coordinate name of a column of colData (if `type(ojb) == SingleCellExperiemnt`).
#' @param y_coord y coordinate name of a column of colData (if `type(ojb) == SingleCellExperiemnt`).
#' @param ptSize Size of point, representing spatial cell like objects (segmented cells, spots etc).
#' @param ptShape Shape of point.
#' @param groupBy Name of a metadata column to colour points by.
#' @param feature Name of a feature (eg gene) to colour points by.
#' @param assay2use Which assay to use when feature is specified.
#' @param reducedDim Name of reducedDim column to colour points by.
#' @param reducedDim2use Name of reducedDim layer to use when reducedDim is specified.
#' @param censor Logical. Whether to set smaller colour values to zero. Default is `FALSE`. (Recommend to set to `TRUE` when visualise PhiSpace scores,)
#' @param quant A value between 0 and 1. Colour values below `quant` will be set to zero if `censor=TRUE`.
#' @param legend.position Position of legened, one of "right", "bottome", "left" and "top".
#' @param legend.symb.size Legend symbol size (for discrete legend).
#' @param legend.title Alternative legend name to replace the original.
#' @param fsize Base font size of figure.
#' @param reOrder Logical. Whether to reorder points according to their values (ascending) or not. Set to be `TRUE` to avoid overplotting.
#' @param ... Other arguments to pass to `basicPlotFormat`.
#'
#' @return A ggplot2 object.
#' @export
#'
VizSpatial <- function(
    obj,
    x_coord = "x",
    y_coord = "y",
    ptSize = 2,
    ptShape = 16,
    groupBy = NULL,
    feature = NULL,
    assay2use = NULL,
    reducedDim = NULL,
    reducedDim2use = "PhiSpace",
    censor = FALSE, # censor tail PhiSpace scores (see cell2location)
    quant = 0.5,
    ## Plot format
    legend.position = "right",
    legend.symb.size = NULL,
    legend.title = NULL,
    fsize = 14,
    reOrder = F,
    ...
){

  if(class(obj) == "SpatialExperiment"){

    coordNames <- spatialCoordsNames(obj)

    if(is.null(coordNames)){
      warning("Coordinate names nonexistent; will use generic names.")
      coordNames <- paste0("coord", 1:length(coordNames))
      spatialCoordsNames(obj) <- coordNames
    }

    plot_dat <- colData(obj) %>% as.data.frame()

    # If coordNames already present in colData(obj), then delete those columns from colData
    if(any(coordNames %in% colnames(colData(obj)))){

      plot_dat[,coordNames[which(coordNames %in% colnames(colData(obj)))]] <- NULL
    }
    plot_dat <- cbind(plot_dat, spatialCoords(obj))

    if(ncol(spatialCoords(obj)) > 2) message("There are more than 2 spatial coordinates in obj, only using the first two as x and y coordinates.")

    x_coord = coordNames[1]
    y_coord = coordNames[2]

  } else {

    coordNames <- c(x_coord, y_coord)

    if(class(obj) != "SingleCellExperiment") stop("The input obj has to be either SpatialExperiment object or SingleCellExperiment.")
    plot_dat <- colData(obj) %>% as.data.frame()

    if(!all(coordNames %in% colnames(colData(obj)))) stop("Not all spatial coordinates are present in colData of obj.")
  }



  if(!is.null(feature)){

    if(!(feature %in% rownames(obj))) stop("Feature not in obj.")

    if(is.null(assay2use)) assay2use <- "counts"

    plot_dat[[feature]] <- as.numeric(assay(obj, assay2use)[feature, ])

    if(reOrder){
      plot_dat <- dplyr::arrange(plot_dat, !!sym(feature))

    }

  } else if (!is.null(reducedDim)) {

    if(!(reducedDim2use %in% reducedDimNames(obj))) stop("Non-existent reducedDim2use.")

    if(censor){

      plot_dat[[reducedDim]] <- censor(
        as.numeric(reducedDim(obj, reducedDim2use)[,reducedDim]),
        quant = quant
      )
    } else {

      plot_dat[[reducedDim]] <- as.numeric(reducedDim(obj, reducedDim2use)[,reducedDim]) %>% matrix()
    }

    if(reOrder){
      plot_dat <-  dplyr::arrange(plot_dat, !!sym(reducedDim))
    }
  }


  if(!is.null(groupBy)){

    if(reOrder){
      plot_dat <- dplyr::arrange(plot_dat, !!sym(groupBy))
    }
  }

  p <- plot_dat %>%
    ggplot(
      aes(
        x = !!sym(x_coord),
        y = !!sym(y_coord)
      )
    )

  if(!is.null(groupBy)){

    p <- p +
      geom_point(
        aes(
          colour = !!sym(groupBy)
        ),
        size = ptSize, shape = ptShape,
        stroke = 0
      )

    if(is.numeric(plot_dat[,groupBy])){
      p <- p + scale_colour_gradientn(colours = MATLAB_cols)
    }

  } else if (!is.null(feature)) {

    p <- p +
      geom_point(
        aes(
          colour = !!sym(feature)
        ),
        size = ptSize, shape = ptShape,
        stroke = 0
      ) + scale_colour_gradientn(
        colours = MATLAB_cols
      )
  } else if (!is.null(reducedDim)){

    p <- p +
      geom_point(
        aes(
          colour = !!sym(reducedDim)
        ),
        size = ptSize, shape = ptShape,
        stroke = 0
      ) + scale_colour_gradientn(
        colours = MATLAB_cols
      )
  } else {

    p <- p + geom_point(size = ptSize, shape = ptShape)
  }

  p <- basicPlotFormat(
    p, legend.position = legend.position,
    x.text.blank = T, y.text.blank = T,
    axis.title.x.blank = T, axis.title.y.blank = T, base_size = fsize, ...
  )

  if(!is.null(legend.symb.size)){

    p <- p +
      guides(
        colour = guide_legend(
          override.aes = list(size = legend.symb.size)
        )
      )
  }

  if(!is.null(legend.title)){

    p <- p + guides(colour=guide_legend(title=legend.title))
  }

  return(p)
}
