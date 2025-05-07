#' Plot PhiSpace scores as a heatmap
#'
#' @param PhiSpaceScore Matrix.
#' @param phenoDict A named list, categorising columns (phenotypes) in `PhiSpaceScore`. Total number of elements have to be identical to number of columns of `PhiSpaceScore`.
#' @param phenotypes Optional. Types of phenotypes. If provided, has to be of the same length as phenoDict.
#' @param queryLabs Character vector.
#' @param refLvls Character vector.
#' @param queryLvls Character vector.
#' @param cluster_columns Logic.
#' @param cluster_rows Logic.
#' @param show_row_names Logic.
#' @param show_column_names Logic.
#' @param row_title_rot Rotation (angle) of row titles, inherited from `ComplexHeatmap::Heatmap`
#' @param ... Additional parameters adjusting the `ComplexHeatmap` object
#' @param phenoDict
#'
#' @return A heatmap of PhiSpace scores.
#' @export
plotPhiSpaceHeatMap <- function(PhiSpaceScore,
                                phenoDict = NULL,
                                phenotypes = NULL,
                                queryLabs = NULL,
                                refLvls = NULL,
                                queryLvls = NULL,
                                # parameters for ComplexHeatmap
                                cluster_columns = FALSE,
                                cluster_rows = FALSE,
                                show_row_names = FALSE,
                                show_column_names = TRUE,
                                row_title_rot = 0,
                                ...){

  if(!is.null(phenoDict)){

    if(!is.null(phenotypes)){

      if(length(phenoDict) != length(phenoDict)) stop("Length of phenotypes has to be the same as phenoDict.")
      names(phenoDict) <- phenotypes
    }

    phenoDict <- data.frame(
      labs = unlist(phenoDict),
      phenotypeCategory = rep(names(phenoDict), sapply(phenoDict, length) )
    )

  } else {

    if(is.null(phenotypes)) phenotypes <- "phenotypes"
    phenoDict <- data.frame(
      labs = colnames(PhiSpaceScore),
      phenotypeCategory = rep(phenotypes[1], ncol(PhiSpaceScore))
    )

  }

  if(is.null(refLvls)) refLvls <- phenoDict$labs

  plot_dat <- PhiSpaceScore[,refLvls] %>% as.data.frame()

  if(!is.null(queryLabs)){

    if(is.null(queryLvls)){

       plot_dat[,"label"] <- factor(queryLabs)
       plot_dat <- dplyr::arrange(plot_dat, plot_dat$label)

    } else {

      plot_dat <- plot_dat %>% dplyr::mutate(label = factor(queryLabs, levels = queryLvls))
      plot_dat <- dplyr::arrange(plot_dat, plot_dat$label)
    }

  }

  if(is.null(queryLabs)){
    rowSplit <- NULL
  } else {
    rowSplit <- plot_dat$label
  }

  ComplexHeatmap::Heatmap(as.matrix(plot_dat[,1:length(refLvls)]),
                          cluster_columns = cluster_columns,
                          cluster_rows = cluster_rows,
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          column_split = phenoDict$phenotypeCategory,
                          row_split = rowSplit, row_title_rot = row_title_rot,
                          ...)


}


