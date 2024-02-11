#' Plot PhiSpace scores as a heatmap
#'
#' @param PhiSpaceScore Matrix.
#' @param reference `SingleCellExperiment` object.
#' @param phenotypes Character vector.
#' @param queryLabs Character vector.
#' @param refLvls Character vector.
#' @param queryLvls Character vector.
#' @param cluster_columns Logic.
#' @param cluster_rows Logic.
#' @param show_row_names Logic.
#' @param show_column_names Logic.
#' @param ... Additional parameters adjusting the `ComplexHeatmap` object
#'
#' @return A heatmap of PhiSpace scores.
#' @export
plotPhiSpaceHeatMap <- function(PhiSpaceScore,
                                reference,
                                phenotypes,
                                queryLabs = NULL,
                                refLvls = NULL,
                                queryLvls = NULL,
                                # parameters for ComplexHeatmap
                                cluster_columns = FALSE,
                                cluster_rows = FALSE,
                                show_row_names = FALSE,
                                show_column_names = TRUE,
                                ...){

  phenoDict <-
    data.frame(
      labs = colnames(codeY(reference, phenotypes)),
      phenotypeCategory =
        rep(phenotypes,
            apply(colData(reference)[,phenotypes], 2, function(x) length(unique(x)))
        )
    )

  if(is.null(refLvls)) refLvls <- phenoDict$labs

  plot_dat0 <- PhiSpaceScore[,refLvls]

  if(!is.null(queryLabs)){
    if(is.null(queryLvls)){
      plot_dat <-
        plot_dat0 %>%
        as.data.frame() %>%
        dplyr::arrange(factor(queryLabs))
    } else {
      plot_dat <-
        plot_dat0 %>%
        as.data.frame() %>%
        dplyr::arrange(factor(queryLabs, levels = queryLvls))
    }
  } else {
    plot_dat <- plot_dat0 %>%
      as.data.frame()
  }

  if(is.null(queryLabs)){
    rowSplit <- NULL
  } else {
    rowSplit <- plot_dat$label
  }

  ComplexHeatmap::Heatmap(as.matrix(plot_dat[,1:ncol(plot_dat0)]),
                          cluster_columns = cluster_columns,
                          cluster_rows = cluster_rows,
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          column_split = phenoDict$phenotypeCategory,
                          row_split = rowSplit,
                          heatmap_legend_param = list(title = ""),
                          ...)


}


