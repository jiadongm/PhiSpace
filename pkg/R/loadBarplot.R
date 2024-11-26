#' Plot loadings of dimension reduction as bar plot.
#'
#' @param Loadings A data.frame or an object (eg matrix) that's convertable to data.frame: columns are components and rows are features.
#' @param comp Character. Which component to plot.
#' @param showInt Logical. Weather to show cell type interaction (eg in PhiSpace ST cell type co-presence analysis).
#' @param absVal Logical. Rank loadings by absolute values or not.
#' @param showNeg Logical. Show negative loadings or not.
#' @param nfeat Number of top loadings to show.
#' @param fsize Font size.
#' @param xlab x axis title.
#'
#' @return A ggplot2 object.
#' @export
loadBarplot <- function(
    Loadings, comp = "comp1", showInt = F, absVal = T, showNeg = F,
    nfeat = 30, fsize = 6, xlab = ""
){

  if(absVal){

    plot_dat <- Loadings %>% as.data.frame() %>%
      dplyr::arrange(desc(abs(!!sym(comp))))  %>%
      dplyr::slice_head(n = nfeat) %>%
      dplyr::arrange((!!sym(comp)))
  } else {

    if(showNeg){
      plot_dat <- Loadings %>% as.data.frame() %>%
        dplyr::arrange(desc(-(!!sym(comp))))  %>%
        dplyr::slice_head(n = nfeat) %>%
        dplyr::arrange((!!sym(comp)))
    } else {

      plot_dat <- Loadings %>% as.data.frame() %>%
        dplyr::arrange(desc(!!sym(comp)))  %>%
        dplyr::slice_head(n = nfeat) %>%
        dplyr::arrange((!!sym(comp)))
    }
  }

  if(showInt){
    rname_list <- lapply(
      strsplit(rownames(plot_dat), "<->"),
      function(x) sort(x)
    )
    rname_new <- sapply(
      rname_list,
      function(x) paste(x[1], x[2], sep = "<->")
    )
    rownames(plot_dat) <- rname_new
  }

  p <- plot_dat %>%
    dplyr::mutate(
      interaction = factor(
        rownames(plot_dat), levels = rownames(plot_dat)
      )
    ) %>%
    ggplot() +
    geom_bar(
      aes(x = !!sym(comp), y = interaction, fill = !!sym(comp)),
      stat = "identity"
    ) +
    xlab(xlab) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_light(base_size = fsize) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  return(p)
}
