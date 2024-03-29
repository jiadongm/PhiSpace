#' Plot Sankey diagram for three hard classification results.
#'
#' @param class1 Vector.
#' @param class2 Vector.
#' @param class3 Vector.
#' @param add Logic.
#' @param fontsize Numeric.
#'
#' @export
plotSankey3 <- function(class1, class2, class3, add = TRUE, fontsize = 12){


  # er_re12 <- classErr(class1, class2)
  # er_re23 <- classErr(class3, class2)

  if(add){
    class1 <- paste0(class1, "-")
    class2 <- paste0(class2, "_")
  }


  # Define a flow incidence matrix, from rows to cols
  rowNams <- colNams <- unique(c(class1, class2, class3))

  incidMat <- matrix(0, nrow = length(rowNams), ncol = length(colNams),
                     dimnames = list(rowNams, colNams))


  confTab12 <- table(class1, class2)
  confTab23 <- table(class2, class3)

  confTab <- confTab12
  for(i in 1:nrow(confTab)){
    fromClass <- rownames(confTab)[i]
    for(j in 1:ncol(confTab)){
      toClass <- colnames(confTab)[j]
      incidMat[fromClass, toClass] <- confTab[i, j]
    }
  }
  confTab <- confTab23
  for(i in 1:nrow(confTab)){
    fromClass <- rownames(confTab)[i]
    for(j in 1:ncol(confTab)){
      toClass <- colnames(confTab)[j]
      incidMat[fromClass, toClass] <- confTab[i, j]
    }
  }

  # Following code from https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
  sankeyLinks <- as.data.frame(incidMat) %>%
    dplyr::mutate(source = rownames(incidMat)) %>%
    tidyr::pivot_longer(!source, names_to = "target", values_to = "value")
  sankeyLinks <- sankeyLinks %>%
    dplyr::filter(sankeyLinks$value != 0)

  sankeyNodes <- data.frame(
    name = c(as.character(sankeyLinks$source),
             as.character(sankeyLinks$target)) |> unique()
  )

  sankeyLinks$IDsource <- match(sankeyLinks$source, sankeyNodes$name)-1
  sankeyLinks$IDtarget <- match(sankeyLinks$target, sankeyNodes$name)-1

  p <-
    suppressMessages(networkD3::sankeyNetwork(Links = sankeyLinks, Nodes = sankeyNodes,
                                              Source = "IDsource", Target = "IDtarget",
                                              Value = "value", NodeID = "name",
                                              sinksRight = FALSE, fontSize = fontsize))
  print(p)

  # return(list(er12 = er_re12$err, er23 = er_re23$err))
}
