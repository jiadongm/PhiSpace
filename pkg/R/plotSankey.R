#' Plot Sankey diagram for two classification results.
#'
#' @param class1
#' @param class2
#' @param fontsize
#'
#'
#' @export
plotSankey <- function(class1, class2, fontsize = 12){

  confTab <- table(class1, class2)

  # Plot Sankey Network Diagram
  colLinks <- c('source', 'target', 'value')
  sankeyLinks <- matrix(ncol = length(colLinks), nrow = 0)
  colnames(sankeyLinks) <- colLinks
  for(i in 1:nrow(confTab)){
    for(j in 1:ncol(confTab)){
      if(confTab[i, j]>0){
        sankeyLinks <- rbind(sankeyLinks,
                             # - 1 is for zero index
                             c(i - 1, nrow(confTab) + j - 1, confTab[i, j]))
      }
    }
  }
  sankeyNodes <- data.frame(name = c(rownames(confTab), colnames(confTab)))
  p <- networkD3::sankeyNetwork(Links = as.data.frame(sankeyLinks),
                                Nodes = sankeyNodes, NodeID = 'name',
                                Source = 'source', Target = 'target',
                                Value = 'value',
                                fontSize = fontsize)
  print(p)
  # return(classErr(class1, class2)$err)
}
