# Matlab blue, yellow, red continuous colour scale
MATLAB_cols <- c(rgb(54, 70, 157, maxColorValue = 255),
                 rgb(61, 146, 185, maxColorValue = 255),
                 rgb(126, 203, 166, maxColorValue = 255),
                 rgb(204, 234, 156, maxColorValue = 255),
                 rgb(249, 252, 181, maxColorValue = 255),
                 rgb(255, 226, 144, maxColorValue = 255),
                 rgb(253, 164, 93, maxColorValue = 255),
                 rgb(234, 95, 70, maxColorValue = 255),
                 rgb(185, 30, 72, maxColorValue = 255))
colScale <- ggplot2::scale_color_gradientn(colours = MATLAB_cols)


clust_cols <- c(
  "1" = "#FDCDAC",
  "2" = "#F4CAE4",
  "3" = "#CBD5E8",
  "4" = "#E6F5C9",
  "5" = "gray",
  "6" = "#1F78B4",
  "7" = "#B3E2CD",
  "8" = "#FFF2AE",
  "9" = "#F1E2CC",
  "10" = "#A6CEE3",
  background = "grey40"
)

barcodeBiColCode <- c(
  Nonbarcoded =  "#377EB8", Barcoded = "#E41A1C"
)

    


## Customised plot function
tempDenPlot <- function(
    xx, 
    xx_eval, 
    xylimits, 
    ch = NULL, # convext hull defining the spleen section
    printPlot = TRUE, 
    psize = 0.5, 
    bcName = "barcode",  # barcode name
    spleenOnly = F
){
  
  xlimts <- xylimits[,1]
  ylimts <- xylimits[,2]
  
  kde_res <- ks::kde(xx, eval.points = xx_eval)
  
  if(spleenOnly){
    
    p <- xx_eval %>%
      as.data.frame() %>%
      mutate(
        cols =kde_res$estimate,
        isin_data = isin_data
      ) %>%
      filter(
        isin_data
      ) 
  } else {
    
    p <- xx_eval %>%
      as.data.frame() %>%
      mutate(
        cols =kde_res$estimate
      ) 
  }
  
  p <- p %>%
    arrange(cols) %>%
    ggplot(aes(x, y)) +
    geom_point(aes(colour = cols), size = psize, alpha = 1, stroke = 0) +
    theme_void(base_size = 6) +
    ggtitle(bcName) +
    xlim(xlimts) +
    ylim(ylimts) +
    # scale_colour_gradientn(colours = MATLAB_cols) +
    theme(
      legend.position = "none",
      title = element_text(hjust = 0.5, vjust = -1)
    ) 
  
  if(!is.null(ch)){
    p <- p +
      geom_polygon(
        data = ch,
        aes(x, y),
        fill = NA,
        alpha = 1,
        colour = "black"
      )
  }
  
  if(printPlot) print(p)
  
  return(list(
    kde_res = kde_res,
    p = p
  ))
}






tempBarPlot_byClust <- function(clustNum){
  
  pDat <- regCoef %>%
    as.data.frame() %>%
    mutate(
      celltype = celltype
    ) %>%
    pivot_longer(
      `10`:`9`,
      names_to = "cluster",
      values_to = "coef"
    ) 
  celltypeImp <- pDat %>% 
    group_by(
      celltype
    ) %>%
    summarise(
      celltype_imp = sum(abs(coef))
    ) %>%
    as.data.frame() %>%
    arrange(-celltype_imp)
  p <- pDat %>%
    mutate(
      celltype = factor(
        pDat$celltype,
        levels = celltypeImp$celltype
      )
    ) %>%
    ggplot(
      aes(celltype, coef)
    ) +
    geom_bar(
      aes(
        fill = cluster
      ), 
      stat = "identity"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, face = "bold"),
      axis.title.x = element_blank()
    ) +
    scale_fill_manual(values = clust_cols) 
  
  return(p)
}




tempBoxPlot <- function(plot_dat, cType, fsize){
  
  plot_dat %>%
    ggplot(aes(cluster, !!sym(cType))) +
    geom_boxplot(
      aes(
        fill = cluster
      ),
      outlier.size = 0.5, lwd = 0.2, outlier.stroke = 0, outlier.alpha = 0.5
    ) +
    scale_fill_manual(values = clust_cols) +
    theme_bw(
      base_size = fsize
    ) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      plot.title = element_text(size = fsize, hjust = 0.5, vjust = -1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    ggtitle(cType) +
    ylim(-1,1)
}

tempSaveBox <- function(
    plot_dat,
    cTypes,
    tissueName = "DawsonMouseSpleenAML",
    savePlots = TRUE,
    returnPlot = FALSE,
    # Save plots parameters
    width = 7.5,
    height = 7.5,
    fignrow = 7,
    figncol = 10,
    fsize = 6
){
  
  outPlots <- vector("list", length(cTypes)) %>% `names<-`(cTypes)
  for(ii in 1:length(cTypes)){
    
    cType <- cTypes[ii]
    outPlots[[ii]] <- tempBoxPlot(plot_dat, cType, fsize) 
  }
  
  if(savePlots){
    
    npheno <- length(outPlots)
    nFigPerPlot <- fignrow * figncol
    for(ii in 1:ceiling(npheno/nFigPerPlot)){
      
      idx_start <- (ii-1)*nFigPerPlot + 1
      idx_end <- min(ii*nFigPerPlot, npheno)
      ggsave(
        paste0(
          "figs/cTypeEnrichBoxplots/",
          tissueName, "_", ii, ".png"
        ),
        ggpubr::ggarrange(
          plotlist = outPlots[idx_start:idx_end], 
          ncol = figncol, nrow = fignrow
        ), 
        width = width, height = height
      )
    }
  }
  if(returnPlot) return(outPlots)
}


tempClustPlot <- function(plot_dat, spleenOnly){
  
  if(spleenOnly) plot_dat <- plot_dat %>% filter(cell_id %in% colnames(query))
    
  p <- plot_dat %>%
    ggplot(
      aes(
        !!sym(xCoordName), !!sym(yCoordName), colour = clusters
      )
    ) +
    geom_point(
      size = 0.5, stroke = 0
    ) +
    theme_void(base_size = 6) +
    xlim(xylimits[,1]) +
    ylim(xylimits[,2]) +
    scale_colour_manual(values = clust_cols) +
    geom_polygon(
      data = ch,
      aes(x, y),
      fill = NA,
      alpha = 1,
      colour = "black"
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 2)
      )
    ) +
    theme(
      legend.key.spacing = unit(0, "pt"),
      legend.title = element_text(face = "bold")
    )
  return(p)
}




tempPlot <- function(kclust = 2, comps = c(1, 2), filter = NULL){
  
  clusts <- clust_list[[kclust - 1]] %>% as.character()
  names(clusts) <- rownames(barAssay)
  idx <- intersect(
    names(clusts),
    colnames(query)
  )
  
  if(!is.null(filter)){
    
    idx <- idx[(clusts[idx] %in% filter)]
  }
  
  clusts <- clusts[idx]
  
  
  
  comps_1 <- paste0("comp", comps[1])
  comps_2 <- paste0("comp", comps[2])
  
  plot_dat <- pca_res$scores[idx,] %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(
      cluster = clusts
    ) 
  
  plot_load <- pca_res$loadings[, c(comps_1, comps_2)] 
  plot_load <- plot_load[selectFeat(plot_load, nfeat = 5)$selectedFeat, ] %>%
    as.matrix() %>%
    as.data.frame()
  
  plot_dat %>%
    ggplot(
      aes(
        x = !!sym(comps_1),
        y = !!sym(comps_2)
      )
    ) +
    geom_point(
      aes(
        colour = cluster
      ),
      size = 1
    ) +
    scale_colour_manual(values = clust_cols) +
    geom_text(
      data = plot_load,
      label = rownames(plot_load),
      size = 2,
      position = position_jitter(0.1, 0.1)
    )
}



tempSpatPlot <- function(cType){
  
  p <- plot_dat %>%
    arrange(!!sym(cType)) %>%
    ggplot(
      aes(x, y)
    ) +
    geom_point(
      aes(
        colour = !!sym(cType)
      ),
      size = 0.8, stroke = 0
    ) +
    scale_colour_gradientn(colours = MATLAB_cols) +
    theme_void() +
    # ggtitle(cType) +
    theme(
      legend.position = "none"
    )
  
  return(p)
}


# Save PhiSpace annotation results as spatial heatmaps
tempSavePlots <- function(
    sce,
    methodName = "PhiSpace",
    tissueName = "",
    coordNames = c("sdimx", "sdimy"),
    freeColScale = F,
    quants = c(0.1,1),
    psize = 0.5,
    savePlots = TRUE,
    legendPosition = "none",
    censQuant = 1,
    returnPlot = F,
    # Save plots parameters
    width = 10,
    height = 10,
    fignrow = 4,
    figncol = 4,
    fsize_title = 10,
    plot_margin = c(0,0,0,0)
){
  
  scores <- as.matrix(
    reducedDim(sce, methodName)
  )
  
  if(freeColScale){
    
    lmts <- NULL
  } else {
    
    lmts <- quantile(scores, quants)
  }
  
  ctypes <- colnames(scores) %>% sort()
  outPlots <- vector("list", length(ctypes)) 
  names(outPlots) <- ctypes
  for(i in 1:length(ctypes)){
    ctype <- ctypes[i]
    cols <- scores[, ctype]
    cols <- PhiSpace:::censor(
      cols,
      quant = censQuant
    )
    outPlots[[i]] <-
      sce@colData %>%
      as.data.frame() %>%
      mutate(cols = cols) %>%
      arrange(cols) %>%
      ggplot(aes(x = !!sym(coordNames[1]), y = !!sym(coordNames[2]))) +
      geom_point(aes(colour = cols), size = psize, alpha = 1, stroke = 0) +
      theme_void() +
      scale_colour_gradientn(
        colours = MATLAB_cols,
        limits = lmts 
      ) +
      ggtitle(ctype) +
      theme(
        legend.position = legendPosition, 
        plot.title = element_text(size = fsize_title, hjust = 0.5, vjust = -1),
        plot.margin = unit(plot_margin, "points")
      )
  }
  if(savePlots){
    
    npheno <- length(outPlots)
    nFigPerPlot <- fignrow * figncol
    for(ii in 1:ceiling(npheno/nFigPerPlot)){
      
      idx_start <- (ii-1)*nFigPerPlot + 1
      idx_end <- min(ii*nFigPerPlot, npheno)
      ggsave(
        paste0(
          "figs/spatialHeatmap/", tissueName, "/", 
          tissueName, "_", methodName, "_", ii, ".png"
        ),
        ggpubr::ggarrange(
          plotlist = outPlots[idx_start:idx_end], 
          ncol = figncol, nrow = fignrow
        ), 
        width = width, height = height
      )
    }
  }
  if(returnPlot) return(outPlots)
}




tempMDSplot <- function(pca_res){
  
  mdsRes <- pca_res$scores  %>% `colnames<-`(c("PC1", "PC2")) %>% as.data.frame()
  mdsRes <- mdsRes %>% 
    mutate(
      cluster = rownames(mdsRes)
    )
  p <- ggplot(
    mdsRes,
    aes(PC1, PC2, color = cluster)
  ) + 
    geom_text(
      aes(label = cluster), size.unit = "pt", size = 8, key_glyph = "point", fontface = "bold"
    ) +
    scale_colour_manual(values = clust_cols) +
    theme_bw(base_size = 6) + 
    theme(
      legend.title = element_blank(),
      legend.position = "top", 
      legend.key.spacing = unit(-5, "pt")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 2)
      )
    )  
  return(p)
}




tempLoadBarplot <- function(
    pca_res, comp = "comp1", nfeat = 30, fsize = 6, xlab = "", rankABS = T
){
  
  if(rankABS){
    plot_dat <- pca_res$loadings %>% as.data.frame() %>% 
      arrange(desc(abs(!!sym(comp))))  %>%
      slice_head(n = nfeat) %>%
      arrange((!!sym(comp)))
  } else {
    plot_dat <- pca_res$loadings %>% as.data.frame() %>% 
      arrange(desc(!!sym(comp)))  %>%
      slice_head(n = nfeat) %>%
      arrange((!!sym(comp)))
  }
  rname_list <- lapply(
    strsplit(rownames(plot_dat), "<->"),
    function(x) sort(x)
  )
  rname_new <- sapply(
    rname_list,
    function(x) paste(x[1], x[2], sep = "<->")
  )
  rownames(plot_dat) <- rname_new
  p <- plot_dat %>%
    mutate(
      interaction = factor(
        rownames(.),
        levels = rownames(.)
      )
    ) %>% 
    ggplot() +
    geom_bar(
      aes(x = !!sym(comp), y = interaction, fill = !!sym(comp)),
      stat = "identity"
    ) +
    xlab(xlab) +
    scale_fill_gradientn(colours = c("blue", "white", "red")) +
    theme_pubr(base_size = fsize) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_blank(), 
      legend.position = "none"
    )
  return(p)
}

# For plotting kde contours
kde_plot_dat <- function(xx, cutoff = 0.1, bgridsize = 500){
  den_res <- ks::kde(xx, bgridsize = bgridsize)
  plot_dat <- data.frame(
    x = rep(den_res$eval.points[[1]], length(den_res$eval.points[[1]])),
    y = rep(den_res$eval.points[[1]], rep(length(den_res$eval.points[[1]]), length(den_res$eval.points[[1]])) ),
    z = as.vector(den_res$estimate)
  )
  out <- plot_dat[plot_dat$z > quantile(plot_dat$z, cutoff), ]
  return(out)
}