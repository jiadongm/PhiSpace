tissueNames_Visium <- c(
  "D2_1", "D2_2", 
  "P10_T1", "P10_B1", 
  "P11_T3", "P11_B1", 
  "P15_B1", "P15_T2", 
  "P16_B1", "P16_T1",
  "P17_B1", "P17_T1", 
  "P19_B2", "P19_T1", 
  "P24_B1", "P24_T1", 
  "P25_B2", "P25_T1",
  "D1_1", "D1_2",
  "P10_B2", "P10_T2", "P10_T3", "P10_T4",
  "P11_B2", "P11_T1", "P11_T2", "P11_T4",
  "P15_B2", "P15_T1",
  "P16_B2", "P16_T2",
  "P17_B2", "P17_T2",
  "P19_B1", "P19_T2",
  "P24_B2", "P24_T2",
  "P25_B1", "P25_T2"
)

Visium_cancerTypes <- c(
  "Healthy", "Healthy", # D2
  "LUAD", "Background", # P10
  "LUSC", "Background", # P11
  "Background", "LUAD", # P15
  "Background", "LUAD", # P16
  "Background", "LUSC", # P17
  "Background", "LUSC", # P19
  "Background", "LUAD", # P24
  "Background", "LUAD", # P25
  "Healthy", "Healthy", # D1
  "Background", "LUAD", "LUAD", "LUAD", # P10
  "Background", "LUSC", "LUSC", "LUSC", # P11
  "Background", "LUAD", # P15
  "Background", "LUAD", # P16
  "Background", "LUSC", # P17
  "Background", "LUSC", # P19
  "Background", "LUAD", # P24
  "Background", "LUAD" # P25
)

tissueNames_CosMx <- c(
  "Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3",
  "Lung6",
  "Lung9_Rep1", "Lung9_Rep2",
  "Lung12",     "Lung13"
)

VisiumCancerCols <- c(
  LUSC = "#E41A1C",
  LUAD = "#984EA3",
  Healthy = "#4DAF4A",
  Background = "#377EB8"
)

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

# Colour scheme for spatial niches defined by data generators
nicheCol <- c(
  `tumor interior` = "black",
  tumour = "black",
  `tumor-stroma boundary` = "#999999",
  immune = "#E41A1C",
  `myeloid-enriched stroma` = "#377EB8",
  neutrophils = "#FF7F00",
  macrophages = "#984EA3",
  stroma = "#4DAF4A",
  `plasmablast-enriched stroma` = "#F781BF",
  `lymphoid structure` = "#A65628"
)



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





# Venn diagram
tempVenn <- function(
    deg.results, 
    AzimuthMarkers,
    # PhiSpace = TRUE, 
    DEG_max = 100, 
    plotVenn = TRUE,
    title = ""
  ){
  
  # suppressWarnings(
  #   EnhancedVolcano::EnhancedVolcano(
  #     deg.results %>% as.data.frame(),
  #     lab = rownames(deg.results),
  #     x = "avg_log2FC",
  #     y = "p_val_adj",
  #     FCcutoff = 2, 
  #     pCutoff = 1e-100
  #   )
  # )
  
  # gene imporatnce is measured by Pi score of Xiao et al. (2014)
  if(is.null(deg.results$p_val_adj)){
    
    gImportance <- vector("numeric", length = 0)
  } else {
    
    gImportance <- -log10(deg.results$p_val_adj) * deg.results$avg_log2FC
  }
  
  deg_num <- min(
    length(gImportance),
    DEG_max
  )
  
  # Top DEG obtained by MAST method
  dgeSet <- rownames(deg.results)[rank(-gImportance) <= deg_num] %>% sort
  # PhiSpace marker genes
  # if(PhiSpace) mkGenes <- orderedFeat[1:100,cType]
  # Azimuth marker genes (ground truth)
  AzimuthGenes <- AzimuthMarkers[[cType]]
  
  # if(PhiSpace){
  #   
  #   plist <- list(
  #     SpatialDEG = dgeSet,
  #     AzimuthMarkers = AzimuthGenes,
  #     PhiSpaceMarkers = mkGenes
  #   )
  # } else {
  #   
  #   plist <- list(
  #     DEG_MAST = dgeSet,
  #     Azimuth = AzimuthGenes
  #   )
  # }
  
  plist <- list(
    DEG_MAST = dgeSet,
    Azimuth = AzimuthGenes
  )
  
  
  if(plotVenn){
    
    p <- ggvenn::ggvenn(
      plist,
      fill_color = c(
        "#66C2A5",
        "#FC8D62",
        "#8DA0CB"
      )[1:length(plist)],
      stroke_size = 1,  
      show_percentage = F, 
      text_size = 8, 
      set_name_size = 0
    ) +
      ggtitle(title)
  } else {
    
    p <- NULL
  }
  
  return(
    list(
      plot = p,
      intersection = length(intersect(dgeSet, AzimuthGenes)),
      dgeSet = dgeSet,
      AzimuthGenes = AzimuthGenes
    )
  )
}


fitZINB <- function(x){
  
  # Initial fit (without ZI)
  suppressWarnings(fitRes <- MASS::fitdistr(x, "negative binomial"))
  size <- fitRes$estimate["size"]
  prob <- size/sum(fitRes$estimate)
  # ZINB fit
  fitResZI <- iZID::nb.zihmle(
    x, r = size, p = prob, type = "zi"
  )
  mu <- unname(fitResZI[,"r"]*(1-fitResZI[,"p"])/fitResZI[,"p"])
  ZIprob <- unname(fitResZI[,"phi"])
  return(
    list(
      mu = mu,
      ZIprob = ZIprob
    )
  )
}



## Clustering purity
purity <- function(clust, origin){
  
  tab <- table(clust, origin)
  col_sums <- colSums(tab) # Sizes of referenec ('true') clusters
  row_sums <- rowSums(tab) # Sizes of predicted clusters
  tot_sum <- sum(col_sums) # total sample size
  
  col_maxes <- apply(tab, 2, max) 
  row_maxes <- apply(tab, 1, max)
  
  Pur <- sum(col_maxes)/tot_sum
  invPur <- sum(row_maxes)/tot_sum
  
  # Van Rijsbergen's F measure
  sumSizes <- 
    tcrossprod(
      rep(1, ncol(tab)),
      col_sums
    ) +
    tcrossprod(
      row_sums,
      rep(1, nrow(tab))
    )
  
  Fvalues <- 2 * tab/sumSizes
  
  Fmeasure <- sum(apply(Fvalues, 2, max) * col_sums) / tot_sum
  
  
  return(
    list(
      Pur = Pur,
      invPur = invPur,
      Fmeasure = Fmeasure
    )
  )
}




sparse.cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans
  
  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}




# For plotting PCs of vectorised correlation matrices
tempPCAplot <- function(pca_res, PCidx = c(1,2), fsize = 6){
  
  
  mdsRes <- pca_res$scores  %>% `colnames<-`(paste0("PC", 1:ncol(pca_res$scores))) %>% as.data.frame()
  mdsRes <- mdsRes %>%
    mutate(
      niche = corMatVec$niche,
      sample = gsub("Lung", "L", corMatVec$sample)
    )
  p <- ggplot(
    mdsRes,
    aes(
      !!sym(colnames(mdsRes)[PCidx[1]]), 
      !!sym(colnames(mdsRes)[PCidx[2]]), 
      colour = niche
    )
  ) + 
    geom_text(
      aes(label = sample), size.unit = "pt", size = fsize, key_glyph = "point"
    ) +
    scale_colour_manual(values = nicheCol) +
    theme_bw(base_size = fsize) + 
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



tempHeatmap <- function(
    interaction_matrix, 
    fsize = 6, 
    binColour = T,  
    useSeriate = T,
    Nbins = 10,
    preDefOrder = NULL
  ){
  
  if(binColour){
    
    brks <- quantile(as.vector(interaction_matrix), seq(0,1,length.out = Nbins))
    interaction_matrix_binned <- cut(
      interaction_matrix,
      breaks = brks,
      labels = as.character(brks[2:length(brks)])
    ) |>
      as.character() |>
      as.numeric() |>
      matrix(nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix)) |>
      `dimnames<-`(
        dimnames(interaction_matrix)
      )
  } else {
    
    interaction_matrix_binned <- interaction_matrix
  }
  
  if(useSeriate){
    
    o <- seriate(interaction_matrix)
    hm <- Heatmap(
      interaction_matrix_binned, 
      name = "Correlation",
      show_row_names = T, show_column_names = F,
      show_heatmap_legend = F,
      row_order = get_order(o, 1), column_order = get_order(o, 2), 
      row_names_gp = gpar(fontsize = fsize), column_names_gp = gpar(fontsize = fsize)
      # rect_gp = gpar(type = "none"),
      # cell_fun = function(j, i, x, y, w, h, fill) {
      #   if(as.numeric(x) <= 1 - as.numeric(y)) {
      #     grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
      #   }
      # }
    )
  } else {
    
    if(is.null(preDefOrder)) stop("preDefOrder cannot be NULL.")
    hm <- Heatmap(
      interaction_matrix_binned, show_column_names = T,
      show_heatmap_legend = F, show_row_dend = F, show_column_dend = F,
      row_order = preDefOrder, column_order = preDefOrder, 
      row_names_gp = gpar(fontsize = fsize), column_names_gp = gpar(fontsize = fsize)
    )
  }
  
  return(hm)
}





tempClustSignature <- function(
    sce, clustVarName, clustIdx, whichReduced = "PhiSpace",  
    topK = 10, byAbs = F, fsize = 6
    ){
  
  plotDat <- reducedDim(sce, whichReduced) %>% as.data.frame() %>%
    mutate(clusters = colData(sce)[,clustVarName]) %>%
    pivot_longer(!clusters, names_to = "cellType", values_to = "scores") %>%
    filter(clusters == clustIdx) %>%
    group_by(cellType) %>% mutate(MedianValue = median(scores)) %>% 
    ungroup() %>% mutate(cellType = reorder(cellType, MedianValue))
  
  # plotDat$broadType <- category_lookup[plotDat$cellType %>% as.character()]
  
  if(byAbs){
    topBoxes <- plotDat %>%
      distinct(cellType, MedianValue) %>%
      arrange(desc(abs(MedianValue))) %>%
      slice(1:topK) %>%
      pull(cellType)
  } else {
    topBoxes <- plotDat %>%
      distinct(cellType, MedianValue) %>%
      arrange(desc(MedianValue)) %>%
      slice(1:topK) %>%
      pull(cellType)
  }
  
  
  p <-  plotDat %>%
    filter(cellType %in% topBoxes) %>%
    ggplot(
      aes(
        cellType, scores
        , fill = cellType
      )
    ) + geom_boxplot(
      outlier.size = 0.3, outlier.stroke = 0, linewidth = 0.1
    ) +
    scale_fill_manual(values = tempCTypeCols) +
    theme_bw(base_size = fsize) +
    theme(
      axis.text.x = element_text(hjust = 1, angle = 90)
    ) +
    ylim(c(-0.5,1)) +
    ggtitle(
      paste("Cluster", clustIdx)
    )
  return(p)
}



