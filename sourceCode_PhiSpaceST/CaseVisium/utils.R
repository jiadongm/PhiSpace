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
  "Background", "LUAD" # P25
  # "Healthy", "Healthy", # D1
  # "Background", "LUAD", "LUAD", "LUAD", # P10
  # "Background", "LUSC", "LUSC", "LUSC", # P11
  # "Background", "LUAD", # P15
  # "Background", "LUAD", # P16
  # "Background", "LUSC", # P17
  # "Background", "LUSC", # P19
  # "Background", "LUAD", # P24
  # "Background", "LUAD" # P25
)

VisiumCancerCols <- c(
  Background = "#377EB8",
  Healthy = "#4DAF4A",
  LUAD = "#984EA3",
  LUSC = "#E41A1C"
)

tissueNames_CosMx <- c(
  "Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3",
  "Lung6",
  "Lung9_Rep1", "Lung9_Rep2",
  "Lung12",     "Lung13"
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

## Binning subcellular ST into square bins
# Based on SubcellularSpatialData::allocateSquare
squareBin <- function(x, y, bw = 50){
  
  stopifnot(length(x) == length(y))
  stopifnot(bw > 0)
  
  # bw = diff(range(x))/bins
  x = ceiling((x - min(x))/bw)
  y = ceiling((y - min(y))/bw)
  x[x == 0] = 1
  y[y == 0] = 1
  max_x = max(x)
  sqix = data.frame(
    bin_id = x + (y - 1) * max_x, 
    bin_x = x, 
    bin_y = y
  )
  return(sqix)
}

## Binning transcripts into sparse count matrix
# a simplified version of SubcellularSpatialData::tx2spe
binTranscripts <- function(x, bw = 50){
  
  required_cols = c("gene", "x", "y", "counts")
  missing_cols = setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop(paste0("The following columns are missing from the provided transcript table 'x': ", 
                paste(missing_cols, collapse = ", ")))
  }
  
  x = dplyr::ungroup(
    dplyr::group_by(
      cbind(x, squareBin(x$x, x$y, bw = bw)),
      bin_id, 
      .add = TRUE
    )
  )
  
  
  bin_annot = dplyr::ungroup(
    dplyr::summarise(
      dplyr::group_by(
        dplyr::select(
          dplyr::select(
            dplyr::filter(
              x, !is.na(bin_id)
            ), 
            !c(gene)
          ), 
          bin_id, 
          dplyr::where(is.character)
        ), 
        bin_id
      ), 
      dplyr::across(
        everything(), 
        ~ifelse(
          all(is.na(.x)), 
          NA_character_, 
          tail(names(sort(table(na.omit(.x)))),1)
        )
      )
    )
  )
  
  
  bin_annot = as.data.frame(
    dplyr::mutate(
      dplyr::full_join(
        dplyr::summarise(
          dplyr::group_by(
            dplyr::select(
              dplyr::filter(x,!is.na(bin_id)), 
              !c(gene, counts)
            ), 
            bin_id
          ), 
          dplyr::across(dplyr::where(is.numeric), ~mean(.x, na.rm = TRUE))
        ), 
        bin_annot, 
        by = dplyr::join_by(bin_id)
      ), 
      rname = bin_id
    )
  )
  
  rownames(bin_annot) = bin_annot$rname
  bin_annot = dplyr::select(bin_annot, !rname)
  
  gene_annot = dplyr::filter(
    dplyr::select(dplyr::filter(x, !is.na(bin_id)), gene), 
    !duplicated(paste(gene))
  )
  
  counts = with(
    dplyr::summarise(
      dplyr::group_by(
        dplyr::select(
          dplyr::mutate(
            dplyr::filter(
              x, 
              !is.na(bin_id)
              ), 
            rname = as.factor(bin_id), 
            gene = as.factor(gene)
          ), 
          c(rname, gene, counts)
        ), 
        rname, gene
      ), 
      counts = sum(counts)
    ), 
    Matrix::sparseMatrix(
      i = as.numeric(gene), 
      j = as.numeric(rname), 
      x = counts, 
      dimnames = list(levels(gene), levels(rname))
    )
  )
  
  counts = counts[gene_annot$gene, rownames(bin_annot)]
  
  sce = SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts), 
    colData = bin_annot
  )
  
  return(sce)
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


tempPCAplot <- function(pca_res){
  
  mdsRes <- pca_res$scores  %>% `colnames<-`(c("PC1", "PC2")) %>% as.data.frame()
  mdsLabs <- rep("Tumour-background", nrow(mdsRes))
  mdsLabs[grepl("T", rownames(mdsRes))] <- "Tumour"
  mdsLabs[grepl("D", rownames(mdsRes))] <- "Control"
  mdsRes <- mdsRes %>% 
    mutate(
      `Sample type` = mdsLabs,
      sample = rownames(mdsRes),
      cancerType = Visium_cancerTypes
    )
  p <- ggplot(
    mdsRes,
    aes(PC1, PC2, color = `cancerType`)
  ) + 
    geom_text(
      aes(label = sample), size.unit = "pt", size = 6, key_glyph = "point"
    ) +
    scale_colour_manual(values = VisiumCancerCols) +
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


# Venn diagram
tempVenn <- function(
    deg.results, 
    AzimuthMarkers,
    # PhiSpace = TRUE, 
    DEG_max = 100, 
    plotVenn = TRUE,
    title = "",
    fsize = 6
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
        "#80B1D3",
        "#FB8072",
        "#8DA0CB"
      )[1:length(plist)],
      stroke_size = 1,  
      show_percentage = F, 
      text_size = fsize, 
      set_name_size = 0
    ) +
      ggtitle(title) +
      theme(
        title = element_text(size = 6)
      )
    
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


# For calcualting prognosis scores
meanEsts <- function(mkNames){
  
  mkLvls <- sapply(
    1:18,
    function(x){
      
      query <- query_list[[x]]
      mean(rowSums(assay(query, assay2use)[mkNames, ,drop=F]))/meanGeneLvls[x]
    }
  )
  return(mkLvls)
}















