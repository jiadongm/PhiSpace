adjPlots <- function(
    p,
    theme_overall = c("void", "classic", "ggplot2", "bw", "linedraw", "light", "dark", "minimal"),
    # font sizes
    default.fsize = 6,
    axis.text.x.fsize = NULL,
    axis.text.y.fsize = NULL,
    axis.title.x.fsize = NULL,
    axis.title.y.fsize = NULL,
    legend.text.fsize = NULL,
    legend.title.fsize = NULL,
    # margins
    axis.title.x.margin = 0,
    axis.title.y.margin = 0,
    # positions
    legend.position = "top",
    legend.title.position = "left",
    # spacings
    legend.spacing.x = -8,
    legend.key.spacing = -5,
    legend.key.spacing.y = -10,
    legend.box.margin = -0,
    legend.box.spacing = 0,
    # sizes
    legend.obj.size = 2
){
  
  
  theme_overall <- match.arg(theme_overall)
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
  
  ## Font sizes
  if(is.null(axis.text.x.fsize)) axis.text.x.fsize <- default.fsize
  if(is.null(axis.text.y.fsize)) axis.text.y.fsize <- default.fsize
  if(is.null(axis.title.x.fsize)) axis.title.x.fsize <- default.fsize
  if(is.null(axis.title.y.fsize)) axis.title.y.fsize <- default.fsize
  if(is.null(legend.text.fsize)) legend.text.fsize <- default.fsize
  if(is.null(legend.title.fsize)) legend.title.fsize <- default.fsize
  
  p_adj <- p +
    theme_fun() +
    theme(
      # axis.text = element_text(size = axis.text.x.fsize),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.text = element_text(size = legend.text.fsize),
      legend.title = element_text(size = legend.title.fsize, face = "bold"),
      legend.position = legend.position,
      legend.title.position = legend.title.position,
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
      legend.background = element_blank()
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = legend.obj.size
        )
      ),
      shape = guide_legend(
        override.aes = list(
          size = legend.obj.size
        )
      )
    )
  
  return(p_adj)
  
}


alignmentScore <-
  function(
    data, 
    batch, 
    ncomp = NULL, 
    scale = F, 
    k = round(0.1 * nrow(data))
  ){
    if(!is.null(ncomp)){
      
      data <- getPC(X = data, ncomp = ncomp, scale = scale)$scores
    } 
    
    
    knn_mat <- FNN::get.knn(data, k=k)$nn.index
    
    x <- sapply(
      1:nrow(data),
      function(i){
        sum(
          batch[knn_mat[i,]] == batch[i]
        )
      }
    )
    
    result <- 1 - (mean(x) - k/nrow(data))/(k - k/nrow(data))
    return(result)
  }


tempPlotUMAP <- function(umap_res, 
                         idx, 
                         cols = "disease",
                         ylimts = NULL){
  
  plot_dat <- umap_res$layout %>%
    as.data.frame() %>%
    `colnames<-`(c("umap1", "umap2")) %>% 
    mutate(
      celltype = majCellType,
      disease = disease,
      batch = query$Batch[idx]
    )
  p1 <- plot_dat %>%
    ggplot(
      aes(
        umap1, umap2
      )
    ) +
    geom_point(
      aes(
        colour = !!sym(cols)
      ),
      size = 0.3,
      shape = 16,
      stroke = 0
    ) 
  
  if(cols == "celltype"){
    p1 <- p1 + scale_colour_manual(values = cols_type) 
  }
  
  if(cols == "disease"){
    p1 <- p1 + scale_colour_manual(values = cols_disease) 
  }
  
  if(!is.null(ylimts)){
    p1 <- p1 + ylim(ylimts)
  }
  
  return(p1)
  
}

mk_sc <- function(
    sc,
    tNames,
    dis_mat
){
  
  # sc <- cbind(
  #   sc[,tNames],
  #   sc[,sev_lvls]
  # ) %>%
  #   `colnames<-`(
  #     c(tNames, sev_lvls_simp)
  # )
  
  # if(is.null(dis_mat)){
  #   corr <- cor(
  #     sc[,tNames],
  #     sc[,sev_lvls_simp]
  #   )
  # } else {
  #   corr <- cor(
  #     sc[,tNames],
  #     dis_mat
  #   )
  # }
  
  tNames <- intersect(
    tNames,
    colnames(sc)
  )
  
  corr <- as.data.frame(
    sc[,tNames]
  ) %>%
    mutate(
      cond = as.character(dis_mat)
    ) %>%
    group_by(cond) %>%
    summarise(
      across(
        everything(),
        median
      )
    ) %>%
    data.frame(row.names = 1) %>%
    t
  
  
  return(corr)
}


tempCorHeat <- function(
    corMat, 
    row_title = "RNA derived",
    clust_rows = NULL
    # clust_cols = NULL
){
  
  if(is.null(clust_rows)){
    clust_rows <- as.dendrogram(seriate(dist(corMat), method = "GW")[[1]]) 
  } 
  
  # if(is.null(clust_cols)){
  #   clust_cols <- as.dendrogram(seriate(dist(t(corMat)), method = "GW")[[1]])
  # }
  
  p <- Heatmap(
    corMat,
    # col = aurora_cols,
    cluster_rows = clust_rows,
    cluster_columns = F,
    show_row_names = TRUE,
    show_column_names = TRUE,
    name = "Median score",
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_row_dend = F,
    show_column_dend = T,
    column_names_rot = 90,
    row_title = row_title, 
    row_title_gp = gpar(fontsize = 6),
    heatmap_legend_param = list(
      title_position = "leftcenter",
      title_gp = gpar(fontsize = 6, fontface = "bold"),
      grid_height = unit(2, "mm"),
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = 5),
      legend_direction = "horizontal"
    )
  )
  return(
    list(
      plot = p,
      clust_rows = clust_rows
    )
  )
}



adjPlots2 <- function(
    p,
    theme_overall = c("void", "classic", "ggplot2", "bw", "linedraw", "light", "dark", "minimal"),
    # font sizes
    default.fsize = 6,
    axis.text.x.fsize = NULL,
    axis.text.y.fsize = NULL,
    axis.title.x.fsize = NULL,
    axis.title.y.fsize = NULL,
    legend.text.fsize = NULL,
    legend.title.fsize = NULL,
    # margins
    axis.title.x.margin = 0,
    axis.title.y.margin = 0,
    # positions
    legend.position = "top",
    legend.title.position = "left",
    # spacings
    legend.spacing.x = 0,
    legend.key.spacing = 0,
    legend.key.spacing.y = -10,
    legend.box.margin = -0,
    legend.box.spacing = 0,
    # sizes
    legend.obj.size = 1
){
  
  
  theme_overall <- match.arg(theme_overall)
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
  
  ## Font sizes
  if(is.null(axis.text.x.fsize)) axis.text.x.fsize <- default.fsize
  if(is.null(axis.text.y.fsize)) axis.text.y.fsize <- default.fsize
  if(is.null(axis.title.x.fsize)) axis.title.x.fsize <- default.fsize
  if(is.null(axis.title.y.fsize)) axis.title.y.fsize <- default.fsize
  if(is.null(legend.text.fsize)) legend.text.fsize <- default.fsize
  if(is.null(legend.title.fsize)) legend.title.fsize <- default.fsize
  
  p_adj <- p +
    theme_fun() +
    theme(
      axis.text.y = element_text(size = axis.text.x.fsize),
      axis.text.x = element_text(size = axis.text.x.fsize, angle = 90),
      axis.title = element_blank(),
      legend.text = element_text(size = legend.text.fsize),
      # legend.title = element_text(size = legend.title.fsize, face = "bold"),
      legend.title = element_blank(),
      legend.position = legend.position,
      legend.title.position = legend.title.position,
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
      strip.text = element_text(
        size = default.fsize
      )
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = legend.obj.size
        )
      ),
      shape = guide_legend(
        override.aes = list(
          size = legend.obj.size
        )
      )
    )
  
  return(p_adj)
  
}
plotBmat <- function(Bhat){
  
  diseaseConds <- c(
    "Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical"
  )
  diseaseConds <- c(
    paste0(diseaseConds, "(ADT)"),
    paste0(diseaseConds, "(RNA)")
  )
  
  diseaseConds_simp <- c(
    "Healthy", "Asymp", "Mild", "Moderate", "Severe", "Critical"
  )
  diseaseConds_simp <- factor(
    c(
      diseaseConds_simp,
      diseaseConds_simp
    ),
    levels = diseaseConds_simp
  )
  
  modality <- rep(
    c("ADT", "RNA"), c(6, 6)
  )
  
  pdata <- as.data.frame(Bhat[diseaseConds,]) %>%
    mutate(
      phenotype = diseaseConds_simp,
      modality = modality
    ) %>%
    pivot_longer(
      !c(phenotype, modality),
      names_to = "disease",
      values_to = "coef"
    )
  
  
  p <- pdata %>%
    ggplot(
      aes(phenotype, coef, fill = modality)
    ) +
    geom_bar(
      stat = "identity"
    ) +
    facet_wrap(~disease, ncol = 2) 
  
  return(p)
}