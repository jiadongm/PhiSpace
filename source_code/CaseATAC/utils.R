# Global viz parameters
fsize <- 6
# margins
axis.title.x.margin = 0
axis.title.y.margin = 0
# positions
legend.position = "top"
legend.title.position = "left"
# spacings
legend.spacing.x = -2
legend.key.spacing = -0
legend.key.spacing.y = -10
legend.box.margin = -0
legend.box.spacing = 0
# sizes
legend.obj.size = 2




tempPlot <- function(
    umap_res,
    labs,
    labName = "celltype",
    colSc = "Set1",
    ## Theme
    theme_overall = c("void", "classic", "ggplot2", "bw", "linedraw", "light", "dark", "minimal"),
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
    legend.obj.size = 2,
    seed = 2039
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
  
  
  plot_dat <- umap_res %>%
    mutate(labs = labs) 
  set.seed(seed)
  idx <- sample(1:nrow(plot_dat))
  plot_dat <- plot_dat[idx,]
  
  plot_dat %>%
    ggplot(aes(comp1, comp2)) +
    geom_point(aes(colour = labs), size = 0.7, stroke = 0) +
    scale_colour_brewer(
      palette = colSc
    ) +
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
        ),
        title = labName
      ),
      shape = guide_legend(
        override.aes = list(
          size = legend.obj.size
        )
      )
    )
}



# Definition of purity, for evaluating clustering.
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