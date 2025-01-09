xyRange <- 2000
binIDs <- unique(bcFilt$cell_id)
barAssay <- sapply(
  1:length(binIDs),
  function(x){
    
    binID <- binIDs[x]
    coords <- bcFilt %>%
      filter(cell_id == binID)
    coords <- coords[1, c("x_center", "y_center")]
    
    xcond <- (abs(bcFilt$x_center - coords[1,1]) < xyRange)
    ycond <- (abs(bcFilt$y_center - coords[1,2]) < xyRange)
    idx <- xcond & ycond
    
    table(bcFilt[idx,"barcode"])
    
  }
) %>%
  t() %>%
  `rownames<-`(binIDs)
plot_dat <- unique(
  bcFilt[,c("cell_id", xCoordName, yCoordName)]
) %>%
  as.data.frame()

clust_cols <- RColorBrewer::brewer.pal(10, "Set3")
names(clust_cols) <- as.character(1:10)

set.seed(523423)

outClusts <- lapply(
  2:10,
  function(kclust){
    
    clust_res <- kmeans(
      barAssay,
      centers = kclust,
      iter.max = 100,
      nstart = 50
    )
    
    return(clust_res)
  }
)
clust_list <- lapply(
  outClusts,
  function(x){
    factor(x$cluster, levels = sort(unique(x$cluster)))
  }
)

outPlots <- vector("list", length(clust_list))
for(x in length(clust_list):1){
  
  clust <- clust_list[[x]]
  
  if(x < length(clust_list)){
    clust_old <-  clust_list[[x+1]]
    sol <- clue::solve_LSAP(
      table(clust, clust_old),
      T
    )
    kclust <- length(unique(clust_old))
    adj <- (1:kclust)[sol]
    
    temp <- rep(NA, length(clust))
    
    for(jj in 1:length(adj)){
      
      lvls <- levels(clust)
      lvls_old <- levels(clust_old)
      temp[clust == lvls[jj]] <- lvls_old[adj[jj]]
    }
    
    clust <- as.factor(temp)
    clust_list[[x]] <- clust
  }
  
  plot_dat <- plot_dat %>%
    mutate(
      clusters = clust
    )
  
  outPlots[[x]] <- plot_dat %>%
    ggplot(
      aes(
        !!sym(xCoordName),
        !!sym(yCoordName),
        colour = clusters
      )
    ) +
    geom_point(
      size = 0.5
    ) +
    theme_void() +
    xlim(xylimits[,1]) +
    ylim(xylimits[,2]) +
    scale_colour_manual(values = clust_cols) +
    geom_polygon(
      data = ch,
      aes(x, y),
      fill = NA,
      alpha = 1,
      colour = "black"
    )
  
  
}

suppressWarnings(
  ggarrange(
    plotlist = outPlots,
    nrow = 3,
    ncol = 3
  )
)

ggsave(
  paste0("fig/clusts/spleen_clusts_xyRange=", xyRange, ".pdf"), 
  width = 12, 
  height = 12
)

plot_dat_wide <- barAssay %>%
  as.data.frame() %>%
  mutate(
    clust = clust_list[[7]]
  ) %>%
  group_by(clust) %>%
  summarise(
    across(
      starts_with("GFP"),
      sum
    )
  ) 
plot_prop <- plot_dat_wide %>%
  select(!clust) %>%
  as.matrix()
plot_prop <- (plot_prop/rowSums(plot_prop)) %>%
  as.data.frame() %>%
  mutate(
    clust = plot_dat_wide$clust
  ) %>%
  pivot_longer(
    !clust,
    values_to = "count",
    names_to = "barcode"
  ) %>%
  mutate(
    barcode = factor(
      barcode,
      levels = names(sortedBC)
    )
  ) 
plot_prop %>%
  ggplot(
    aes(barcode, count)
  ) +
  geom_bar(
    aes(fill = clust),
    stat = "identity"
  ) +
  scale_fill_manual(
    values = clust_cols
  ) +
  facet_wrap(~ clust) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(
  paste0("fig/clusts/clustComposition_xyRange=", xyRange, ".pdf"), 
  width = 12, 
  height = 12
)

# Raw
plot_dat_wideRaw <- barAssayRaw %>%
  as.data.frame() %>%
  mutate(
    clust = clust_list[[7]]
  ) %>%
  group_by(clust) %>%
  summarise(
    across(
      starts_with("GFP"),
      sum
    )
  ) 
plot_propRaw <- plot_dat_wideRaw %>%
  select(!clust) %>%
  as.matrix()
plot_propRaw <- (plot_propRaw/rowSums(plot_propRaw)) %>%
  as.data.frame() %>%
  mutate(
    clust = plot_dat_wideRaw$clust
  ) %>%
  pivot_longer(
    !clust,
    values_to = "count",
    names_to = "barcode"
  ) %>%
  mutate(
    barcode = factor(
      barcode,
      levels = names(sortedBC)
    )
  ) 
plot_propRaw %>%
  ggplot(
    aes(barcode, count)
  ) +
  geom_bar(
    aes(fill = clust),
    stat = "identity"
  ) +
  scale_fill_manual(
    values = clust_cols
  ) +
  facet_wrap(~ clust) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(
  paste0("fig/clusts/clustCompositionRaw_xyRange=", xyRange, ".pdf"), 
  width = 12, 
  height = 12
)