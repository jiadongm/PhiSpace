library(tidyr)

Errs_list <- readRDS("output/benchPB/PBMC_Errs_list.rds")
Errs <- do.call(cbind, Errs_list)
Errs <- Errs %>% 
  t() %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "method", values_to = "error") 
# Errs_over <- Errs %>%
#   filter(grepl("ver", method)) %>% 
#   mutate(method = factor(
#     method, levels = c("NP1nnOver", "NPover", "NPnormOver", "SRover") 
#   ))
Errs_bal <- Errs %>%
  filter(grepl("al", method))
# p1 <-
#   Errs_over %>%
#   ggplot(aes(x = method, y = error)) +
#   geom_boxplot()
p2 <-
  Errs_bal %>%
  ggplot(aes(x = error, y = method)) +
  geom_boxplot()
p2


Errs_list <- readRDS("output/benchPB/BM_Errs_list.rds")
Errs <- do.call(cbind, Errs_list)
Errs <- Errs %>% 
  t() %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "method", values_to = "error") 
# Errs_over <- Errs %>%
#   filter(grepl("ver", method)) %>% 
#   mutate(method = factor(
#     method, levels = c("NP1nnOver", "NPover", "NPnormOver", "SRover") 
#   ))
Errs_bal <- Errs %>%
  filter(grepl("al", method))
# p1 <-
#   Errs_over %>%
#   ggplot(aes(x = method, y = error)) +
#   geom_boxplot()
p2 <-
  Errs_bal %>%
  ggplot(aes(x = error, y = method)) +
  geom_boxplot()
p2


### Bridge integration results -------------------------------------------------

## PBMC
Bint_errs <- readRDS("output/benchBridge/PBMC_BridgeInt_Errs.rds")
Bint_errs <- do.call(cbind, Bint_errs)
NP_errs <- readRDS("output/benchBridge/PBMC_NPints_Errs.rds")
NP_errs <- do.call(cbind, NP_errs)
rbind(Bint_errs, NP_errs) %>%
  t() %>%
  as.data.frame() %>%
  select(ends_with("bal", ignore.case = T)) %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "err") %>%
  ggplot(aes(x= err, y = method)) +
  geom_boxplot()


## Multi-batch BMMC
NP_errs <- readRDS("output/benchBridge/BM_NPints_Errs.rds")
Bint_errs <- readRDS("output/benchBridge/BM_BridgeInt_Errs.rds")
NP_errs[[1]] %>%
  t() %>%
  as.data.frame() %>%
  select(ends_with("bal", ignore.case = T)) %>%
  mutate(BridgeIntBal = Bint_errs[[1]][2,], .before = SingleRbal) %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "err") %>%
  ggplot(aes(x= err, y = method)) +
  geom_boxplot()

NP_errs[[2]] %>%
  t() %>%
  as.data.frame() %>%
  select(ends_with("bal", ignore.case = T)) %>%
  mutate(BridgeIntBal = Bint_errs[[2]][2,], .before = SingleRbal) %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "err") %>%
  ggplot(aes(x= err, y = method)) +
  geom_boxplot()









if(F){
  ### Old ------------------------------------------------------------------------
  library(umap)
  library(plotly)
  setwd("~/Dropbox/Research projects/Sincast2.0/Theme_bridge")
  Errs_list <- readRDS("output/benchPB/PBMC_Errs_list_balancedAgg.rds")
  Errs_pbmc <- do.call(cbind, Errs_list)
  Errs_list <- readRDS("output/benchPB/BM_Errs_list_balancedAgg.rds")
  Errs_bm <- do.call(cbind, Errs_list)
  # Errs_list <- readRDS("output/benchPB/multi_Errs_list.rds")
  # Errs_bm <- do.call(cbind, Errs_list)
  # boxplot(t(Errs[1:2,]))
  # boxplot(Errs[1,] - Errs[2,])
  
  
  
  Errs <- Errs_bm %>% 
    t() %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = "method", values_to = "error") 
  Errs_over <- Errs %>%
    filter(grepl("ver", method)) %>% 
    mutate(method = factor(
      method, levels = c("NP1nnOver", "NPover", "NPnormOver", "SRover") 
    ))
  Errs_bal <- Errs %>%
    filter(grepl("al", method))
  p1 <-
    Errs_over %>%
    ggplot(aes(x = method, y = error)) +
    geom_boxplot()
  p2 <-
    Errs_bal %>%
    ggplot(aes(x = method, y = error)) +
    geom_boxplot()
  p1+p2
  
  
  ## plotly to render violin plot
  Errs <- data.frame(ers = c(as.vector(t(Errs_pbmc[1:2,])), as.vector(t(Errs_bm))),
                     labels_fine = rep(c("SingleR PBMC", "NPints PBMC", "SingleR BM", "NPints BM"), rep(100, 4)),
                     labels_dataset = rep(c("PBMC", "BM"), c(200,200)),
                     labels_method = rep(c("SingleR", "NPints", "SingleR", "NPints"), rep(100, 4))
  )
  fig <- Errs %>%
    plot_ly(type = "violin") %>%
    add_trace(x = ~labels_dataset[labels_method == "SingleR"],
              y = ~ers[labels_method == "SingleR"],
              legendgroup = "SingleR",
              scalegroup = "SingleR",
              name = "SingleR",
              side = "negative",
              box = list(visible = T)) %>%
    add_trace(x = ~labels_dataset[labels_method == "NPints"],
              y = ~ers[labels_method == "NPints"],
              legendgroup = "NPints",
              scalegroup = "NPints",
              name = "NPints",
              side = "positive",
              box = list(visible = T)) %>%
    layout(xaxis = list(title = ""),
           yaxis = list(title = "Overal error"))
  fig  
  
  
  wilcox.test(Errs_pbmc[1,], Errs_pbmc[2,], "greater")
  wilcox.test(Errs_bm[1,], Errs_bm[2,], "greater")
  
  fig %>% add_annotations(x = c("PBMC", "BM"), y = c("0.18", "0.14"),
                          text = c("P-value = 0.1498", "P-value = 0"),
                          showarrow = FALSE,
                          font = list(size = 24))
  
  
  
  
  
  ## Visualising soft classification results
  custom.settings <- umap.defaults
  custom.settings$metric <- "pearson2"
  # custom.settings$n_neighbors <- 15
  # custom.settings$n_components <- 3
  # toUMAP <- getPC(t(assay(query, "rank")), ncomp = 30)$scores
  toUMAP <- NPintsScore
  # toUMAP <- sr_re$scores
  umap.re <- umap(toUMAP, config = custom.settings)
  fig <- plotly::plot_ly() %>%
    plotly::add_markers(data = data.frame(umap.re$layout),
                        x = ~X1, 
                        y = ~X2, 
                        # z = ~X3,
                        color = colData(query)[, YtrainName],
                        # colors = referenceColors,
                        marker = list(size = 5)) %>%
    plotly::layout(legend = list(orientation = "h"))
  fig
}





