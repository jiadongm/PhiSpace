## 10_analyse_performance.R
## Compare PhiSpace (normalised & unnormalised), RCTD, Cell2location, TACCO
## on individual datasets.
##
## Outputs:
##   figures/performance_by_dataset.pdf/png  — line chart per metric
##   figures/performance_delta.pdf/png       — delta from best method

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(scales)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  base_dir <- dirname(normalizePath(script_path))
} else {
  base_dir <- getwd()
}

fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- paste0("dataset", 1:32)

# ============================================================================
# 1. Load dataset sizes
# ============================================================================
cat("=== Loading dataset sizes ===\n")
sizes <- read.csv(file.path(fig_dir, "dataset_sizes.csv"),
                  stringsAsFactors = FALSE)
# Ensure column name consistency
if (!"dataset" %in% colnames(sizes)) colnames(sizes)[1] <- "dataset"
sizes$dataset <- as.integer(sizes$dataset)

# ============================================================================
# 2. Load per-dataset mean metrics (same loop for all four methods)
# ============================================================================
cat("=== Loading per-dataset metrics ===\n")

method_paths <- list(
  PhiSpace        = file.path(base_dir, "results_scran_log1p_filtered",
                               "all_metrics", datasets, "summary.csv"),
  PhiSpace_unnorm = file.path(base_dir, "results_scran_log1p_filtered",
                               "all_metrics_unnorm", datasets, "summary.csv"),
  PhiSpace_minmax = file.path(base_dir, "results_scran_log1p_filtered",
                               "all_metrics_minmax", datasets, "summary.csv"),
  RCTD            = file.path(base_dir, "results_rctd",
                               datasets, "metrics", "summary.csv"),
  Cell2location   = file.path(base_dir, "results_cell2location",
                               "all_metrics", datasets, "summary.csv"),
  TACCO           = file.path(base_dir, "results_tacco",
                               "all_metrics", datasets, "summary.csv")
)

load_method <- function(paths, method_name) {
  rows <- list()
  for (i in seq_along(paths)) {
    f <- paths[i]
    if (!file.exists(f)) next
    s <- read.csv(f, stringsAsFactors = FALSE)
    ds_num <- as.integer(sub("dataset", "", datasets[i]))
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds_num,
      Method  = method_name,
      PCC     = s$Mean[s$Metric == "PCC"],
      SSIM    = s$Mean[s$Metric == "SSIM"],
      RMSE    = s$Mean[s$Metric == "RMSE"],
      JSD     = s$Mean[s$Metric == "JSD"],
      stringsAsFactors = FALSE
    )
  }
  rbindlist(rows)
}

all_methods <- rbindlist(lapply(names(method_paths), function(m) {
  cat("  Loading", m, "...\n")
  load_method(method_paths[[m]], m)
}))

cat("  Loaded", nrow(all_methods), "dataset-method rows\n")

# Composite AS
all_methods[, AS := (PCC + SSIM + (1 - RMSE) + (1 - JSD)) / 4]

# Pivot to long format
df <- all_methods %>%
  pivot_longer(cols = c("PCC", "SSIM", "RMSE", "JSD", "AS"),
               names_to = "Metric", values_to = "Value") %>%
  as.data.table()

# ============================================================================
# 3. Join dataset sizes
# ============================================================================
df <- merge(df, sizes[, c("dataset", "n_celltypes", "n_cells",
                           "n_common_genes")],
            by = "dataset", all.x = TRUE)

# Order datasets by n_celltypes
ds_order <- sizes %>% arrange(n_celltypes) %>% pull(dataset)
df$dataset <- factor(df$dataset, levels = ds_order)

df$Metric <- factor(df$Metric, levels = c("PCC", "SSIM", "RMSE", "JSD", "AS"))
df$Method <- factor(df$Method,
                    levels = c("PhiSpace", "PhiSpace_unnorm",
                               "PhiSpace_minmax",
                               "RCTD", "Cell2location", "TACCO"))

# ============================================================================
# 4. Line chart — performance by dataset
# ============================================================================
cat("=== Plotting performance by dataset ===\n")

method_colours <- c(
  PhiSpace        = "#e74c3c",
  PhiSpace_unnorm = "#e67e22",
  PhiSpace_minmax = "#f39c12",
  RCTD            = "#2ecc71",
  Cell2location   = "#9b59b6",
  TACCO           = "#3498db"
)

p1 <- ggplot(df, aes(x = dataset, y = Value, group = Method)) +
  geom_line(aes(colour = Method), linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(fill = n_cells, size = n_common_genes, shape = Method),
             alpha = 0.8, stroke = 0.3) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_shape_manual(values = c(PhiSpace = 21, PhiSpace_unnorm = 25,
                                PhiSpace_minmax = 4,
                                RCTD = 22, Cell2location = 23,
                                TACCO = 24)) +
  scale_fill_viridis_c(name = "n_cells", option = "plasma") +
  scale_size_continuous(name = "n_common_genes", range = c(1, 4)) +
  scale_colour_manual(values = method_colours, guide = "none") +
  guides(
    fill  = guide_colourbar(order = 2),
    size  = guide_legend(order = 3),
    shape = guide_legend(order = 1, override.aes = list(size = 3))
  ) +
  labs(x = "Dataset (ordered by number of cell types)",
       y = "Score") +
  theme_bw(base_size = 11) +
  theme(
    strip.text         = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      size = 6),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "performance_by_dataset.pdf"), p1,
       width = 18, height = 5)
ggsave(file.path(fig_dir, "performance_by_dataset.png"), p1,
       width = 18, height = 5, dpi = 300)
cat("  Saved performance_by_dataset.pdf/png\n")

# ============================================================================
# 5. Delta-from-best panel
# ============================================================================
cat("=== Plotting delta from best ===\n")

# For PCC, SSIM, AS: best = max; for RMSE, JSD: best = min
higher_better <- c("PCC", "SSIM", "AS")

delta <- df %>%
  group_by(dataset, Metric) %>%
  mutate(
    best = ifelse(Metric %in% higher_better, max(Value, na.rm = TRUE),
                  min(Value, na.rm = TRUE)),
    delta = ifelse(Metric %in% higher_better, best - Value, Value - best)
  ) %>%
  ungroup()

# delta = 0 means this method is the best for that dataset+metric
p2 <- ggplot(delta, aes(x = dataset, y = delta, group = Method)) +
  geom_line(aes(colour = Method), linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(colour = Method, shape = Method), size = 1.5, alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = c(PhiSpace = 16, PhiSpace_unnorm = 8,
                                PhiSpace_minmax = 4,
                                RCTD = 17, Cell2location = 15,
                                TACCO = 18)) +
  labs(x = "Dataset (ordered by number of cell types)",
       y = "Delta from best (0 = best)") +
  theme_bw(base_size = 11) +
  theme(
    strip.text         = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      size = 6),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "performance_delta.pdf"), p2,
       width = 18, height = 5)
ggsave(file.path(fig_dir, "performance_delta.png"), p2,
       width = 18, height = 5, dpi = 300)
cat("  Saved performance_delta.pdf/png\n")

# ============================================================================
# 6. Print summary
# ============================================================================
cat("\n=== Times each method is best (per dataset+metric) ===\n")
best_counts <- delta %>%
  filter(delta == 0) %>%
  count(Method, Metric) %>%
  pivot_wider(names_from = Metric, values_from = n, values_fill = 0)
print(as.data.frame(best_counts))

cat("\n=== Mean delta from best (lower = closer to best) ===\n")
mean_delta <- delta %>%
  group_by(Method, Metric) %>%
  summarise(mean_delta = round(mean(delta, na.rm = TRUE), 4),
            .groups = "drop") %>%
  pivot_wider(names_from = Metric, values_from = mean_delta)
print(as.data.frame(mean_delta))

cat("\nDone. Figures saved to", fig_dir, "\n")
