## 09_visualise_results.R
## 1. Runtime boxplots (PhiSpace, TACCO, Seurat, SPOTlight, RCTD, Cell2location)
## 2. Performance boxplots across 7 methods (PCC, SSIM, RMSE, JSD, AS)
##
## Methods included (7): PhiSpace, RCTD, Cell2location, TACCO, Seurat, SPOTlight, DSTG
## Re-run results used for: PhiSpace, RCTD, Cell2location, TACCO, Seurat, SPOTlight
## Paper values used only for: DSTG (no re-run available)
## Runtime plotted for 6 methods (excluding DSTG — no runtime data).
##
## PhiSpace uses scran+log1p with cell-type filtering (<25 cells).

library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

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

# ===========================================================================
# 1. Runtime boxplot (6 methods — all except DSTG)
# ===========================================================================
cat("=== Collecting runtime data ===\n")

collect_runtime <- function() {
  rows <- list()
  for (ds in datasets) {
    ds_num <- as.integer(sub("dataset", "", ds))

    # PhiSpace (scran + log1p filtered)
    phi_file <- file.path(base_dir, "results_scran_log1p_filtered", ds,
                          "PhiSpace_runtime.txt")
    if (file.exists(phi_file)) {
      phi <- read.csv(phi_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "PhiSpace",
        runtime_seconds = phi$runtime_seconds
      )
    }

    # TACCO
    tacco_file <- file.path(base_dir, "results_tacco", ds,
                            "TACCO_runtime.txt")
    if (file.exists(tacco_file)) {
      tacco <- read.csv(tacco_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "TACCO",
        runtime_seconds = tacco$runtime_seconds
      )
    }

    # Seurat
    seurat_file <- file.path(base_dir, "results_seurat", ds,
                             "Seurat_runtime.txt")
    if (file.exists(seurat_file)) {
      seurat <- read.csv(seurat_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "Seurat",
        runtime_seconds = seurat$runtime_seconds
      )
    }

    # SPOTlight
    spotlight_file <- file.path(base_dir, "results_spotlight", ds,
                                "SPOTlight_runtime.txt")
    if (file.exists(spotlight_file)) {
      spotlight <- read.csv(spotlight_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "SPOTlight",
        runtime_seconds = spotlight$runtime_seconds
      )
    }

    # RCTD
    rctd_file <- file.path(base_dir, "results_rctd", ds,
                           "RCTD_runtime.txt")
    if (file.exists(rctd_file)) {
      rctd <- read.csv(rctd_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "RCTD",
        runtime_seconds = rctd$runtime_seconds
      )
    }

    # Cell2location
    c2l_file <- file.path(base_dir, "results_cell2location", ds,
                          "Cell2location_runtime.txt")
    if (file.exists(c2l_file)) {
      c2l <- read.csv(c2l_file)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds_num, method = "Cell2location",
        runtime_seconds = c2l$runtime_seconds
      )
    }
  }
  rbindlist(rows)
}

runtime_df <- collect_runtime()

cat("Methods found:", paste(unique(runtime_df$method), collapse = ", "), "\n")
cat("Datasets per method:\n")
print(runtime_df[, .N, by = method])

if (nrow(runtime_df) > 0) {
  method_order <- c("PhiSpace", "TACCO", "Seurat", "SPOTlight",
                    "RCTD", "Cell2location")
  method_order <- method_order[method_order %in% unique(runtime_df$method)]
  runtime_df$method <- factor(runtime_df$method, levels = method_order)

  method_counts <- runtime_df[, .N, by = method]
  method_labels <- setNames(
    paste0(method_counts$method, "\n(n=", method_counts$N, ")"),
    method_counts$method
  )

  runtime_cols <- c(
    "PhiSpace"      = "#e74c3c",
    "TACCO"         = "#3498db",
    "Seurat"        = "#f39c12",
    "SPOTlight"     = "#1abc9c",
    "RCTD"          = "#2ecc71",
    "Cell2location" = "#9b59b6"
  )

  p_rt <- ggplot(runtime_df,
                 aes(x = method, y = runtime_seconds, fill = method)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, width = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
    scale_x_discrete(labels = method_labels) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    scale_fill_manual(values = runtime_cols) +
    labs(x = NULL, y = "Runtime (seconds, log scale)",
         title = "Computation time across 32 simulated datasets") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "runtime_boxplot.pdf"), p_rt,
         width = 8, height = 5)
  ggsave(file.path(fig_dir, "runtime_boxplot.png"), p_rt,
         width = 8, height = 5, dpi = 300)
  cat("Runtime boxplot saved.\n")

  cat("\n=== Runtime summary (seconds) ===\n")
  print(runtime_df[, .(
    n      = .N,
    mean   = round(mean(runtime_seconds), 1),
    median = round(median(runtime_seconds), 1),
    sd     = round(sd(runtime_seconds), 1),
    min    = round(min(runtime_seconds), 1),
    max    = round(max(runtime_seconds), 1)
  ), by = method])
}

# ===========================================================================
# 2. Performance boxplots (7 methods)
# ===========================================================================
cat("\n=== Loading performance metrics ===\n")

# --- Li et al. paper methods (per-cell-type CSVs) -------------------------
# Only used for DSTG (not re-run by us)
paper_dir <- file.path(dirname(base_dir),
                       "FigureData", "Figure4", "32SimulationData")

load_paper_metric <- function(file, metric_name) {
  d <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  methods <- setdiff(colnames(d), c("", "dataset"))
  colnames(d)[1] <- "celltype"
  d %>%
    select(-celltype) %>%
    group_by(dataset) %>%
    summarise(across(all_of(methods), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    pivot_longer(-dataset, names_to = "Method", values_to = "Value") %>%
    mutate(Metric = metric_name)
}

paper_pcc  <- load_paper_metric(file.path(paper_dir, "PCC.csv"),  "PCC")
paper_ssim <- load_paper_metric(file.path(paper_dir, "SSIM.csv"), "SSIM")
paper_rmse <- load_paper_metric(file.path(paper_dir, "RMSE.csv"), "RMSE")
paper_jsd  <- load_paper_metric(file.path(paper_dir, "JS.csv"),   "JSD")

paper <- bind_rows(paper_pcc, paper_ssim, paper_rmse, paper_jsd)

# Keep only DSTG from paper values (all other methods are re-run by us)
paper <- paper %>% filter(Method == "DSTG")
paper$dataset <- as.integer(paper$dataset)

# --- Helper: load per-dataset summaries from results directory -------------
load_rerun_summaries <- function(results_dir, method_name,
                                 summary_subpath = "summary.csv") {
  rows <- list()
  for (ds in datasets) {
    ds_num <- as.integer(sub("dataset", "", ds))
    summary_file <- file.path(results_dir, ds, summary_subpath)
    if (!file.exists(summary_file)) next
    s <- read.csv(summary_file, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds_num,
      PCC  = s$Mean[s$Metric == "PCC"],
      SSIM = s$Mean[s$Metric == "SSIM"],
      RMSE = s$Mean[s$Metric == "RMSE"],
      JSD  = s$Mean[s$Metric == "JSD"]
    )
  }
  if (length(rows) == 0) return(NULL)
  dt <- rbindlist(rows)
  data.frame(
    dataset = dt$dataset,
    Method  = method_name,
    stringsAsFactors = FALSE
  ) %>% bind_cols(dt[, c("PCC","SSIM","RMSE","JSD")]) %>%
    pivot_longer(c("PCC","SSIM","RMSE","JSD"),
                 names_to = "Metric", values_to = "Value")
}

# --- PhiSpace scran+log1p filtered ----------------------------------------
phi_long <- load_rerun_summaries(
  file.path(base_dir, "results_scran_log1p_filtered", "all_metrics"),
  "PhiSpace"
)

# --- RCTD -----------------------------------------------------------------
rctd_long <- load_rerun_summaries(
  file.path(base_dir, "results_rctd"),
  "RCTD",
  summary_subpath = file.path("metrics", "summary.csv")
)

# --- Cell2location --------------------------------------------------------
c2l_long <- load_rerun_summaries(
  file.path(base_dir, "results_cell2location", "all_metrics"),
  "Cell2location"
)

# --- TACCO (dataset-level summary) ----------------------------------------
tacco <- read.csv(file.path(base_dir,
    "results_tacco/TACCO_all_metrics.csv"), stringsAsFactors = FALSE)

tacco_long <- data.frame(
    dataset = tacco$dataset,
    Method  = "TACCO",
    stringsAsFactors = FALSE
) %>% bind_cols(tacco[, c("PCC","SSIM","RMSE","JSD")]) %>%
  pivot_longer(c("PCC","SSIM","RMSE","JSD"),
               names_to = "Metric", values_to = "Value")

# --- Seurat (re-run results) ---------------------------------------------
seurat_long <- load_rerun_summaries(
  file.path(base_dir, "results_seurat", "all_metrics"),
  "Seurat"
)

# --- SPOTlight (re-run results) -------------------------------------------
spotlight_long <- load_rerun_summaries(
  file.path(base_dir, "results_spotlight", "all_metrics"),
  "SPOTlight"
)

# --- Combine all -----------------------------------------------------------
all_long <- list(paper, phi_long, rctd_long, c2l_long, tacco_long,
                 seurat_long, spotlight_long)
all_long <- all_long[!sapply(all_long, is.null)]

df <- bind_rows(all_long)
df$dataset <- as.integer(df$dataset)
df <- df[!is.na(df$Value), ]

# Filter to our 7 methods only
target_methods <- c("PhiSpace", "RCTD", "Cell2location", "TACCO",
                    "Seurat", "SPOTlight", "DSTG")
df <- df %>% filter(Method %in% target_methods)

# --- Composite AS = (PCC + SSIM + (1-RMSE) + (1-JSD)) / 4 ----------------
df_wide <- df %>%
  group_by(dataset, Method, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metric, values_from = Value)

df_as <- df_wide %>%
  mutate(AS = (PCC + SSIM + (1 - RMSE) + (1 - JSD)) / 4) %>%
  select(dataset, Method, AS) %>%
  filter(!is.na(AS)) %>%
  mutate(Metric = "AS") %>%
  rename(Value = AS)

df <- bind_rows(df, df_as)
df$Metric <- factor(df$Metric, levels = c("PCC","SSIM","RMSE","JSD","AS"))

# Method ordering: PhiSpace first, then alphabetical
all_methods <- sort(unique(df$Method))
all_methods <- c("PhiSpace", setdiff(all_methods, "PhiSpace"))
df$Method <- factor(df$Method, levels = all_methods)

# Colours: PhiSpace highlighted red, others distinct
method_cols <- c(
  "PhiSpace"      = "#e74c3c",
  "Cell2location" = "#9b59b6",
  "DSTG"          = "#95a5a6",
  "RCTD"          = "#2ecc71",
  "Seurat"        = "#f39c12",
  "SPOTlight"     = "#1abc9c",
  "TACCO"         = "#3498db"
)

# --- Plot ------------------------------------------------------------------
p <- ggplot(df, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(width = 0.65, outlier.shape = 21, outlier.size = 1) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = method_cols) +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 13),
    axis.title.x       = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Score")

ggsave(file.path(fig_dir, "boxplot_all_methods.pdf"), p,
       width = 14, height = 4.5)
ggsave(file.path(fig_dir, "boxplot_all_methods.png"), p,
       width = 14, height = 4.5, dpi = 300)
cat("Performance boxplot saved.\n")

# --- Print ranking table ---------------------------------------------------
cat("\n=== Method ranking (median across datasets) ===\n")
rank_df <- df_wide %>%
  group_by(Method) %>%
  summarise(PCC  = median(PCC,  na.rm = TRUE),
            SSIM = median(SSIM, na.rm = TRUE),
            RMSE = median(RMSE, na.rm = TRUE),
            JSD  = median(JSD,  na.rm = TRUE),
            .groups = "drop") %>%
  mutate(AS = (PCC + SSIM + (1 - RMSE) + (1 - JSD)) / 4) %>%
  arrange(desc(AS))
print(as.data.frame(rank_df))

cat("\nDone. Figures saved to", fig_dir, "\n")
