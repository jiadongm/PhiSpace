# Result checks for vignette builds. Each metric function summarises the files a
# vignette writes to output/. check_results() compares the summary with the
# reviewed values in scripts/<slug>-results.tsv: a metric with a tolerance must
# be within that absolute tolerance; a metric without one must match exactly.

result_metrics <- list()

result_metrics$StereoSeq <- function(output_dir) {
  read <- function(name) qs2::qs_read(file.path(output_dir, name), validate_checksum = TRUE)
  m <- list()

  # Clone density estimates; the vignette uses their range for the colour scale.
  clones <- read("cloneKDEres.qs2")
  kde <- unlist(lapply(clones, function(x) x$kde_res$estimate))
  m$kde_n_clones <- length(clones)
  m$kde_estimate_min <- min(kde)
  m$kde_estimate_max <- max(kde)
  m$kde_estimate_sum <- sum(kde)

  # Barcode k-means for k = 2, ..., 10.
  bc <- read("barcode-seq_clustering.qs2")
  for (i in seq_along(bc)) m[[sprintf("barcode_k%d_tot_withinss", length(bc[[i]]$size))]] <- bc[[i]]$tot.withinss

  # Bridge annotation of the Stereo-seq bins.
  phi <- read("PhiRes.qs2")
  score <- phi$PhiSpaceScore
  m$phi_n_bins <- nrow(score)
  m$phi_n_celltypes <- ncol(score)
  m$phi_celltypes_sha256 <- digest::digest(colnames(score), algo = "sha256")
  m$phi_bins_sha256 <- digest::digest(rownames(score), algo = "sha256")
  m$phi_n_selected_features <- length(phi$selectedFeat)
  m$phi_score_mean <- mean(score)
  m$phi_score_sd <- sd(as.vector(score))
  for (ct in c("Neutro(Spleen)", "Granulo(BM)", "T2(Neutro)", "ErythBla(BM)", "Macro(Spleen)", "Naive B(BM)")) {
    m[[paste0("phi_sd_", ct)]] <- sd(score[, ct])
  }

  # Niche clustering of PhiSpace scores and of gene expression. Cluster labels
  # are arbitrary, so sizes are sorted.
  cluster_files <- c(phi_cluster = "PhiClustRes.qs2", gex_cluster = "GexClustRes.qs2")
  for (name in names(cluster_files)) {
    km <- read(cluster_files[[name]])
    sizes <- sort(km$size)
    for (i in seq_along(sizes)) m[[sprintf("%s_size_rank%d", name, i)]] <- sizes[[i]]
    m[[paste0(name, "_tot_withinss")]] <- km$tot.withinss
  }

  # Niche enrichment scores.
  sig <- read("sigScores.qs2")
  m$sig_n_celltypes <- nrow(sig)
  m$sig_niches <- paste(colnames(sig), collapse = ",")
  m$sig_min <- min(sig)
  m$sig_max <- max(sig)
  m$sig_abs_sum <- sum(abs(sig))
  m
}

result_metrics$Visium <- function(output_dir) {
  res <- qs2::qs_read(file.path(output_dir, "combo_PhiRes.qs2"), validate_checksum = TRUE)
  m <- list()
  m$samples <- paste(names(res), collapse = ",")
  m$n_celltypes <- unique(vapply(res, ncol, integer(1)))
  m$celltypes_sha256 <- digest::digest(colnames(res[[1]]), algo = "sha256")
  if (!all(vapply(res, function(x) identical(colnames(x), colnames(res[[1]])), logical(1)))) {
    stop("Visium samples have different cell type columns.")
  }
  for (s in names(res)) {
    m[[paste0(s, "_n_spots")]] <- nrow(res[[s]])
    m[[paste0(s, "_score_mean")]] <- mean(res[[s]])
    m[[paste0(s, "_score_sd")]] <- sd(as.vector(res[[s]]))
  }
  m$P11_T3_sd_B_cells <- sd(res[["P11_T3"]][, "B cells"])
  m
}

result_metrics$BridgeAnnotation <- function(output_dir) {
  m <- list()
  top <- list()
  files <- c(rna = "PhiResRNA.rds", peaks = "PhiResATACpeaks.rds", ga = "PhiResATAC_GA.rds")
  for (name in names(files)) {
    res <- readRDS(file.path(output_dir, files[[name]]))
    score <- res$PhiSpaceScore
    m[[paste0(name, "_n_cells")]] <- nrow(score)
    m[[paste0(name, "_n_celltypes")]] <- ncol(score)
    m[[paste0(name, "_celltypes_sha256")]] <- digest::digest(colnames(score), algo = "sha256")
    m[[paste0(name, "_cells_sha256")]] <- digest::digest(rownames(score), algo = "sha256")
    m[[paste0(name, "_ncomp")]] <- res$ncomp
    m[[paste0(name, "_n_selected_features")]] <- length(res$selectedFeat)
    m[[paste0(name, "_score_mean")]] <- mean(score)
    m[[paste0(name, "_score_sd")]] <- sd(as.vector(score))
    norm <- normPhiScores(score)
    m[[paste0(name, "_norm_rowmax_mean")]] <- mean(apply(norm, 1, max))
    top[[name]] <- colnames(norm)[max.col(norm, ties.method = "first")]
  }
  # Agreement of the top labels from peaks and gene activity, the two query annotations.
  m$peaks_ga_top_label_agreement <- mean(top$peaks == top$ga)
  m
}

result_metrics$getting_started <- function(output_dir) {
  read <- function(name) qs2::qs_read(file.path(output_dir, name), validate_checksum = TRUE)
  m <- list()

  # Subsampled and annotated query.
  query <- read("quPhi.qs2")
  m$query_n_cells <- ncol(query)
  m$query_n_genes <- nrow(query)
  types <- table(query$mainTypes)
  for (ct in names(types)) m[[paste0("query_n_", ct)]] <- types[[ct]]
  score <- SingleCellExperiment::reducedDim(query, "PhiSpace")
  m$query_phenotypes <- paste(colnames(score), collapse = ",")
  m$query_cells_sha256 <- digest::digest(rownames(score), algo = "sha256")
  m$query_score_mean <- mean(score)
  m$query_score_sd <- sd(as.vector(score))
  for (ph in colnames(score)) m[[paste0("query_mean_", ph)]] <- mean(score[, ph])

  # Reference predictions used for the phenotype-space PCA.
  ref_score <- SingleCellExperiment::reducedDim(read("refPhi.qs2"), "PhiSpace")
  m$reference_n_samples <- nrow(ref_score)
  m$reference_score_sd <- sd(as.vector(ref_score))

  # Feature importance; the vignette reports the number of genes for nfeat = 300.
  imp <- read("tuneRes.qs2")$impScores
  m$imp_n_genes <- nrow(imp)
  m$imp_n_phenotypes <- ncol(imp)
  m$imp_abs_sum <- sum(abs(imp))
  ordered <- selectFeat(imp)$orderedFeatMat
  m$n_genes_top300 <- length(unique(as.character(ordered[1:300, ])))
  m
}

result_metrics[["CITE-seq"]] <- function(output_dir) {
  m <- list()
  top <- list()
  files <- c(adt = "doublePhenoADTPhiRes.rds", rna = "doublePhenoRNAPhiRes.rds")
  for (name in names(files)) {
    res <- readRDS(file.path(output_dir, files[[name]]))
    for (part in c("PhiSpaceScore", "YrefHat")) {
      score <- res[[part]]
      key <- paste0(name, if (part == "YrefHat") "_ref" else "_query")
      m[[paste0(key, "_n_cells")]] <- nrow(score)
      m[[paste0(key, "_n_phenotypes")]] <- ncol(score)
      m[[paste0(key, "_phenotypes_sha256")]] <- digest::digest(colnames(score), algo = "sha256")
      m[[paste0(key, "_cells_sha256")]] <- digest::digest(rownames(score), algo = "sha256")
      m[[paste0(key, "_score_mean")]] <- mean(score)
      m[[paste0(key, "_score_sd")]] <- sd(as.vector(score))
    }
    m[[paste0(name, "_ncomp")]] <- res$ncomp
    m[[paste0(name, "_n_selected_features")]] <- length(res$selectedFeat)
  }
  # UMAP coordinates depend on floating-point details, so only their cells are checked.
  umaps <- c(umap_adt = "doublePhenoADTPhiUMAPquery.rds", umap_rna = "doublePhenoRNAPhiUMAPquery.rds",
             umap_combo = "doublePhenoComboPhiUMAPquery.rds")
  for (name in names(umaps)) {
    layout <- readRDS(file.path(output_dir, umaps[[name]]))$layout
    m[[paste0(name, "_n_cells")]] <- nrow(layout)
    m[[paste0(name, "_cells_sha256")]] <- digest::digest(rownames(layout), algo = "sha256")
  }
  m
}

result_metrics$PerturbSeq <- function(output_dir) {
  # Normalised PhiSpace scores from PhiSpace(): reducedDim(query, "PhiSpace").
  norm <- qs2::qs_read(file.path(output_dir, "PhiSpaceScores.qs2"), validate_checksum = TRUE)
  m <- list()
  m$query_n_cells <- nrow(norm)
  m$query_n_celltypes <- ncol(norm)
  m$query_celltypes_sha256 <- digest::digest(colnames(norm), algo = "sha256")
  m$query_cells_sha256 <- digest::digest(rownames(norm), algo = "sha256")
  m$norm_score_mean <- mean(norm)
  m$norm_score_sd <- sd(as.vector(norm))
  # Scores behind the activation and Th1/Th2 analyses.
  for (ct in c("T cells, CD4+, naive, stimulated", "T cells, CD8+, naive, stimulated",
               "T cells, CD4+, Th1", "T cells, CD4+, Th2")) {
    m[[paste0("norm_mean_", ct)]] <- mean(norm[, ct])
    m[[paste0("norm_sd_", ct)]] <- sd(norm[, ct])
  }
  m
}

result_metrics$CosMx <- function(output_dir) {
  read <- function(name) qs2::qs_read(file.path(output_dir, name), validate_checksum = TRUE)
  m <- list()
  # Lung5_Rep1 annotation by the four lineage references (normalised scores).
  sc_list <- read("CosMxLung5Rep1PhiRes4Refs.qs2")
  m$lineages <- paste(names(sc_list), collapse = ",")
  for (lineage in names(sc_list)) {
    sc <- sc_list[[lineage]]
    key <- paste0("lung5rep1_", lineage)
    m[[paste0(key, "_n_cells")]] <- nrow(sc)
    m[[paste0(key, "_n_celltypes")]] <- ncol(sc)
    m[[paste0(key, "_celltypes_sha256")]] <- digest::digest(colnames(sc), algo = "sha256")
    m[[paste0(key, "_cells_sha256")]] <- digest::digest(rownames(sc), algo = "sha256")
    m[[paste0(key, "_score_mean")]] <- mean(sc)
    m[[paste0(key, "_score_sd")]] <- sd(as.vector(sc))
  }
  # PhiSpace niches of Lung5_Rep1. Labels are arbitrary, so sizes are sorted;
  # sizes allow 1% because k-means assignments of borderline cells depend on BLAS.
  km <- read("Lung5_Rep1_PhiClusts4Refs.qs2")
  sizes <- sort(km$size)
  for (i in seq_along(sizes)) m[[sprintf("niche_size_rank%d", i)]] <- sizes[[i]]
  m$niche_tot_withinss <- km$tot.withinss
  m
}

format_metrics <- function(metrics) {
  vapply(metrics, function(x) {
    if (length(x) != 1) stop("Each result metric must be a single value.")
    if (is.numeric(x)) format(x, digits = 15) else as.character(x)
  }, character(1))
}

check_results <- function(vignette, output_dir,
                          slug = gsub("_", "-", tolower(vignette))) {
  observed <- format_metrics(result_metrics[[vignette]](output_dir))
  out <- Sys.getenv("PHISPACE_RESULTS_OUT")
  if (nzchar(out)) {
    write.table(data.frame(metric = names(observed), value = observed), out,
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  expected_file <- file.path("scripts", paste0(slug, "-results.tsv"))
  expected <- read.delim(expected_file, colClasses = "character", na.strings = "")
  missing <- setdiff(expected$metric, names(observed))
  extra <- setdiff(names(observed), expected$metric)
  if (length(missing) || length(extra)) {
    stop("Result metrics do not match ", expected_file, ". Missing: ",
         paste(missing, collapse = ", "), ". Unexpected: ", paste(extra, collapse = ", "))
  }
  obs <- observed[expected$metric]
  tol <- as.numeric(expected$tolerance)
  ok <- obs == expected$value
  num <- !is.na(tol)
  ok[num] <- abs(as.numeric(obs[num]) - as.numeric(expected$value[num])) <= tol[num]
  ok[is.na(ok)] <- FALSE
  if (!all(ok)) {
    stop("Results differ from ", expected_file, ":\n",
         paste0("  ", expected$metric[!ok], ": expected ", expected$value[!ok],
                ", observed ", obs[!ok], collapse = "\n"),
         "\nReview the change before updating the expected values.")
  }
  message("Checked ", length(ok), " ", vignette, " result metrics against ", expected_file)
  invisible(observed)
}
