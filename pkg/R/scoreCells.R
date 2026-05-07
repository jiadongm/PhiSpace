#' Score query cells against reference classes
#'
#' `scoreCells()` provides a lightweight correlation- and signature-based
#' fallback for regimes where a full PhiSpace model is not meaningful, such as
#' single-class references or references with very few classes. It stores score
#' matrices in `reducedDim()` slots with one column per reference class.
#'
#' @param reference A `SingleCellExperiment` or `SpatialExperiment` object.
#' @param query A `SingleCellExperiment` or `SpatialExperiment` object to score.
#' @param class_col Character or `NULL`. Column in `colData(reference)` defining
#'   reference classes. If `NULL`, all reference cells are treated as one class.
#' @param assay_ref Character. Assay in `reference` used for scoring.
#' @param assay_query Character. Assay in `query` used for scoring.
#' @param signature_method Character. Signature selection method. `"findMarkers"`
#'   uses `scran::findMarkers()` for multi-class references and falls back to
#'   mean-expression ranking for single-class references. `"mean"` uses
#'   mean-expression ranking for every class.
#' @param marker_test Character. Test type passed to `scran::findMarkers()`.
#' @param marker_pval_type Character. `pval.type` passed to
#'   `scran::findMarkers()`.
#' @param marker_direction Character. Direction used when selecting markers:
#'   `"up"`, `"down"`, or `"any"`.
#' @param n_top_genes Integer. Candidate gene pool per class before removing
#'   housekeeping genes.
#' @param n_signature_genes Integer. Number of final signature genes per class.
#' @param remove_housekeeping Logical. Whether to remove housekeeping genes from
#'   class signatures.
#' @param housekeeping_patterns Character vector of regex patterns identifying
#'   housekeeping genes.
#' @param housekeeping_genes Character vector of additional housekeeping genes.
#'   If `NULL`, a small built-in list is used.
#' @param shared_genes_only Logical. If `TRUE`, use the intersection of
#'   reference and query genes. If `FALSE`, reference and query row names must be
#'   identical.
#' @param scoring Character. Which score matrices to compute: `"both"`,
#'   `"correlation"`, or `"signature"`.
#' @param cor_method Character. Correlation method: `"spearman"` or `"pearson"`.
#' @param zscore_signature Logical. Whether to z-score signature scores. For
#'   multi-class references, z-scoring is across classes within each cell. For a
#'   single-class reference, z-scoring is across cells.
#' @param score_prefix Character. Prefix for output `reducedDim()` names.
#' @param BPPARAM BiocParallel parameter object passed to
#'   `scran::findMarkers()`.
#' @param verbose Logical. Whether to emit progress messages.
#'
#' @return The modified `query` object with score matrices stored in
#'   `reducedDim()` slots and run details stored in `metadata(query)$scoreCells`.
#'
#' @export
scoreCells <- function(reference,
                       query,
                       class_col = NULL,
                       assay_ref = "logcounts",
                       assay_query = "logcounts",
                       signature_method = c("findMarkers", "mean"),
                       marker_test = c("t", "wilcox", "binom"),
                       marker_pval_type = "any",
                       marker_direction = c("up", "down", "any"),
                       n_top_genes = 300,
                       n_signature_genes = 100,
                       remove_housekeeping = TRUE,
                       housekeeping_patterns = c("^Rpl", "^Rps", "^mt-", "^MT-", "^Mrpl", "^Mrps"),
                       housekeeping_genes = NULL,
                       shared_genes_only = TRUE,
                       scoring = c("both", "correlation", "signature"),
                       cor_method = c("spearman", "pearson"),
                       zscore_signature = TRUE,
                       score_prefix = "scoreCells",
                       BPPARAM = BiocParallel::SerialParam(),
                       verbose = TRUE) {

  signature_method <- match.arg(signature_method)
  marker_test <- match.arg(marker_test)
  marker_direction <- match.arg(marker_direction)
  scoring <- match.arg(scoring)
  cor_method <- match.arg(cor_method)

  .validate_score_cells_inputs(
    reference = reference,
    query = query,
    class_col = class_col,
    assay_ref = assay_ref,
    assay_query = assay_query,
    marker_pval_type = marker_pval_type,
    n_top_genes = n_top_genes,
    n_signature_genes = n_signature_genes,
    remove_housekeeping = remove_housekeeping,
    housekeeping_patterns = housekeeping_patterns,
    housekeeping_genes = housekeeping_genes,
    shared_genes_only = shared_genes_only,
    zscore_signature = zscore_signature,
    score_prefix = score_prefix,
    verbose = verbose
  )

  classes <- .get_score_cells_classes(reference, class_col)

  if (shared_genes_only) {
    feature_names <- intersect(rownames(reference), rownames(query))
  } else {
    if (!identical(rownames(reference), rownames(query))) {
      stop("When shared_genes_only = FALSE, reference and query row names must be identical.")
    }
    feature_names <- rownames(reference)
  }
  if (length(feature_names) == 0) {
    stop("Reference and query do not share any genes.")
  }

  reference <- reference[feature_names, , drop = FALSE]
  query <- query[feature_names, , drop = FALSE]

  ref_expr <- assay(reference, assay_ref)
  query_expr <- assay(query, assay_query)

  build <- .buildClassSignatures(
    reference = reference,
    ref_expr = ref_expr,
    classes = classes,
    signature_method = signature_method,
    marker_test = marker_test,
    marker_pval_type = marker_pval_type,
    marker_direction = marker_direction,
    n_top_genes = n_top_genes,
    n_signature_genes = n_signature_genes,
    remove_housekeeping = remove_housekeeping,
    housekeeping_patterns = housekeeping_patterns,
    housekeeping_genes = housekeeping_genes,
    BPPARAM = BPPARAM,
    verbose = verbose
  )

  computed_scores <- list()

  if (scoring %in% c("both", "correlation")) {
    if (verbose) {
      message("Computing ", cor_method, " correlation scores.")
    }
    correlation_scores <- .scoreCorrelation(
      query_expr = query_expr,
      centroids = build$centroids,
      method = cor_method
    )
    reducedDim(query, paste0(score_prefix, "_correlation")) <- correlation_scores
    computed_scores$correlation <- paste0(score_prefix, "_correlation")
  }

  if (scoring %in% c("both", "signature")) {
    if (verbose) {
      message("Computing signature scores.")
    }
    signature_scores <- .scoreSignature(
      query_expr = query_expr,
      signatures = build$signatures,
      zscore = zscore_signature
    )
    reducedDim(query, paste0(score_prefix, "_signature")) <- signature_scores
    computed_scores$signature <- paste0(score_prefix, "_signature")
  }

  score_metadata <- list(
    signatures = build$signatures,
    marker_stats = build$marker_stats,
    centroids = build$centroids,
    classes = names(build$signatures),
    params = list(
      class_col = class_col,
      assay_ref = assay_ref,
      assay_query = assay_query,
      signature_method = signature_method,
      marker_test = marker_test,
      marker_pval_type = marker_pval_type,
      marker_direction = marker_direction,
      n_top_genes = n_top_genes,
      n_signature_genes = n_signature_genes,
      remove_housekeeping = remove_housekeeping,
      housekeeping_patterns = housekeeping_patterns,
      shared_genes_only = shared_genes_only,
      scoring = scoring,
      cor_method = cor_method,
      zscore_signature = zscore_signature,
      score_prefix = score_prefix,
      reducedDims = computed_scores,
      timestamp = Sys.time()
    )
  )

  S4Vectors::metadata(query)$scoreCells <- score_metadata

  return(query)
}

.validate_score_cells_inputs <- function(reference,
                                         query,
                                         class_col,
                                         assay_ref,
                                         assay_query,
                                         marker_pval_type,
                                         n_top_genes,
                                         n_signature_genes,
                                         remove_housekeeping,
                                         housekeeping_patterns,
                                         housekeeping_genes,
                                         shared_genes_only,
                                         zscore_signature,
                                         score_prefix,
                                         verbose) {

  valid_object <- function(x) {
    inherits(x, "SingleCellExperiment") || inherits(x, "SpatialExperiment")
  }

  if (!valid_object(reference)) {
    stop("reference must be a SingleCellExperiment or SpatialExperiment object.")
  }
  if (!valid_object(query)) {
    stop("query must be a SingleCellExperiment or SpatialExperiment object.")
  }
  if (!(assay_ref %in% assayNames(reference))) {
    stop("assay_ref is not present in reference.")
  }
  if (!(assay_query %in% assayNames(query))) {
    stop("assay_query is not present in query.")
  }
  if (!is.null(class_col) && !(class_col %in% colnames(colData(reference)))) {
    stop("class_col is not a column of colData(reference).")
  }
  if (!is.character(marker_pval_type) || length(marker_pval_type) != 1 ||
      is.na(marker_pval_type)) {
    stop("marker_pval_type must be a single character value.")
  }
  if (!.is_positive_integer(n_top_genes)) {
    stop("n_top_genes must be a positive integer.")
  }
  if (!.is_positive_integer(n_signature_genes)) {
    stop("n_signature_genes must be a positive integer.")
  }
  if (!is.logical(remove_housekeeping) || length(remove_housekeeping) != 1 ||
      is.na(remove_housekeeping)) {
    stop("remove_housekeeping must be TRUE or FALSE.")
  }
  if (!is.character(housekeeping_patterns)) {
    stop("housekeeping_patterns must be a character vector.")
  }
  if (!is.null(housekeeping_genes) && !is.character(housekeeping_genes)) {
    stop("housekeeping_genes must be NULL or a character vector.")
  }
  if (!is.logical(shared_genes_only) || length(shared_genes_only) != 1 ||
      is.na(shared_genes_only)) {
    stop("shared_genes_only must be TRUE or FALSE.")
  }
  if (!is.logical(zscore_signature) || length(zscore_signature) != 1 ||
      is.na(zscore_signature)) {
    stop("zscore_signature must be TRUE or FALSE.")
  }
  if (!is.character(score_prefix) || length(score_prefix) != 1 ||
      is.na(score_prefix) || !nzchar(score_prefix)) {
    stop("score_prefix must be a non-empty character string.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
}

.get_score_cells_classes <- function(reference, class_col) {
  if (is.null(class_col)) {
    classes <- factor(rep("reference", ncol(reference)), levels = "reference")
  } else {
    values <- colData(reference)[[class_col]]
    if (any(is.na(values))) {
      stop("class_col cannot contain NA values.")
    }
    class_levels <- unique(as.character(values))
    classes <- factor(as.character(values), levels = class_levels)
  }

  if (nlevels(classes) == 0) {
    stop("No reference classes found.")
  }

  return(classes)
}

.buildClassSignatures <- function(reference,
                                  ref_expr,
                                  classes,
                                  signature_method,
                                  marker_test,
                                  marker_pval_type,
                                  marker_direction,
                                  n_top_genes,
                                  n_signature_genes,
                                  remove_housekeeping,
                                  housekeeping_patterns,
                                  housekeeping_genes,
                                  BPPARAM,
                                  verbose) {

  class_names <- levels(classes)
  centroids <- .compute_class_centroids(ref_expr, classes)

  if (is.null(housekeeping_genes)) {
    housekeeping_genes <- .default_housekeeping_genes()
  }
  housekeeping_flag <- .is_housekeeping_gene(
    genes = rownames(reference),
    patterns = housekeeping_patterns,
    extra_genes = housekeeping_genes
  )

  if (signature_method == "findMarkers" && length(class_names) > 1) {
    if (verbose) {
      message("Selecting class signatures with scran::findMarkers().")
    }
    markers <- scran::findMarkers(
      x = ref_expr,
      groups = classes,
      test.type = marker_test,
      pval.type = marker_pval_type,
      direction = marker_direction,
      BPPARAM = BPPARAM
    )
    signature_res <- .signatures_from_find_markers(
      markers = markers,
      class_names = class_names,
      n_top_genes = n_top_genes,
      n_signature_genes = n_signature_genes,
      remove_housekeeping = remove_housekeeping,
      housekeeping_flag = housekeeping_flag,
      marker_direction = marker_direction
    )
  } else {
    if (verbose) {
      message("Selecting class signatures by mean expression ranking.")
    }
    signature_res <- .signatures_from_mean(
      centroids = centroids,
      n_top_genes = n_top_genes,
      n_signature_genes = n_signature_genes,
      remove_housekeeping = remove_housekeeping,
      housekeeping_flag = housekeeping_flag
    )
  }

  return(
    list(
      signatures = signature_res$signatures,
      marker_stats = signature_res$marker_stats,
      centroids = centroids
    )
  )
}

.compute_class_centroids <- function(ref_expr, classes) {
  class_names <- levels(classes)
  centroids <- vapply(
    class_names,
    function(cl) {
      rowMeans(ref_expr[, classes == cl, drop = FALSE], na.rm = TRUE)
    },
    numeric(nrow(ref_expr))
  )
  rownames(centroids) <- rownames(ref_expr)
  colnames(centroids) <- class_names
  return(centroids)
}

.signatures_from_find_markers <- function(markers,
                                          class_names,
                                          n_top_genes,
                                          n_signature_genes,
                                          remove_housekeeping,
                                          housekeeping_flag,
                                          marker_direction) {

  signatures <- vector("list", length(class_names))
  names(signatures) <- class_names
  stats_list <- vector("list", length(class_names))
  names(stats_list) <- class_names

  for (cl in class_names) {
    tab <- as.data.frame(markers[[cl]])
    tab$gene <- rownames(markers[[cl]])

    if ("summary.logFC" %in% colnames(tab)) {
      tab <- piScore(
        tab,
        logfc_col = "summary.logFC",
        pval_col = if ("FDR" %in% colnames(tab)) "FDR" else "p.value",
        rank = FALSE
      )
    } else if ("FDR" %in% colnames(tab)) {
      tab <- tab[order(tab$FDR, na.last = TRUE), , drop = FALSE]
      tab$pi_score <- NA_real_
    } else {
      tab$pi_score <- NA_real_
    }

    if ("summary.logFC" %in% colnames(tab)) {
      if (marker_direction == "up") {
        tab <- tab[is.na(tab$summary.logFC) | tab$summary.logFC > 0, , drop = FALSE]
        tab <- tab[order(-tab$pi_score, na.last = TRUE), , drop = FALSE]
      } else if (marker_direction == "down") {
        tab <- tab[is.na(tab$summary.logFC) | tab$summary.logFC < 0, , drop = FALSE]
        tab <- tab[order(tab$pi_score, na.last = TRUE), , drop = FALSE]
      } else {
        tab <- tab[order(-abs(tab$pi_score), na.last = TRUE), , drop = FALSE]
      }
    }

    tab <- tab[seq_len(min(nrow(tab), n_top_genes)), , drop = FALSE]
    tab$class <- cl
    tab$rank <- seq_len(nrow(tab))
    tab$pi_rank <- tab$rank
    tab$is_housekeeping <- housekeeping_flag[tab$gene]
    tab$is_housekeeping[is.na(tab$is_housekeeping)] <- FALSE

    eligible <- if (remove_housekeeping) !tab$is_housekeeping else rep(TRUE, nrow(tab))
    kept_idx <- which(eligible)[seq_len(min(sum(eligible), n_signature_genes))]
    tab$kept <- seq_len(nrow(tab)) %in% kept_idx
    signatures[[cl]] <- tab$gene[tab$kept]
    stats_list[[cl]] <- tab
  }

  marker_stats <- .rbind_fill_data_frames(stats_list)
  rownames(marker_stats) <- NULL

  return(list(signatures = signatures, marker_stats = marker_stats))
}

.signatures_from_mean <- function(centroids,
                                  n_top_genes,
                                  n_signature_genes,
                                  remove_housekeeping,
                                  housekeeping_flag) {

  class_names <- colnames(centroids)
  signatures <- vector("list", length(class_names))
  names(signatures) <- class_names
  stats_list <- vector("list", length(class_names))
  names(stats_list) <- class_names

  for (cl in class_names) {
    mean_expr <- centroids[, cl]
    ord <- order(mean_expr, decreasing = TRUE, na.last = TRUE)
    genes <- rownames(centroids)[ord]
    tab <- data.frame(
      gene = genes,
      class = cl,
      mean_expr = mean_expr[ord],
      rank = seq_along(genes),
      is_housekeeping = housekeeping_flag[genes],
      stringsAsFactors = FALSE
    )
    tab$is_housekeeping[is.na(tab$is_housekeeping)] <- FALSE
    tab <- tab[seq_len(min(nrow(tab), n_top_genes)), , drop = FALSE]

    eligible <- if (remove_housekeeping) !tab$is_housekeeping else rep(TRUE, nrow(tab))
    kept_idx <- which(eligible)[seq_len(min(sum(eligible), n_signature_genes))]
    tab$kept <- seq_len(nrow(tab)) %in% kept_idx
    signatures[[cl]] <- tab$gene[tab$kept]
    stats_list[[cl]] <- tab
  }

  marker_stats <- .rbind_fill_data_frames(stats_list)
  rownames(marker_stats) <- NULL

  return(list(signatures = signatures, marker_stats = marker_stats))
}

.scoreCorrelation <- function(query_expr, centroids, method) {
  common <- intersect(rownames(query_expr), rownames(centroids))
  query_expr <- query_expr[common, , drop = FALSE]
  centroids <- centroids[common, , drop = FALSE]

  scores <- stats::cor(
    x = as.matrix(query_expr),
    y = as.matrix(centroids),
    use = "pairwise.complete.obs",
    method = method
  )
  scores <- as.matrix(scores)
  rownames(scores) <- colnames(query_expr)
  colnames(scores) <- colnames(centroids)
  return(scores)
}

.scoreSignature <- function(query_expr, signatures, zscore) {
  class_names <- names(signatures)
  scores <- vapply(
    class_names,
    function(cl) {
      genes <- intersect(signatures[[cl]], rownames(query_expr))
      if (length(genes) == 0) {
        warning("No signature genes for class '", cl, "' are present in query.")
        return(rep(NA_real_, ncol(query_expr)))
      }
      colMeans(query_expr[genes, , drop = FALSE], na.rm = TRUE)
    },
    numeric(ncol(query_expr))
  )
  scores <- as.matrix(scores)
  rownames(scores) <- colnames(query_expr)
  colnames(scores) <- class_names

  if (zscore) {
    if (ncol(scores) == 1) {
      scores[, 1] <- .zscore_vector(scores[, 1])
    } else {
      scores <- t(apply(scores, 1, .zscore_vector))
      colnames(scores) <- class_names
      rownames(scores) <- colnames(query_expr)
    }
  }

  return(scores)
}

.zscore_vector <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    out <- x - x_mean
  } else {
    out <- (x - x_mean) / x_sd
  }
  return(out)
}

.is_housekeeping_gene <- function(genes, patterns, extra_genes = NULL) {
  out <- rep(FALSE, length(genes))
  names(out) <- genes

  if (length(patterns) > 0) {
    pattern_hits <- vapply(
      patterns,
      function(pattern) grepl(pattern, genes),
      logical(length(genes))
    )
    if (is.null(dim(pattern_hits))) {
      out <- out | pattern_hits
    } else {
      out <- out | rowSums(pattern_hits) > 0
    }
  }

  if (!is.null(extra_genes)) {
    out <- out | names(out) %in% extra_genes
  }

  return(out)
}

.rbind_fill_data_frames <- function(x) {
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(
    x,
    function(tab) {
      missing_names <- setdiff(all_names, names(tab))
      for (nm in missing_names) {
        tab[[nm]] <- NA
      }
      tab[, all_names, drop = FALSE]
    }
  )
  do.call(rbind, x)
}

.default_housekeeping_genes <- function() {
  c(
    "ACTB", "Actb", "GAPDH", "Gapdh", "B2M", "B2m", "MALAT1", "Malat1",
    "TMSB4X", "Tmsb4x", "TMSB10", "Tmsb10", "HSP90AA1", "Hsp90aa1",
    "HSP90AB1", "Hsp90ab1", "HSPA1A", "Hspa1a", "HSPA1B", "Hspa1b",
    "HSPA8", "Hspa8", "HSPB1", "Hspb1", "FTL", "Ftl", "FTH1", "Fth1",
    "JUN", "Jun", "FOS", "Fos", "EGR1", "Egr1", "HBB", "Hbb",
    "HBA1", "Hba1", "HBA2", "Hba2"
  )
}

.is_positive_integer <- function(x) {
  length(x) == 1 && is.numeric(x) && !is.na(x) &&
    x == as.integer(x) && x > 0
}
