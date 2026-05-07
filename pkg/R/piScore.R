#' Compute pi-score for ranking differentially expressed genes
#'
#' Combines log fold change and p-value into a single significance score
#' (Xiao et al. 2014). By default returns the signed variant
#' \code{log2FC * -log10(p)} for directional ranking; set \code{signed = FALSE}
#' for the paper's original \code{|log2FC| * -log10(p)}.
#'
#' @param de_results A data.frame-like object containing differential expression
#'   results.
#' @param logfc_col Character. Column containing log2 fold change. If the
#'   default is absent, common alternatives are detected automatically.
#' @param pval_col Character. Column containing p-values. If the default is
#'   absent, common alternatives are detected automatically.
#' @param signed Logical. If `TRUE`, preserve the direction of the log fold
#'   change. If `FALSE`, use the unsigned Xiao et al. definition.
#' @param pval_floor Numeric. Minimum p-value used before applying
#'   `-log10(p)`, avoiding infinite scores when p-values are zero.
#' @param add_column Logical. If `TRUE`, append `pi_score` to `de_results`.
#'   If `FALSE`, return a numeric vector.
#' @param rank Logical. If `TRUE` and `add_column = TRUE`, append `pi_rank`
#'   where rank 1 is the largest pi-score.
#'
#' @return If `add_column = TRUE`, a data.frame with `pi_score` and optionally
#'   `pi_rank` appended. Otherwise, a numeric vector aligned to the input rows.
#'
#' @references
#' Xiao Y, Hsiao T-H, Suresh U, et al. (2014) A novel significance score for
#'   gene selection and ranking. \emph{Bioinformatics} 30(6):801-807.
#'   \doi{10.1093/bioinformatics/btr671}
#'
#' @export
piScore <- function(de_results,
                    logfc_col = "avg_log2FC",
                    pval_col = "p_val_adj",
                    signed = TRUE,
                    pval_floor = .Machine$double.xmin,
                    add_column = TRUE,
                    rank = FALSE) {

  logfc_supplied <- !missing(logfc_col)
  pval_supplied <- !missing(pval_col)

  if (!is.data.frame(de_results)) {
    de_results <- as.data.frame(de_results)
  }

  if (!is.logical(signed) || length(signed) != 1 || is.na(signed)) {
    stop("signed must be TRUE or FALSE.")
  }
  if (!is.numeric(pval_floor) || length(pval_floor) != 1 ||
      is.na(pval_floor) || pval_floor <= 0) {
    stop("pval_floor must be a single positive number.")
  }
  if (!is.logical(add_column) || length(add_column) != 1 || is.na(add_column)) {
    stop("add_column must be TRUE or FALSE.")
  }
  if (!is.logical(rank) || length(rank) != 1 || is.na(rank)) {
    stop("rank must be TRUE or FALSE.")
  }

  detected <- .detect_de_columns(
    de_results = de_results,
    logfc_col = logfc_col,
    pval_col = pval_col,
    logfc_supplied = logfc_supplied,
    pval_supplied = pval_supplied
  )
  logfc_col <- detected$logfc_col
  pval_col <- detected$pval_col

  logfc <- de_results[[logfc_col]]
  pval <- de_results[[pval_col]]

  if (!is.numeric(logfc)) {
    stop("logfc_col must identify a numeric column.")
  }
  if (!is.numeric(pval)) {
    stop("pval_col must identify a numeric column.")
  }

  finite_p <- is.finite(pval)
  if (any(finite_p & (pval < 0 | pval > 1), na.rm = TRUE)) {
    stop("p-values must be between 0 and 1.")
  }

  pval_clamped <- pval
  clamp_idx <- !is.na(pval_clamped) & pval_clamped < pval_floor
  if (any(clamp_idx)) {
    warning("Some p-values were below pval_floor and were clamped before scoring.")
    pval_clamped[clamp_idx] <- pval_floor
  }

  logfc_term <- if (signed) logfc else abs(logfc)

  # Xiao et al. (2014): pi = |log2(FC)| * -log10(p); signed mode preserves logFC direction.
  score <- logfc_term * (-log10(pval_clamped))
  names(score) <- rownames(de_results)

  if (!add_column) {
    return(score)
  }

  de_results$pi_score <- score
  if (rank) {
    de_results$pi_rank <- base::rank(-score, ties.method = "min", na.last = "keep")
  }

  return(de_results)
}

.detect_de_columns <- function(de_results,
                               logfc_col,
                               pval_col,
                               logfc_supplied,
                               pval_supplied) {

  columns <- colnames(de_results)

  pairs <- list(
    c("avg_log2FC", "p_val_adj"),
    c("avg_log2FC", "p_val"),
    c("summary.logFC", "FDR"),
    c("logFC", "adj.P.Val"),
    c("logFC", "FDR"),
    c("log2FoldChange", "padj"),
    c("log2FoldChange", "pvalue")
  )

  if (logfc_supplied && !(logfc_col %in% columns)) {
    stop("logfc_col '", logfc_col, "' is not present in de_results.")
  }
  if (pval_supplied && !(pval_col %in% columns)) {
    stop("pval_col '", pval_col, "' is not present in de_results.")
  }

  if (logfc_supplied && pval_supplied) {
    return(list(logfc_col = logfc_col, pval_col = pval_col))
  }

  if (!logfc_supplied && logfc_col %in% columns &&
      !pval_supplied && pval_col %in% columns) {
    return(list(logfc_col = logfc_col, pval_col = pval_col))
  }

  if (logfc_supplied && !pval_supplied) {
    for (candidate in pairs) {
      if (identical(candidate[1], logfc_col) && candidate[2] %in% columns) {
        return(list(logfc_col = logfc_col, pval_col = candidate[2]))
      }
    }
  }

  if (!logfc_supplied && pval_supplied) {
    for (candidate in pairs) {
      if (candidate[1] %in% columns && identical(candidate[2], pval_col)) {
        return(list(logfc_col = candidate[1], pval_col = pval_col))
      }
    }
  }

  if (!logfc_supplied && !pval_supplied) {
    for (candidate in pairs) {
      if (all(candidate %in% columns)) {
        return(list(logfc_col = candidate[1], pval_col = candidate[2]))
      }
    }
  }

  stop(
    "Could not detect log fold change and p-value columns. ",
    "Provide logfc_col and pval_col explicitly."
  )
}
