test_that("generated signatures drive both score types", {
  fixture <- make_scorecells_fixture()

  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = "celltype",
    signature_method = "mean",
    n_top_genes = 4,
    n_signature_genes = 2,
    verbose = FALSE
  )

  cor_scores <- SingleCellExperiment::reducedDim(out, "scoreCells_correlation")
  sig_scores <- SingleCellExperiment::reducedDim(out, "scoreCells_signature")
  meta <- S4Vectors::metadata(out)$scoreCells
  ref_expr <- SummarizedExperiment::assay(fixture$reference, "logcounts")
  query_expr <- SummarizedExperiment::assay(fixture$query, "logcounts")

  expect_equal(dim(cor_scores), c(3, 2))
  expect_equal(dim(sig_scores), c(3, 2))
  expect_equal(rownames(cor_scores), colnames(fixture$query))
  expect_equal(colnames(cor_scores), c("A", "B"))
  expect_identical(meta$signature_source, "generated")
  expect_identical(meta$correlation_genes, meta$signatures)
  expect_false(meta$params$signature_genes_supplied)

  for (cl in c("A", "B")) {
    genes <- meta$signatures[[cl]]
    centroid <- rowMeans(
      ref_expr[, fixture$reference$celltype == cl, drop = FALSE]
    )
    query_name <- if (cl == "A") "query1" else "query2"
    manual <- stats::cor(
      query_expr[genes, query_name],
      centroid[genes],
      method = "spearman"
    )
    expect_equal(unname(cor_scores[query_name, cl]), unname(manual))
  }

  raw_signature <- vapply(
    meta$signatures,
    function(genes) colMeans(query_expr[genes, , drop = FALSE]),
    numeric(ncol(query_expr))
  )
  manual_signature <- t(apply(raw_signature, 1, PhiSpace:::.zscore_vector))
  dimnames(manual_signature) <- dimnames(sig_scores)
  expect_equal(sig_scores, manual_signature)
  expect_true(all(c("GeneA1", "GeneA2") %in% meta$signatures$A))
  expect_false(any(c("Rpl1", "ACTB") %in% unlist(meta$signatures, use.names = FALSE)))
})

test_that("generated correlation ignores non-signature genes", {
  fixture <- make_scorecells_fixture()
  args <- list(
    reference = fixture$reference,
    class_col = "celltype",
    signature_method = "mean",
    n_top_genes = 4,
    n_signature_genes = 2,
    scoring = "correlation",
    verbose = FALSE
  )

  baseline <- do.call(scoreCells, c(args, list(query = fixture$query)))
  changed_query <- fixture$query
  SummarizedExperiment::assay(changed_query, "logcounts")["Rpl1", ] <-
    c(-1e6, 0, 1e6)
  changed <- do.call(scoreCells, c(args, list(query = changed_query)))

  expect_identical(
    SingleCellExperiment::reducedDimNames(baseline),
    "scoreCells_correlation"
  )
  expect_equal(
    SingleCellExperiment::reducedDim(baseline, "scoreCells_correlation"),
    SingleCellExperiment::reducedDim(changed, "scoreCells_correlation")
  )
})

test_that("single-class generated correlation uses generated signature genes", {
  fixture <- make_scorecells_fixture()
  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = NULL,
    signature_method = "mean",
    n_top_genes = 6,
    n_signature_genes = 3,
    scoring = "correlation",
    verbose = FALSE
  )

  meta <- S4Vectors::metadata(out)$scoreCells
  genes <- meta$signatures$reference
  ref_expr <- SummarizedExperiment::assay(fixture$reference, "logcounts")
  query_expr <- SummarizedExperiment::assay(fixture$query, "logcounts")
  centroid <- rowMeans(ref_expr)
  manual <- stats::cor(
    query_expr[genes, "query1"],
    centroid[genes],
    method = "spearman"
  )
  observed <- SingleCellExperiment::reducedDim(
    out, "scoreCells_correlation"
  )

  expect_equal(unname(observed["query1", "reference"]), unname(manual))
})

test_that("a character vector supplies a one-class signature", {
  fixture <- make_scorecells_fixture()
  supplied <- c("ACTB", "GeneA1", "GeneA2", "ACTB")

  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = NULL,
    signature_genes = supplied,
    n_top_genes = 1,
    n_signature_genes = 1,
    remove_housekeeping = TRUE,
    scoring = "both",
    zscore_signature = FALSE,
    verbose = FALSE
  )

  meta <- S4Vectors::metadata(out)$scoreCells
  expect_identical(meta$signatures, list(reference = unique(supplied)))
  expect_identical(meta$correlation_genes, meta$signatures)
  expect_identical(meta$signature_source, "user")
  expect_true(meta$params$signature_genes_supplied)
  expect_null(meta$marker_stats)

  query_expr <- SummarizedExperiment::assay(fixture$query, "logcounts")
  expected_signature <- colMeans(query_expr[unique(supplied), , drop = FALSE])
  observed_signature <- SingleCellExperiment::reducedDim(
    out, "scoreCells_signature"
  )
  expect_equal(
    unname(observed_signature[, "reference"]),
    unname(expected_signature)
  )
})

test_that("named multi-class signatures are used class by class", {
  fixture <- make_scorecells_fixture()
  supplied <- list(
    A = c("GeneA1", "GeneA2", "Rpl1"),
    B = c("GeneB1", "GeneB2", "ACTB")
  )

  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = "celltype",
    signature_genes = supplied,
    scoring = "both",
    zscore_signature = FALSE,
    score_prefix = "custom",
    verbose = FALSE
  )

  expect_setequal(
    SingleCellExperiment::reducedDimNames(out),
    c("custom_correlation", "custom_signature")
  )
  scores <- SingleCellExperiment::reducedDim(out, "custom_correlation")
  ref_expr <- SummarizedExperiment::assay(fixture$reference, "logcounts")
  query_expr <- SummarizedExperiment::assay(fixture$query, "logcounts")

  for (cl in names(supplied)) {
    centroid <- rowMeans(
      ref_expr[, fixture$reference$celltype == cl, drop = FALSE]
    )
    manual <- apply(
      query_expr[supplied[[cl]], , drop = FALSE],
      2,
      stats::cor,
      y = centroid[supplied[[cl]]],
      method = "spearman"
    )
    expect_equal(unname(scores[, cl]), unname(manual))
  }

  meta <- S4Vectors::metadata(out)$scoreCells
  expect_identical(meta$signatures, supplied)
  expect_identical(meta$correlation_genes, supplied)
  expect_true("ACTB" %in% meta$signatures$B)
})

test_that("supplied signatures work for signature-only scoring", {
  fixture <- make_scorecells_fixture()
  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = "celltype",
    signature_genes = list(
      B = c("GeneB1", "GeneB2"),
      A = c("GeneA1", "GeneA2")
    ),
    scoring = "signature",
    verbose = FALSE
  )

  expect_identical(
    SingleCellExperiment::reducedDimNames(out),
    "scoreCells_signature"
  )
  expect_identical(
    names(S4Vectors::metadata(out)$scoreCells$signatures),
    c("A", "B")
  )
  expect_null(S4Vectors::metadata(out)$scoreCells$correlation_genes)
})

test_that("signature validation rejects malformed inputs", {
  fixture <- make_scorecells_fixture()
  call_score <- function(signature_genes, scoring = "signature") {
    scoreCells(
      fixture$reference,
      fixture$query,
      class_col = "celltype",
      signature_genes = signature_genes,
      scoring = scoring,
      verbose = FALSE
    )
  }

  expect_error(
    call_score(c("GeneA1", "GeneA2")),
    "named list for a multi-class"
  )
  expect_error(
    call_score(list(c("GeneA1"), c("GeneB1"))),
    "named list"
  )
  expect_error(
    call_score(list(A = "GeneA1")),
    "missing: B"
  )
  expect_error(
    call_score(list(A = "GeneA1", B = "GeneB1", C = "GeneA2")),
    "unknown: C"
  )
  duplicated_names <- list("GeneA1", "GeneB1")
  names(duplicated_names) <- c("A", "A")
  expect_error(call_score(duplicated_names), "duplicated class names")
  expect_error(
    call_score(list(A = 1, B = "GeneB1")),
    "non-empty character vector"
  )
  expect_error(
    call_score(list(A = character(), B = "GeneB1")),
    "non-empty character vector"
  )
  expect_error(
    call_score(list(A = c("GeneA1", NA), B = "GeneB1")),
    "NA or blank"
  )
  expect_error(
    call_score(list(A = c("GeneA1", " "), B = "GeneB1")),
    "NA or blank"
  )
})

test_that("unavailable and insufficient correlation genes are handled clearly", {
  fixture <- make_scorecells_fixture()

  expect_warning(
    out <- scoreCells(
      fixture$reference,
      fixture$query,
      class_col = NULL,
      signature_genes = c("GeneA1", "missing", "GeneA1"),
      scoring = "signature",
      verbose = FALSE
    ),
    "1 gene.*unavailable"
  )
  expect_identical(
    S4Vectors::metadata(out)$scoreCells$signatures,
    list(reference = "GeneA1")
  )

  expect_error(
    scoreCells(
      fixture$reference,
      fixture$query,
      class_col = NULL,
      signature_genes = "GeneA1",
      scoring = "correlation",
      verbose = FALSE
    ),
    "at least two effective"
  )
  expect_error(
    scoreCells(
      fixture$reference,
      fixture$query,
      class_col = NULL,
      signature_genes = c("Rpl1", "ACTB"),
      scoring = "correlation",
      verbose = FALSE
    ),
    "nonzero centroid variance"
  )
})

test_that("constant query profiles yield NA correlations", {
  fixture <- make_scorecells_fixture()
  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = NULL,
    signature_genes = c("GeneA1", "GeneA2"),
    scoring = "correlation",
    verbose = FALSE
  )

  scores <- SingleCellExperiment::reducedDim(out, "scoreCells_correlation")
  expect_true(is.na(scores["query3", "reference"]))
  expect_equal(dimnames(scores), list(colnames(fixture$query), "reference"))
})

test_that("score dimensions are query cells by reference classes", {
  fixture <- make_scorecells_fixture()
  one_query <- fixture$query[, "query1", drop = FALSE]

  out <- scoreCells(
    fixture$reference,
    one_query,
    class_col = "celltype",
    signature_genes = list(
      A = c("GeneA1", "GeneA2", "Rpl1"),
      B = c("GeneB1", "GeneB2", "ACTB")
    ),
    scoring = "both",
    verbose = FALSE
  )

  expect_equal(
    dim(SingleCellExperiment::reducedDim(out, "scoreCells_correlation")),
    c(1, 2)
  )
  expect_equal(
    dim(SingleCellExperiment::reducedDim(out, "scoreCells_signature")),
    c(1, 2)
  )
  expect_equal(
    dimnames(SingleCellExperiment::reducedDim(out, "scoreCells_signature")),
    list("query1", c("A", "B"))
  )
})

test_that("vectorised correlation matches stats::cor column by column", {
  set.seed(1)
  X <- matrix(stats::rpois(40 * 30, 1), 40, 30)  # many ties
  X[, 5] <- 2                                     # a constant column
  y <- stats::rpois(40, 3) + stats::runif(40)
  expect_equal(.colRanksAverage(X), apply(X, 2, rank))
  for (method in c("spearman", "pearson")) {
    expected <- suppressWarnings(apply(X, 2, stats::cor, y = y, method = method))
    expected[5] <- NA_real_
    expect_equal(.corWithVector(X, y, method), expected, tolerance = 1e-12)
  }
  expect_length(.corWithVector(X[, 0, drop = FALSE], y, "pearson"), 0)
})

test_that("correlation scores of cells with missing values are unchanged", {
  set.seed(2)
  genes <- paste0("gene", 1:30)
  centroids <- matrix(stats::runif(60), 30, 2, dimnames = list(genes, c("A", "B")))
  centroids[3, "A"] <- NA
  query <- matrix(stats::rnorm(30 * 12), 30, 12, dimnames = list(genes, paste0("c", 1:12)))
  query[1:4, 2] <- NA
  query[, 7] <- 1
  sigs <- list(A = genes[1:20], B = genes[5:30])
  out <- .scoreCorrelation(query, centroids, sigs, "spearman")
  for (cl in c("A", "B")) {
    g <- sigs[[cl]]
    for (i in seq_len(ncol(query))) {
      ok <- is.finite(query[g, i]) & is.finite(centroids[g, cl])
      v <- query[g, i][ok]
      expected <- if (stats::sd(v) == 0) NA_real_ else
        stats::cor(v, centroids[g, cl][ok], method = "spearman")
      expect_equal(out[i, cl], expected, tolerance = 1e-12)
    }
  }
})
