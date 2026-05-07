make_scorecells_fixture <- function() {
  genes <- c("GeneA1", "GeneA2", "GeneB1", "GeneB2", "Rpl1", "ACTB")
  ref_expr <- matrix(
    c(
      8, 7, 8, 1, 1, 1,
      7, 8, 7, 1, 1, 1,
      1, 1, 1, 8, 7, 8,
      1, 1, 1, 7, 8, 7,
      20, 20, 20, 20, 20, 20,
      20, 20, 20, 20, 20, 20
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, paste0("ref", seq_len(6)))
  )
  query_expr <- matrix(
    c(
      8, 1, 4,
      7, 1, 4,
      1, 8, 4,
      1, 7, 4,
      20, 20, 20,
      20, 20, 20
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, paste0("query", seq_len(3)))
  )

  reference <- SingleCellExperiment::SingleCellExperiment(
    list(logcounts = ref_expr),
    colData = S4Vectors::DataFrame(celltype = rep(c("A", "B"), each = 3))
  )
  query <- SingleCellExperiment::SingleCellExperiment(list(logcounts = query_expr))

  list(reference = reference, query = query)
}

test_that("scoreCells writes correlation and signature reducedDims", {
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

  expect_true("scoreCells_correlation" %in% SingleCellExperiment::reducedDimNames(out))
  expect_true("scoreCells_signature" %in% SingleCellExperiment::reducedDimNames(out))

  cor_scores <- SingleCellExperiment::reducedDim(out, "scoreCells_correlation")
  sig_scores <- SingleCellExperiment::reducedDim(out, "scoreCells_signature")

  expect_equal(dim(cor_scores), c(3, 2))
  expect_equal(dim(sig_scores), c(3, 2))
  expect_equal(colnames(cor_scores), c("A", "B"))
  expect_gt(cor_scores["query1", "A"], cor_scores["query1", "B"])
  expect_gt(cor_scores["query2", "B"], cor_scores["query2", "A"])

  signatures <- S4Vectors::metadata(out)$scoreCells$signatures
  expect_true(all(c("GeneA1", "GeneA2") %in% signatures$A))
  expect_false(any(c("Rpl1", "ACTB") %in% unlist(signatures, use.names = FALSE)))
})

test_that("scoreCells handles single-class references", {
  fixture <- make_scorecells_fixture()
  out <- scoreCells(
    fixture$reference,
    fixture$query,
    class_col = NULL,
    signature_method = "findMarkers",
    n_top_genes = 4,
    n_signature_genes = 2,
    scoring = "signature",
    verbose = FALSE
  )

  sig_scores <- SingleCellExperiment::reducedDim(out, "scoreCells_signature")
  expect_equal(dim(sig_scores), c(3, 1))
  expect_equal(colnames(sig_scores), "reference")
  expect_true(isTRUE(stats::sd(sig_scores[, 1]) > 0))
})
