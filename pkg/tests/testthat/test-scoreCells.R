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
