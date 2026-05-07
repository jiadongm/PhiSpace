test_that("piScore computes signed and unsigned scores", {
  de <- data.frame(
    avg_log2FC = c(1, -2, NA),
    p_val_adj = c(0.01, 0.001, 0.05),
    row.names = c("gene1", "gene2", "gene3")
  )

  signed <- piScore(de, add_column = FALSE)
  expect_equal(signed[["gene1"]], 2)
  expect_equal(signed[["gene2"]], -6)
  expect_true(is.na(signed[["gene3"]]))

  unsigned <- piScore(de, signed = FALSE, add_column = FALSE)
  expect_equal(unsigned[["gene2"]], 6)
})

test_that("piScore detects common DE result columns", {
  scran_like <- data.frame(
    summary.logFC = c(2, 1),
    FDR = c(0.01, 0.1)
  )
  scored <- piScore(scran_like, rank = TRUE)

  expect_named(scored, c("summary.logFC", "FDR", "pi_score", "pi_rank"))
  expect_equal(scored$pi_rank, c(1, 2))

  explicit <- piScore(
    scran_like,
    logfc_col = "summary.logFC",
    pval_col = "FDR",
    add_column = FALSE
  )
  expect_equal(unname(explicit), c(4, 1))
})

test_that("piScore clamps zero p-values and validates p-values", {
  de <- data.frame(avg_log2FC = 1, p_val_adj = 0)
  expect_warning(
    out <- piScore(de, pval_floor = 1e-10, add_column = FALSE),
    "clamped"
  )
  expect_equal(unname(out), 10)

  bad <- data.frame(avg_log2FC = 1, p_val_adj = -0.1)
  expect_error(piScore(bad), "between 0 and 1")
})
