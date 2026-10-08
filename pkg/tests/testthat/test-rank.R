dense_ranks <- function(X) {
  X <- as.matrix(X)
  (apply(X, 2, rank, ties.method = "min") - 1)/(nrow(X) - 1)
}

sim_counts <- function(n_genes = 50, n_cells = 40, seed = 1) {
  set.seed(seed)
  X <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3,
                             rand.x = function(m) stats::rpois(m, 2) + 1)
  dimnames(X) <- list(paste0("gene", seq_len(n_genes)), paste0("cell", seq_len(n_cells)))
  X
}

test_that("RTassay ranks within columns, with ties and explicit zeros", {
  X <- sim_counts()
  X[, 3] <- 0                                   # empty column
  X@x[c(1, 5)] <- 0                             # stored (explicit) zeros
  out <- RTassay(X)
  expect_s4_class(out, "dgCMatrix")
  expect_equal(as.matrix(out), dense_ranks(X), ignore_attr = TRUE)
  expect_true(all(out[as.matrix(X) == 0] == 0))
})

test_that("RTassay gives the same ranks for dense input and negative values", {
  X <- sim_counts()
  expect_equal(as.matrix(RTassay(as.matrix(X))), dense_ranks(X), ignore_attr = TRUE)
  Xneg <- X
  Xneg@x[1:10] <- -Xneg@x[1:10]
  expect_equal(as.matrix(RTassay(Xneg)), dense_ranks(Xneg), ignore_attr = TRUE)
})

test_that(".rankWithinCells ranks the genes within each cell", {
  X <- sim_counts()                             # genes x cells
  cells_by_genes <- Matrix::t(X)
  out <- .rankWithinCells(cells_by_genes)
  expect_equal(dim(out), dim(cells_by_genes))
  expect_equal(dimnames(out), dimnames(cells_by_genes))
  expect_equal(as.matrix(out), t(dense_ranks(X)), ignore_attr = TRUE)

  # A cell's values do not depend on the other cells
  sub <- .rankWithinCells(cells_by_genes[1:5, ])
  expect_equal(as.matrix(sub), as.matrix(out[1:5, ]))
})

test_that(".rankWithinCells matches RankTransf on the selected genes", {
  X <- sim_counts()
  selected <- rownames(X)[seq(1, 50, by = 3)]
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = X[selected, ]))
  ranked <- RankTransf(sce, "counts")
  expect_equal(as.matrix(.rankWithinCells(Matrix::t(X)[, selected])),
               t(as.matrix(SummarizedExperiment::assay(ranked, "rank"))))
})
