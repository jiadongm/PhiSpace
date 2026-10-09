make_pb_fixture <- function() {
  counts <- matrix(c(1, 0, 2,
                     3, 5, 0,
                     0, 4, 7), nrow = 3, byrow = TRUE,
                   dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3)))
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::Matrix(counts, sparse = TRUE))
  )
  sce$type <- c("A", "B", "C")
  sce
}

test_that("pseudoBulk sums nPool draws of a cluster's cells", {
  # Two identical cells per cluster, so every pseudo-bulk is nPool copies of
  # that profile
  sce <- make_pb_fixture()
  X <- as.matrix(SummarizedExperiment::assay(sce, "counts"))
  sce <- cbind(sce, sce)
  colnames(sce) <- paste0("cell", 1:6)
  pb <- pseudoBulk(sce, phenotypes = "type", resampSizes = 2, nPool = 4)
  expect_true(is.matrix(SummarizedExperiment::assay(pb, "counts")))
  expect_equal(unname(SummarizedExperiment::assay(pb, "counts")),
               unname(4 * X[, rep(1:3, each = 2)]))
  expect_equal(rownames(pb), rownames(sce))
  expect_equal(colnames(pb), paste0("PB", 1:6))
  expect_equal(pb$type, rep(c("A", "B", "C"), each = 2))
  Y <- SingleCellExperiment::reducedDim(pb, "response")
  expect_equal(colnames(Y), c("A", "B", "C"))
  # codeY() codes membership as 1 and non-membership as -1
  expect_equal(unname(Y), unname(2 * diag(3)[rep(1:3, each = 2), ] - 1))
})

test_that("pseudoBulk with calcMean divides the sums by nPool", {
  sce <- make_pb_fixture()
  sce <- cbind(sce, sce, sce)
  sums <- pseudoBulk(sce, phenotypes = "type", resampSizes = 3, nPool = 5)
  means <- pseudoBulk(sce, phenotypes = "type", resampSizes = 3, nPool = 5,
                      calcMean = TRUE)
  expect_equal(SummarizedExperiment::assay(means, "counts"),
               SummarizedExperiment::assay(sums, "counts") / 5)
})
