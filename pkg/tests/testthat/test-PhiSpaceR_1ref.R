sim_sce <- function(genes, n_cells, seed, prefix) {
  set.seed(seed)
  X <- Matrix::rsparsematrix(length(genes), n_cells, density = 0.3,
                             rand.x = function(m) stats::rpois(m, 3) + 1)
  dimnames(X) <- list(genes, paste0(prefix, seq_len(n_cells)))
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = X, log1p = log1p(X))
  )
  sce$type <- sample(c("A", "B", "C"), n_cells, replace = TRUE)
  sce
}

genes <- paste0("g", 1:60)
ref <- sim_sce(genes, 80, 1, "r")
# The query lacks g1-g5, has extra genes and lists its genes in another order
qry <- sim_sce(c(rev(genes[6:60]), paste0("x", 1:5)), 50, 2, "q")
shared <- genes[6:60]
# The same query reduced to the shared genes, in reference order, one assay
qry_min <- SingleCellExperiment::SingleCellExperiment(
  list(log1p = SummarizedExperiment::assay(qry, "log1p")[shared, ])
)

run <- function(query, ...) {
  PhiSpaceR_1ref(ref, query, phenotypes = "type", refAssay = "log1p",
                 ncomp = 3, ...)
}

test_that("PhiSpaceR_1ref uses only the shared genes, whatever their order", {
  for (args in list(list(), list(nfeat = 10),
                    list(selectedFeat = c("g1", genes[seq(6, 60, by = 3)])))) {
    a <- do.call(run, c(list(qry), args))
    b <- do.call(run, c(list(qry_min), args))
    expect_equal(a$PhiSpaceScore, b$PhiSpaceScore)
    expect_equal(a$YrefHat, b$YrefHat)
    expect_true(all(a$selectedFeat %in% shared))
  }
})

test_that("PhiSpaceR_1ref matches mvr() and phenotype() on the shared genes", {
  sel <- genes[seq(6, 60, by = 2)]
  res <- run(qry, selectedFeat = sel)
  Xref <- Matrix::t(SummarizedExperiment::assay(ref, "log1p")[sel, ])
  fit <- mvr(Xref, codeY(ref, "type"), ncomp = 3, method = "PLS")
  atlas <- list(ncomp = 3, selectedFeat = sel, center = TRUE, scale = FALSE,
                reg_re = fit)
  Xq <- Matrix::t(SummarizedExperiment::assay(qry, "log1p")[sel, ])
  expect_equal(res$PhiSpaceScore, phenotype(Xq, atlas, assayName = "log1p"))
  expect_equal(res$YrefHat, phenotype(Xref, atlas, assayName = "log1p"))
})

test_that("PhiSpaceR_1ref scores a list of queries as separate runs", {
  qry2 <- sim_sce(genes, 30, 3, "s")
  res <- run(list(qry, qry2), nfeat = 10)
  expect_length(res$PhiSpaceScore, 2)
  # Genes shared by all queries define the model, so compare with a run on
  # the second query restricted to those genes
  single <- run(qry2[shared, ], nfeat = 10)
  expect_equal(res$PhiSpaceScore[[2]], single$PhiSpaceScore)
})

test_that("PhiSpaceR_1ref gives the same scores on the rank assay", {
  ranked <- function(sce) RankTransf(sce, "counts")
  a <- PhiSpaceR_1ref(ranked(ref), ranked(qry), phenotypes = "type",
                      refAssay = "rank", ncomp = 3, nfeat = 10)
  b <- PhiSpaceR_1ref(ranked(ref)[shared, ], ranked(qry)[shared, ],
                      phenotypes = "type", refAssay = "rank", ncomp = 3,
                      nfeat = 10)
  expect_equal(a$PhiSpaceScore, b$PhiSpaceScore)
})
