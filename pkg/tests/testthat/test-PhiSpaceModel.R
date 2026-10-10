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

test_that("trainPhiSpace and predict give the PhiSpaceR_1ref scores", {
  for (args in list(list(), list(nfeat = 10),
                    list(selectedFeat = genes[seq(6, 60, by = 3)]))) {
    one <- do.call(PhiSpaceR_1ref,
                   c(list(ref, qry, phenotypes = "type", refAssay = "log1p",
                          ncomp = 3), args))
    model <- do.call(trainPhiSpace,
                     c(list(ref, phenotypes = "type", refAssay = "log1p",
                            ncomp = 3, genes = shared), args))
    expect_s3_class(model, "PhiSpaceModel")
    expect_identical(predict(model, qry), one$PhiSpaceScore)
    expect_identical(model$atlas_re, one$atlas_re)
    expect_identical(model$selectedFeat, one$selectedFeat)
    expect_identical(model$impScores, one$impScores)
    expect_identical(one$model$atlas_re, one$atlas_re)
  }
})

test_that("trainPhiSpace and predict match PhiSpaceR_1ref on the rank assay", {
  rref <- RankTransf(ref, "counts")
  rqry <- RankTransf(qry, "counts")
  one <- PhiSpaceR_1ref(rref, rqry, phenotypes = "type", refAssay = "rank",
                        ncomp = 3, nfeat = 10)
  model <- trainPhiSpace(rref, phenotypes = "type", refAssay = "rank",
                         ncomp = 3, nfeat = 10, genes = shared)
  expect_identical(predict(model, rqry), one$PhiSpaceScore)
})

test_that("a saved model gives the same scores", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                         ncomp = 3, genes = shared, referenceName = "sim",
                         species = "none", geneIdType = "symbol")
  path <- tempfile(fileext = ".rds")
  saveRDS(model, path)
  expect_identical(predict(readRDS(path), qry), predict(model, qry))
  expect_output(print(model), "3 phenotypes, 55 genes, 80 reference cells")
})

test_that("predict accepts a gene by cell matrix and another assay", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                         ncomp = 3, genes = shared)
  X <- SummarizedExperiment::assay(qry, "log1p")
  expect_identical(predict(model, X), predict(model, qry))
  q2 <- qry
  SummarizedExperiment::assay(q2, "lognorm") <- X
  expect_identical(predict(model, q2, assay = "lognorm"), predict(model, qry))
  expect_error(predict(model, qry, assay = "nope"), "not present")
})

test_that("predict stops when the query lacks model genes", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p", ncomp = 3)
  expect_equal(model$selectedFeat, genes)
  expect_error(predict(model, qry), "lacks 5 of the 60 model genes")
})

test_that("trainPhiSpace accepts a response matrix", {
  Y <- codeY(ref, "type")
  a <- trainPhiSpace(ref, response = Y, refAssay = "log1p", ncomp = 3)
  b <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p", ncomp = 3)
  expect_null(a$phenoDict)
  expect_identical(a$atlas_re, b$atlas_re)
})

test_that("predict refuses newer model formats", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                         ncomp = 3, genes = shared)
  model$format_version <- 99L
  expect_error(predict(model, qry), "format version 99")
})

test_that("PhiSpace accepts a model and gives the scores of a reference run", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                         ncomp = 3, nfeat = 10, genes = shared)
  rd <- SingleCellExperiment::reducedDim
  a <- PhiSpace(ref, qry, phenotypes = "type", refAssay = "log1p", ncomp = 3,
                nfeat = 10, storeUnNorm = TRUE)
  b <- PhiSpace(model, qry, storeUnNorm = TRUE)
  expect_identical(rd(b, "PhiSpace"), rd(a, "PhiSpace"))
  expect_identical(rd(b, "PhiSpace_nonNorm"), rd(a, "PhiSpace_nonNorm"))

  # A list of queries
  qry2 <- qry[, 1:20]
  a <- PhiSpace(ref, list(qry, qry2), phenotypes = "type", refAssay = "log1p",
                ncomp = 3, nfeat = 10)
  b <- PhiSpace(model, list(qry, qry2))
  expect_identical(rd(b[[2]], "PhiSpace"), rd(a[[2]], "PhiSpace"))

  expect_error(PhiSpace(model, qry, updateRef = TRUE), "needs the reference cells")
})

test_that("PhiSpace mixes models and references in a list", {
  ref2 <- sim_sce(genes, 70, 4, "t")
  model2 <- trainPhiSpace(ref2, phenotypes = "type", refAssay = "log1p",
                          ncomp = 3, genes = shared)
  rd <- SingleCellExperiment::reducedDim
  a <- PhiSpace(list(one = ref, two = ref2), qry, phenotypes = "type",
                refAssay = "log1p", ncomp = 3)
  b <- PhiSpace(list(one = ref, two = model2), qry, phenotypes = "type",
                refAssay = "log1p", ncomp = 3)
  expect_identical(rd(b, "PhiSpace"), rd(a, "PhiSpace"))
  expect_true("A(two)" %in% colnames(rd(b, "PhiSpace")))
})
