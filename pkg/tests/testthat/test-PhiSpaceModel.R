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

test_that("stored statistics are the centred and scaled cross-products", {
  X <- t(as.matrix(SummarizedExperiment::assay(ref, "log1p")))
  Y <- codeY(ref, "type")
  for (sc in c(FALSE, TRUE)) {
    model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                           ncomp = 3, scale = sc)
    Xc <- scale(X, center = TRUE, scale = sc)
    expect_equal(model$stats$XtX, crossprod(Xc), tolerance = 1e-12,
                 ignore_attr = TRUE)
    expect_equal(model$stats$XtY, crossprod(Xc, Y), tolerance = 1e-12,
                 ignore_attr = TRUE)
    expect_identical(rownames(model$stats$XtX), genes)
  }
})

test_that("a model refitted on the query genes equals a model trained on them", {
  for (args in list(list(center = TRUE, scale = FALSE),
                    list(center = TRUE, scale = TRUE),
                    list(center = FALSE, scale = FALSE),
                    list(center = FALSE, scale = TRUE),
                    list(nfeat = 10),
                    list(regMethod = "PCA", ncomp = 30),
                    list(regMethod = "PCA", ncomp = 30, scale = TRUE))) {
    if (is.null(args$ncomp)) args$ncomp <- 3
    full <- do.call(trainPhiSpace,
                    c(list(ref, phenotypes = "type", refAssay = "log1p"), args))
    keep <- full$selectedFeat[full$selectedFeat %in% shared]
    args$nfeat <- NULL
    direct <- do.call(trainPhiSpace,
                      c(list(ref, phenotypes = "type", refAssay = "log1p",
                             selectedFeat = keep), args))
    expect_message(a <- predict(full, qry), "refitting the model on the")
    expect_equal(a, predict(direct, qry), tolerance = 1e-10)
    refit <- .modelForGenes(full, rownames(qry)) |> suppressMessages()
    expect_equal(refit$atlas_re$reg_re$coefficients,
                 direct$atlas_re$reg_re$coefficients, tolerance = 1e-10)
    expect_equal(refit$stats, direct$stats, tolerance = 1e-10)
  }
})

test_that("a model trained on all genes gives the per-query fit", {
  full <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p", ncomp = 3)
  one <- PhiSpaceR_1ref(ref, qry, phenotypes = "type", refAssay = "log1p",
                        ncomp = 3)
  expect_message(a <- predict(full, qry), "lacks 5 of the 60 model genes")
  expect_equal(a, one$PhiSpaceScore, tolerance = 1e-10)
})

test_that("models without statistics keep the shared coefficients, with a warning", {
  rref <- RankTransf(ref, "counts")
  rqry <- RankTransf(qry, "counts")
  rank_model <- trainPhiSpace(rref, phenotypes = "type", refAssay = "rank",
                              ncomp = 3)
  expect_null(rank_model$stats)
  expect_warning(a <- predict(rank_model, rqry), "\"rank\" model cannot be refitted")
  B <- .coefSlice(rank_model$atlas_re$reg_re$coefficients, 3)
  ar <- rank_model$atlas_re
  ar$selectedFeat <- shared
  ar$reg_re$coefficients <- B[shared, ]
  expect_equal(a, phenotype(t(SummarizedExperiment::assay(rqry, "rank")), ar,
                            "rank"))

  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                         ncomp = 3, keepStats = FALSE)
  expect_null(model$stats)
  lost <- 1 - sum(.coefSlice(model$atlas_re$reg_re$coefficients, 3)[shared, ]^2) /
    sum(.coefSlice(model$atlas_re$reg_re$coefficients, 3)^2)
  expect_warning(predict(model, qry),
                 paste0("no stored statistics.*hold ",
                        format(100 * lost, digits = 2), "%"))
  expect_output(print(model), "approximate")
})

test_that("predict stops when too few genes are shared", {
  model <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p", ncomp = 3)
  expect_error(predict(model, qry[paste0("x", 1:5), ]), "shares none")
  expect_error(suppressMessages(predict(model, qry[c("g6", "g7"), ])),
               "fewer than the 3 components")
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

test_that("PhiSpace refits a model once on the genes of all queries", {
  full <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p", ncomp = 3)
  qry2 <- qry[setdiff(rownames(qry), "g6"), 1:20]
  common <- shared[-1]
  direct <- trainPhiSpace(ref, phenotypes = "type", refAssay = "log1p",
                          ncomp = 3, selectedFeat = common)
  rd <- SingleCellExperiment::reducedDim
  expect_message(b <- PhiSpace(full, list(qry, qry2), storeUnNorm = TRUE),
                 "lacks 6 of the 60")
  expect_equal(rd(b[[1]], "PhiSpace_nonNorm"), predict(direct, qry),
               tolerance = 1e-10)
  expect_equal(rd(b[[2]], "PhiSpace_nonNorm"), predict(direct, qry2),
               tolerance = 1e-10)
})
