test_that("PhiSpace uses scoreCells fallback for one-class references", {
  fixture <- make_scorecells_fixture()
  fixture$reference$celltype <- "tumor"

  out <- PhiSpace(
    reference = fixture$reference,
    query = fixture$query,
    phenotypes = "celltype",
    refAssay = "logcounts",
    queryAssay = "logcounts",
    fallback = "scoreCells",
    fallback_score = "correlation"
  )

  expect_true("PhiSpace" %in% SingleCellExperiment::reducedDimNames(out))
  phi <- SingleCellExperiment::reducedDim(out, "PhiSpace")
  expect_equal(dim(phi), c(ncol(fixture$query), 1))
  expect_equal(colnames(phi), "tumor")
})

test_that("PhiSpace fallback preserves list-query shape", {
  fixture <- make_scorecells_fixture()
  fixture$reference$celltype <- "tumor"

  out <- PhiSpace(
    reference = fixture$reference,
    query = list(q1 = fixture$query, q2 = fixture$query),
    phenotypes = "celltype",
    refAssay = "logcounts",
    queryAssay = "logcounts",
    fallback = "scoreCells",
    fallback_score = "signature"
  )

  expect_equal(length(out), 2)
  expect_true(all(vapply(
    out,
    function(x) "PhiSpace" %in% SingleCellExperiment::reducedDimNames(x),
    logical(1)
  )))
  expect_equal(
    dim(SingleCellExperiment::reducedDim(out[[1]], "PhiSpace")),
    c(ncol(fixture$query), 1)
  )
})

test_that("PhiSpace updateRef is blocked when fallback is used", {
  fixture <- make_scorecells_fixture()
  fixture$reference$celltype <- "tumor"

  expect_error(
    PhiSpace(
      reference = fixture$reference,
      query = fixture$query,
      phenotypes = "celltype",
      refAssay = "logcounts",
      queryAssay = "logcounts",
      updateRef = TRUE,
      fallback = "scoreCells"
    ),
    "YrefHat and YrefHatNorm are not computed"
  )
})
