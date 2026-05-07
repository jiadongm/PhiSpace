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
