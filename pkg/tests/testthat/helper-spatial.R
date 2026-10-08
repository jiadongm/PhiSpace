# Small SpatialExperiment with three spatial groups of 100 cells. Each group sits
# in its own region and has a high score for its own cell type, so clustering
# the scores should recover the groups.
make_spatial_fixture <- function(n_per_group = 100, seed = 1) {
  set.seed(seed)
  n <- 3 * n_per_group
  group <- rep(c("g1", "g2", "g3"), each = n_per_group)
  centres <- rbind(g1 = c(0, 0), g2 = c(10, 0), g3 = c(5, 10))
  coords <- centres[group, ] + matrix(rnorm(2 * n), n, 2)
  colnames(coords) <- c("x", "y")
  cells <- paste0("cell", seq_len(n))
  rownames(coords) <- cells

  scores <- matrix(rnorm(n * 20, sd = 0.1), n, 20,
                   dimnames = list(cells, paste0("type", 1:20)))
  for (i in 1:3) scores[group == paste0("g", i), i] <- scores[group == paste0("g", i), i] + 1

  counts <- matrix(rpois(20 * n, 2), 20, n, dimnames = list(paste0("gene", 1:20), cells))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts, logcounts = log1p(counts)),
    colData = S4Vectors::DataFrame(group = group, row.names = cells),
    spatialCoords = coords
  )
  SingleCellExperiment::reducedDim(spe, "PhiSpace") <- scores
  spe
}

# Run expr and return its value, standard output and messages separately.
capture_all <- function(expr) {
  msgs <- character(0)
  out <- utils::capture.output(
    value <- withCallingHandlers(expr, message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  )
  list(value = value, output = out, messages = msgs)
}

# TRUE when every level of a maps to exactly one level of b and vice versa.
same_partition <- function(a, b) {
  tab <- table(a, b)
  all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
}
