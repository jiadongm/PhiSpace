#' Convert PhiSpace Scores to Probabilities
#'
#' Converts a cell by cell type likelihood matrix to a probability matrix by
#' setting negative values to zero and normalizing each row to sum to 1.
#'
#' @param PhiSc A numeric matrix containing values ranging from -1 to 1.
#'   Rows typically represent cells and columns represent cell types.
#'
#' @return A numeric matrix of the same dimensions as the input, where:
#'   \itemize{
#'     \item All negative values are set to zero
#'     \item Each row sums to 1 (forms a probability distribution)
#'     \item Row and column names are preserved from the input
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets all negative similarity values to zero
#'   \item Calculates the sum of each row
#'   \item Divides each element in a row by the row sum to create probabilities
#'   \item Handles edge cases where rows sum to zero after removing negatives
#' }
#'
#' For rows that contain all negative values (resulting in zero sums after transformation),
#' the function will issue a warning and set those rows to all zeros. Alternative behavior
#' could be implemented to assign uniform probabilities.
#'
#'
#' @export
#'
Score2Prob <- function(PhiSc) {
  # Check if input is a matrix
  if (!is.matrix(PhiSc)) {
    PhiSc <- as.matrix(PhiSc)
  }

  # Check if input is numeric
  if (!is.numeric(PhiSc)) {
    stop("Input must be a numeric matrix")
  }

  # Set negative values to zero
  PhiSc[PhiSc < 0] <- 0

  # Calculate row sums
  row_sums <- rowSums(PhiSc)

  # Handle rows with all zeros (avoid division by zero)
  zero_rows <- row_sums == 0
  if (any(zero_rows)) {
    warning(paste(sum(zero_rows), "row(s) have all zero or negative values.",
                  "These rows will be set to zero probabilities."))
  }

  # Divide each row by its sum to create probabilities
  # Use sweep for efficient row-wise division
  prob_matrix <- sweep(PhiSc, 1, row_sums, FUN = "/")

  # Set NaN values to 0 (from 0/0 division)
  prob_matrix[is.nan(prob_matrix)] <- 0

  # Preserve row and column names
  dimnames(prob_matrix) <- dimnames(PhiSc)

  return(prob_matrix)
}
