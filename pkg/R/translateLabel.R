#' Translate cell type labels from finer to coarser according to a dictionary.
#'
#' @param text
#' @param dictionary
#' @param dFrom
#' @param dTo
#'
#' @return A vector containing translated labels.
#'
#' @export
translateLabel <- function(text, dictionary, dFrom = "label.fine", dTo = "label.main"){
  translated <- plyr::mapvalues(text,
                                from = dictionary[,dFrom],
                                to = dictionary[,dTo], warn_missing = F)
  return(translated)
}
