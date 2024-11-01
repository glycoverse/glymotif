#' Get Available Motifs
#'
#' This function returns the names of the available motifs available in the package.
#'
#' @return A character vector of motif names.
#' @export
available_motifs <- function() {
  glygen_motifs$name
}
