#' Get Available Motifs
#'
#' This function returns the names of the available motifs available in the package.
#'
#' @return A character vector of motif names.
#'
#' @examples
#' available_motifs()[1:5]
#'
#' @export
available_motifs <- function() {
  glygen_motifs$name
}


#' Check if a Motif is Known
#'
#' This function checks if a motif is a known motif in GlycoMotif GlyGen Collection.
#'
#' @param name A character string of the motif name.
#' @return A logical value.
#' @export
is_known_motif <- function(name) {
  name %in% glygen_motifs$name
}


#' Get the Structure Graphs, Alignments, or Aglycons of Known Motifs
#'
#' Given a character vector of motifs names in GlycoMotif GlyGen Collection,
#' these functions return the structure graphs, alignments, or aglycons of the motifs.
#'
#' @param name A character vector of the motif name.
#' @returns `get_motif_graph()` returns a `glycan_graph` when `name` is a character
#' scalar, and a list of `glycan_graph` when `name` is a character vector.
#' `get_motif_alignment()` returns a character vector of motif alignments.
#' `get_motif_aglycon()` returns a character vector of motif aglycons.
#' For all three functions, if `name` has length greater than 1,
#' the return value is named with the motif names.
#'
#' @examples
#' get_motif_graph("N-Glycan core basic")
#' get_motif_alignment("N-Glycan core basic")
#' get_motif_aglycon("N-Glycan core basic")
#'
#' get_motif_graph(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_aglycon(c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' @seealso [availabel_motifs()], [is_known_motif()]
#'
#' @export
get_motif_graph <- function(name) {
  check_names(name)
  res <- glygen_motifs$graph[glygen_motifs$name %in% name]
  if (length(res) > 1) {
    rlang::set_names(res, name)
  } else {
    res[[1]]
  }
}


#' @rdname get_motif_graph
#' @export
get_motif_alignment <- function(name) {
  check_names(name)
  res <- glygen_motifs$alignment[glygen_motifs$name %in% name]
  if (length(res) > 1) res <- rlang::set_names(res, name)
  res
}


#' @rdname get_motif_graph
#' @export
get_motif_aglycon <- function(name) {
  check_names(name)
  res <- glygen_motifs$aglycon[glygen_motifs$name %in% name]
  if (length(res) > 1) res <- rlang::set_names(res, name)
  res
}


check_names <- function(name) {
  if (!all(name %in% glygen_motifs$name)) {
    unknown_names <- name[!name %in% glygen_motifs$name]
    cli::cli_abort("Unknown motif: {.val {unknown_names}}.")
  }
}
