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


#' Get a Motif Structure Graph, Alignment, or Aglycon
#'
#' Given a motif name in the GlycoMotif GlyGen Collection,
#' `get_motif_graph()` returns a `glycan_graph` object of the motif structure,
#' `get_motif_alignment()` returns the alignment of the motif, and
#' `get_motif_aglycon()` returns the aglycon of the motif.
#' For available motifs, use [availabel_motifs()].
#'
#' @param name A character string of the motif name.
#' @return `get_motif_graph()` returns a `glycan_graph` object,
#' `get_motif_alignment()` returns a character string, and
#' `get_motif_aglycon()` returns a character string.
#'
#' @examples
#' get_motif_graph("N-Glycan core basic")
#' get_motif_alignment("N-Glycan core basic")
#' get_motif_aglycon("N-Glycan core basic")
#'
#' @seealso [availabel_motifs()], [is_known_motif()]
#'
#' @export
get_motif_graph <- function(name) {
  if (!name %in% glygen_motifs$name) {
    cli::cli_abort("Unknown motif: {.val {name}}.")
  }
  iupac <- glygen_motifs$iupac[glygen_motifs$name == name]
  graph <- glyparse::parse_iupac_condensed(iupac)
  graph
}


#' @rdname get_motif_graph
#' @export
get_motif_alignment <- function(name) {
  if (!name %in% glygen_motifs$name) {
    cli::cli_abort("Unknown motif: {.val {name}}.")
  }
  glygen_motifs$alignment[glygen_motifs$name == name]
}


#' @rdname get_motif_graph
#' @export
get_motif_aglycon <- function(name) {
  if (!name %in% glygen_motifs$name) {
    cli::cli_abort("Unknown motif: {.val {name}}.")
  }
  glygen_motifs$aglycon[glygen_motifs$name == name]
}
