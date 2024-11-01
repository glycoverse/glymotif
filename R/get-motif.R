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


#' Get a Motif Graph
#'
#' This function returns a motif graph for a given motif name.
#' It is handy if the motif is a known motif in GlycoMotif GlyGen Collection.
#' To create custom motif graphs, use [glyparse::parse_iupac_condensed()].
#'
#' @param name A character string of the motif name.
#' @return A `glycan_graph` object.
#'
#' @examples
#' get_motif("N-Glycan core basic")
#' get_motif("O-Glycan core 1")
#'
#' @export
get_motif <- function(name) {
  if (!name %in% glygen_motifs$name) {
    cli::cli_abort("Unknown motif: {.val {name}}.")
  }
  iupac <- glygen_motifs$iupac[glygen_motifs$name == name]
  graph <- glyparse::parse_iupac_condensed(iupac)
  aglycon <- glygen_motifs$aglycon[glygen_motifs$name == name]
  alignment <- glygen_motifs$alignment[glygen_motifs$name == name]
  list(graph = graph, aglycon = aglycon, alignment = alignment)
}
