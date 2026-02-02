rlang::on_load({
  structures <- glyparse::parse_iupac_condensed(glygen_motifs$iupac)
  glygen_motifs$graph <- glyrepr::get_structure_graphs(structures)
  glygen_motifs <- tibble::as_tibble(glygen_motifs)
})

# `glygen_motifs` is a tibble with the following columns:
# - `accession`: the GlyGen accession number of the motif
# - `name`: the name of the motif
# - `aglycon`: the aglycon of the motif
# - `alignment`: the alignment of the motif
# - `iupac`: the IUPAC condensed string of the motif
# - `graph`: the structure graph (igraph object) of the motif

#' Get All Motifs from the Database
#'
#' This function returns the names of all motifs available in the package.
#' We use GlycoMotif GlyGen Collection (https://glycomotif.glyomics.org/glycomotif/GGM) as the source of the motifs.
#' This function is useful to be integrated with [have_motifs()] and [count_motifs()].
#' For example, use `have_motifs(glycans, db_motifs())` to check against all motifs.
#'
#' @return A character vector of motif names.
#'
#' @examples
#' db_motifs()[1:5]
#'
#' @export
db_motifs <- function() {
  glygen_motifs$name
}


#' Check if a Motif is Known
#'
#' This function checks if motifs are known motifs in GlycoMotif GlyGen Collection.
#'
#' @param name A character vector of the motif name.
#' @return A logical vector.
#' @examples
#' is_known_motif(c("N-Glycan core basic", "O-Glycan core 1", "unknown"))
#' @export
is_known_motif <- function(name) {
  checkmate::assert_character(name)
  name %in% glygen_motifs$name
}


#' Get the Structures, Alignments, or Aglycons of Known Motifs
#'
#' Given a character vector of motifs names in GlycoMotif GlyGen Collection,
#' these functions return the structures, alignments, or aglycons of the motifs.
#'
#' @param name A character vector of the motif name.
#' @returns `get_motif_structure()` returns a `glycan_graph` when `name` is a character
#' scalar, and a list of `glycan_graph` when `name` is a character vector.
#' `get_motif_alignment()` returns a character vector of motif alignments.
#' `get_motif_aglycon()` returns a character vector of motif aglycons.
#' For all three functions, if `name` has length greater than 1,
#' the return value is named with the motif names.
#'
#' @examples
#' get_motif_structure("N-Glycan core basic")
#' get_motif_alignment("N-Glycan core basic")
#' get_motif_aglycon("N-Glycan core basic")
#'
#' get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_aglycon(c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' @seealso [db_motifs()], [is_known_motif()]
#'
#' @export
get_motif_structure <- function(name) {
  check_names(name)
  res <- glygen_motifs$graph[match(name, glygen_motifs$name)]
  glyrepr::as_glycan_structure(res)
}


#' @rdname get_motif_structure
#' @export
get_motif_alignment <- function(name) {
  check_names(name)
  res <- glygen_motifs$alignment[match(name, glygen_motifs$name)]
  if (length(res) > 1) res <- rlang::set_names(res, name)
  res
}


#' @rdname get_motif_structure
#' @export
get_motif_aglycon <- function(name) {
  check_names(name)
  res <- glygen_motifs$aglycon[match(name, glygen_motifs$name)]
  if (length(res) > 1) res <- rlang::set_names(res, name)
  res
}


check_names <- function(name) {
  checkmate::assert_character(name)
  if (!all(name %in% glygen_motifs$name)) {
    unknown_names <- name[!name %in% glygen_motifs$name]
    cli::cli_abort("Unknown motif: {.val {unknown_names}}.")
  }
}
