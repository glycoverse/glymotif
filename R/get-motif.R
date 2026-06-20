# `glygen_motifs` is a tibble with the following columns:
# - `accession`: the GlyGen accession number of the motif
# - `name`: the name of the motif
# - `alignment`: the alignment of the motif
# - `glycan_structure`: the glycan structure (glyrepr::glycan_structure()) of the motif

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


#' Get the Structures or Alignments of Known Motifs
#'
#' Given a character vector of motifs names in GlycoMotif GlyGen Collection,
#' these functions return the structures or alignments of the motifs.
#'
#' @param name A character vector of the motif name.
#' @returns
#'   - `get_motif_structure()`: a [glyrepr::glycan_structure()]
#'   - `get_motif_alignment()`: a character vector of motif alignments.
#'
#' For all three functions, if `name` has length greater than 1,
#' the return value is named with the motif names.
#'
#' @examples
#' get_motif_structure("N-Glycan core basic")
#' get_motif_alignment("N-Glycan core basic")
#'
#' get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' @seealso [db_motifs()], [is_known_motif()]
#'
#' @export
get_motif_structure <- function(name) {
  check_names(name)
  glygen_motifs$glycan_structure[match(name, glygen_motifs$name)]
}


#' @rdname get_motif_structure
#' @export
get_motif_alignment <- function(name) {
  check_names(name)
  res <- glygen_motifs$alignment[match(name, glygen_motifs$name)]
  if (length(res) > 1) {
    res <- rlang::set_names(res, name)
  }
  res
}


check_names <- function(name) {
  checkmate::assert_character(name)
  if (!all(name %in% glygen_motifs$name)) {
    unknown_names <- name[!name %in% glygen_motifs$name]
    cli::cli_abort("Unknown motif: {.val {unknown_names}}.")
  }
}
