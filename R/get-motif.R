# `glygen_motifs` is a tibble with the following columns:
# - `accession`: the GlyGen accession number of the motif
# - `name`: the name of the motif
# - `alignment`: the alignment of the motif
# - `glycan_structure`: the glycan structure (glyrepr::glycan_structure()) of the motif

#' Get All Motifs from the Database
#'
#' This function returns all motifs available in the package as a glycan
#' structure vector.
#' We use GlycoMotif GlyGen Collection (https://glycomotif.glyomics.org/glycomotif/GGM) as the source of the motifs.
#' This function is useful to be integrated with [have_motifs()] and [count_motifs()].
#' For example, use `have_motifs(glycans, db_motifs())` to check against all motifs.
#'
#' The returned vector carries database `alignment` and `name` metadata so
#' motif-matching functions can use GlycoMotif alignment rules directly.
#'
#' @return A [glyrepr::glycan_structure()] vector.
#'
#' @examples
#' db_motifs()[1:5]
#'
#' @export
db_motifs <- function() {
  new_db_motifs(
    glygen_motifs$glycan_structure,
    alignment = glygen_motifs$alignment,
    name = glygen_motifs$name
  )
}


#' Create a Database Motif Vector
#'
#' Adds glymotif database metadata to a glycan structure vector.
#'
#' @param x A [glyrepr::glycan_structure()] vector.
#' @param alignment A character vector of database alignment labels.
#' @param name A character vector of database motif names.
#'
#' @return A `glymotif_db_motifs` vector.
#' @noRd
new_db_motifs <- function(x, alignment, name) {
  checkmate::assert_true(glyrepr::is_glycan_structure(x))
  checkmate::assert_character(alignment, len = length(x))
  checkmate::assert_character(name, len = length(x), any.missing = TRUE)

  vctrs::new_vctr(
    vctrs::vec_data(x),
    graphs = attr(x, "graphs", exact = TRUE),
    alignment = alignment,
    name = name,
    class = c("glymotif_db_motifs", "glyrepr_structure")
  )
}


#' Test if Motifs Came from the Glymotif Database
#'
#' @param x An object.
#'
#' @return `TRUE` when `x` is a glymotif database motif vector.
#' @noRd
is_db_motifs <- function(x) {
  inherits(x, "glymotif_db_motifs")
}


#' Get Database Motif Alignments
#'
#' @param x A `glymotif_db_motifs` vector.
#'
#' @return A character vector of alignment labels.
#' @noRd
db_motif_alignments <- function(x) {
  attr(x, "alignment", exact = TRUE)
}


#' Get Database Motif Names
#'
#' @param x A `glymotif_db_motifs` vector.
#'
#' @return A character vector of motif names.
#' @noRd
db_motif_names <- function(x) {
  attr(x, "name", exact = TRUE)
}


#' Remove the Database Motif Subclass
#'
#' @param x A `glymotif_db_motifs` vector.
#'
#' @return A plain [glyrepr::glycan_structure()] vector.
#' @noRd
strip_db_motifs <- function(x) {
  attr(x, "alignment") <- NULL
  attr(x, "name") <- NULL
  class(x) <- setdiff(class(x), "glymotif_db_motifs")
  x
}


#' Subset Database Motif Vectors
#'
#' Keeps database metadata aligned with the returned motif structures.
#'
#' @param x A `glymotif_db_motifs` vector.
#' @param i Subset index.
#' @param ... Additional arguments passed to the underlying vector method.
#'
#' @return A `glymotif_db_motifs` vector.
#' @export
`[.glymotif_db_motifs` <- function(x, i, ...) {
  if (missing(i)) {
    return(x)
  }

  # Use the subset index directly instead of vctrs::vec_restore().
  # Different database motifs can share the same structure string but have
  # different metadata, so restoring from values cannot recover the right row.
  new_db_motifs(
    strip_db_motifs(x)[i, ...],
    alignment = db_motif_alignments(x)[i],
    name = db_motif_names(x)[i]
  )
}


#' Coerce Database Motif Vectors to Character
#'
#' Delegates to the underlying glycan structure character representation.
#'
#' @param x A `glymotif_db_motifs` vector.
#' @param ... Additional arguments passed to [as.character()].
#'
#' @return A character vector.
#' @export
as.character.glymotif_db_motifs <- function(x, ...) {
  as.character(strip_db_motifs(x), ...)
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
