# `glygen_motifs` is a tibble with the following columns:
# - `accession`: the GlyGen accession number of the motif
# - `name`: the name of the motif
# - `alignment`: the alignment of the motif
# - `glycan_structure`: the glycan structure (glyrepr::glycan_structure()) of the motif

#' Get All Motifs from the Database
#'
#' This function returns a database motif specification.
#' We use GlycoMotif GlyGen Collection (https://glycomotif.glyomics.org/glycomotif/GGM) as the source of the motifs.
#' This function is useful to be integrated with [have_motifs()] and [count_motifs()].
#' For example, use `have_motifs(glycans, db_motifs())` to check against all motifs.
#'
#' Use [db_motif_info()] to inspect the motifs included in the database.
#'
#' @return A `db_motifs_spec` object.
#'
#' @examples
#' db_motifs()
#' db_motif_info()
#'
#' @export
db_motifs <- function() {
  structure(list(), class = "db_motifs_spec")
}


#' @export
print.db_motifs_spec <- function(x, ...) {
  cli::cli_text("<{.cls db_motifs_spec}>")
  cli::cli_text(
    "This object should be passed to the {.arg motifs} argument of {.fn have_motifs}, {.fn count_motifs}, {.fn match_motifs}, {.fn add_motifs_lgl}, or {.fn add_motifs_int}."
  )
  cli::cli_text("Configuration: uses all GlyGen GlycoMotif database motifs")
  invisible(x)
}


#' Get Database Motif Information
#'
#' Returns metadata for all motifs available in the package.
#'
#' @return A tibble with `source`, `source_id`, `accession`, `name`,
#'   `alignment`, and `glycan_structure` columns.
#'
#' @examples
#' db_motif_info()
#'
#' @export
db_motif_info <- function() {
  tibble::tibble(
    source = rep("GlyGen", nrow(glygen_motifs)),
    source_id = rep("GGM", nrow(glygen_motifs)),
    accession = glygen_motifs$accession,
    name = glygen_motifs$name,
    alignment = glygen_motifs$alignment,
    glycan_structure = glygen_motifs$glycan_structure
  )
}


#' @export
resolve_motif_spec.db_motifs_spec <- function(
  glycans,
  spec,
  alignments,
  match_degree,
  strict_sub = TRUE,
  ignore_linkages = FALSE
) {
  info <- db_motif_info()
  motifs <- info$glycan_structure
  names(motifs) <- db_motif_labels(info)

  list(
    motifs = motifs,
    alignments = info$alignment,
    match_degree = NULL,
    allow_duplicate_motifs = TRUE
  )
}


#' Create Database Motif Labels
#'
#' @param info A database motif metadata tibble.
#'
#' @return A character vector of labels.
#' @noRd
db_motif_labels <- function(info) {
  labels <- info$name
  missing_labels <- is.na(labels) | labels == ""
  labels[missing_labels] <- paste0(
    info$source_id[missing_labels],
    ".",
    info$accession[missing_labels]
  )
  labels
}


#' Check if Names are Database Motif Names
#'
#' Internal predicate for character inputs that refer to named database motifs.
#'
#' @param name A character vector of the motif name.
#' @return A logical vector.
#' @noRd
is_db_motif_name <- function(name) {
  checkmate::assert_character(name)
  name %in% db_motif_info()$name
}


#' Get Database Motif Structures or Alignments by Name
#'
#' Internal helpers for resolving legacy character motif-name inputs.
#'
#' @param name A character vector of the motif name.
#' @return A [glyrepr::glycan_structure()] vector.
#' @noRd
db_motif_structure_by_name <- function(name) {
  check_db_motif_names(name)
  info <- db_motif_info()
  info$glycan_structure[match(name, info$name)]
}


#' Get Database Motif Alignments by Name
#'
#' @param name A character vector of the motif name.
#' @return A character vector of motif alignments.
#' @noRd
db_motif_alignment_by_name <- function(name) {
  check_db_motif_names(name)
  info <- db_motif_info()
  res <- info$alignment[match(name, info$name)]
  if (length(res) > 1) {
    res <- rlang::set_names(res, name)
  }
  res
}


#' Check Database Motif Names
#'
#' @param name A character vector of the motif name.
#' @return `NULL`, invisibly.
#' @noRd
check_db_motif_names <- function(name) {
  checkmate::assert_character(name)
  if (!all(is_db_motif_name(name))) {
    unknown_names <- name[!is_db_motif_name(name)]
    cli::cli_abort("Unknown motif: {.val {unknown_names}}.")
  }
  invisible(NULL)
}
