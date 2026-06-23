# `glygen_motifs` is a tibble with the following columns:
# - `source_id`: the collection identifier of the motif
# - `source`: the collection name of the motif
# - `accession`: the accession number of the motif
# - `name`: the name of the motif
# - `alignment`: the alignment of the motif
# - `glycan_structure`: the glycan structure (glyrepr::glycan_structure()) of the motif

#' Get All Motifs from the Database
#'
#' This function returns a database motif specification.
#' We use GlycoMotif collections (https://glycomotif.glyomics.org/glycomotif/)
#' as the source of the motifs.
#' This function is useful to be integrated with [have_motifs()] and [count_motifs()].
#' For example, use `have_motifs(glycans, db_motifs())` to check against the
#' default GlyGen motif collection, or pass `source_id` to use another
#' collection.
#'
#' Use [db_motif_info()] to inspect the motifs included in the database.
#' You can use `dplyr::distinct(db_motif_info(), source_id, source)` to get
#' all available sources.
#'
#' @param source_id A character vector of motif collection identifiers to use.
#'   Defaults to `"GGM"` for backward compatibility.
#'   Use `dplyr::distinct(db_motif_info(), source_id, source)` to get all available sources.
#'   You can use more than one motif collections like `c("GGM", "CCRC")`.
#'   To use all available motifs, use the `"GM"` collection directly.
#' @return A `db_motifs_spec` object.
#'
#' @examples
#' db_motifs()
#'
#' @export
db_motifs <- function(source_id = "GGM") {
  validate_db_motif_source_id(source_id)
  structure(list(source_id = source_id), class = "db_motifs_spec")
}


#' @export
print.db_motifs_spec <- function(x, ...) {
  cli::cli_text("<{.cls db_motifs_spec}>")
  cli::cli_text(
    "This object should be passed to the {.arg motifs} argument of {.fn have_motifs}, {.fn count_motifs}, {.fn match_motifs}, {.fn add_motifs_lgl}, or {.fn add_motifs_int}."
  )
  cli::cli_text(
    "Configuration: uses GlycoMotif database source ID{?s}: {.val {x$source_id}}"
  )
  invisible(x)
}


#' Get Database Motif Information
#'
#' Returns metadata for all motifs available in the package.
#' You can use `dplyr::distinct(db_motif_info(), source_id, source)` to get
#' all available sources.
#'
#' It contains the following columns:
#' - `source_id`: the collection identifier of the motif
#' - `source`: the collection name of the motif
#' - `accession`: the accession number of the motif
#' - `name`: the name of the motif
#' - `alignment`: the alignment of the motif
#' - `glycan_structure`: the glycan structure (glyrepr::glycan_structure()) of the motif
#'
#' @return A tibble.
#'
#' @examples
#' db_motif_info()
#'
#' @export
db_motif_info <- function() {
  tibble::tibble(
    source = glygen_motifs$source,
    source_id = glygen_motifs$source_id,
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
  info <- info[info$source_id %in% spec$source_id, ]
  motifs <- info$glycan_structure
  names(motifs) <- db_motif_labels(info)

  list(
    motifs = motifs,
    alignments = info$alignment,
    match_degree = NULL,
    allow_duplicate_motifs = TRUE
  )
}


#' Validate Database Motif Source IDs
#'
#' @param source_id A character vector of motif collection identifiers.
#' @return `NULL`, invisibly.
#' @noRd
validate_db_motif_source_id <- function(source_id) {
  checkmate::assert_character(source_id, any.missing = FALSE, min.len = 1)
  available_source_ids <- unique(db_motif_info()$source_id)
  unknown_source_ids <- setdiff(source_id, available_source_ids)
  if (length(unknown_source_ids) > 0) {
    cli::cli_abort(c(
      "Unknown motif source{?s}: {.val {unknown_source_ids}}.",
      "i" = "Use {.code dplyr::distinct(db_motif_info(), source_id, source)} to inspect available sources."
    ))
  }
  invisible(NULL)
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


#' Check if a Motif is Known
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `is_known_motif()` was deprecated in glymotif 0.16.0. Use
#' [db_motif_info()] to inspect database motifs instead.
#'
#' @param name A character vector of motif names.
#' @return A logical vector.
#' @examples
#' is_known_motif(c("N-Glycan core basic", "O-Glycan core 1", "unknown"))
#' @export
is_known_motif <- function(name) {
  lifecycle::deprecate_warn("0.16.0", "is_known_motif()", "db_motif_info()")
  checkmate::assert_character(name)
  name %in% db_motif_name_info()$name
}


#' Get the Structures or Alignments of Known Motifs
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `get_motif_structure()` and `get_motif_alignment()` were deprecated in
#' glymotif 0.16.0. Use [db_motif_info()] to inspect database motifs instead.
#'
#' @param name A character vector of motif names.
#' @returns
#' - `get_motif_structure()`: a [glyrepr::glycan_structure()] vector.
#' - `get_motif_alignment()`: a character vector of motif alignments.
#'
#' For `get_motif_alignment()`, if `name` has length greater than 1, the return
#' value is named with the motif names.
#'
#' @examples
#' get_motif_structure("N-Glycan core basic")
#' get_motif_alignment("N-Glycan core basic")
#'
#' get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
#' get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' @seealso [db_motif_info()]
#' @export
get_motif_structure <- function(name) {
  lifecycle::deprecate_warn(
    "0.16.0",
    "get_motif_structure()",
    "db_motif_info()"
  )
  db_motif_structure_by_name(name)
}


#' @rdname get_motif_structure
#' @export
get_motif_alignment <- function(name) {
  lifecycle::deprecate_warn(
    "0.16.0",
    "get_motif_alignment()",
    "db_motif_info()"
  )
  db_motif_alignment_by_name(name)
}


#' Check if Names are Resolvable Database Motif Names
#'
#' Internal predicate for character inputs that refer to named GGM database motifs.
#'
#' @param name A character vector of the motif name.
#' @return A logical vector.
#' @noRd
is_db_motif_name <- function(name) {
  checkmate::assert_character(name)
  name %in% db_motif_name_info()$name
}


#' Suggest Similar Database Motif Names
#'
#' Finds close GGM motif names for diagnostics without changing exact,
#' case-sensitive motif-name resolution.
#'
#' @param name A character vector of unresolved motif names.
#' @param max_distance The maximum edit distance to consider.
#'
#' @return A named list of character vectors.
#' @noRd
suggest_db_motif_names <- function(name, max_distance = NULL) {
  checkmate::assert_character(name)
  candidates <- unique(db_motif_name_info()$name)
  candidates <- candidates[!is.na(candidates) & candidates != ""]

  purrr::map(
    name,
    function(one_name) {
      distances <- utils::adist(one_name, candidates, ignore.case = FALSE)
      distances <- as.integer(distances[1, ])
      threshold <- if (is.null(max_distance)) {
        max(2L, floor(nchar(one_name) * 0.25))
      } else {
        max_distance
      }
      best_distance <- min(distances)

      if (best_distance > threshold) {
        return(character())
      }

      candidates[distances == best_distance]
    }
  ) |>
    rlang::set_names(name)
}


#' Format Similar Database Motif Name Suggestions
#'
#' @param name A character vector of unresolved motif names.
#'
#' @return A named character vector of cli message bullets.
#' @noRd
format_db_motif_name_suggestions <- function(name) {
  suggestions <- suggest_db_motif_names(name)
  suggestions <- suggestions[purrr::map_lgl(suggestions, ~ length(.x) > 0)]

  if (length(suggestions) == 0) {
    return(character())
  }

  formatted <- purrr::imap_chr(
    suggestions,
    function(suggestion, input) {
      suggestion_text <- paste0('"', suggestion, '"', collapse = ", ")

      if (length(suggestions) == 1 && length(suggestion) == 1) {
        return(paste0("Did you mean ", suggestion_text, "?"))
      }

      paste0("For \"", input, "\", did you mean ", suggestion_text, "?")
    }
  )

  rlang::set_names(formatted, rep("i", length(formatted)))
}


#' Guidance for Inspecting Resolvable Database Motif Names
#'
#' @return A named character vector of cli message bullets.
#' @noRd
db_motif_name_info_guidance <- function() {
  c(
    "i" = 'Use `db_motif_info() |> dplyr::filter(source_id == "GGM")` to inspect valid GGM motif names.'
  )
}


#' Get Database Motif Rows Resolvable by Name
#'
#' Only GGM motif names are accepted for legacy character name lookup through
#' the `motif` and `motifs` arguments. Non-GGM collections are available
#' through [db_motifs()] or [db_motif_info()], but their `name` values are
#' labels rather than character inputs for motif resolution.
#'
#' @return A tibble with GGM motif metadata.
#' @noRd
db_motif_name_info <- function() {
  info <- db_motif_info()
  info[info$source_id == "GGM", ]
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
  info <- db_motif_name_info()
  info$glycan_structure[match(name, info$name)]
}


#' Get Database Motif Alignments by Name
#'
#' @param name A character vector of the motif name.
#' @return A character vector of motif alignments.
#' @noRd
db_motif_alignment_by_name <- function(name) {
  check_db_motif_names(name)
  info <- db_motif_name_info()
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
    cli::cli_abort(c(
      "Unknown motif{?s}: {.val {unknown_names}}.",
      format_db_motif_name_suggestions(unknown_names),
      db_motif_name_info_guidance()
    ))
  }
  invisible(NULL)
}
