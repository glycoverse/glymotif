#' Check if a Glycan have the Given Motifs
#'
#' @description
#' This is a vectorized versions of [has_motif()].
#' It checks if a glycan has the given motifs, or all known motifs
#' in the GlycoMotif GlyGen Collection if `motifs` is missing.
#'
#' A logical vector is returned.
#' Names of `motifs` are preserved as names of the output logical vectors.
#' If no names are provided, and the input is a character vector, the names
#' will be the vector itself.
#'
#' @param glycan A 'glycan_graph' object, or an IUPAC-condensed structure string.
#' @param motifs A list of 'glycan_graph' objects, a character vector of
#' IUPAC-condensed structure strings, or a character vector of known motif names.
#' @param ... Ignored.
#' @param alignments A character vector of alignment names for the motifs,
#' either of the same length as `motifs`, or a single character string.
#' Could be "substructure", "core", "terminal", or "whole".
#' Default to "substructure".
#' If not provided and `motifs` are known motif names or missing (using all known motifs),
#' the alignments in the GlycoMotif GlyGen Collection will be used.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical vector, indicating if the glycan has the motifs.
#'
#' @examples
#' library(glyparse)
#'
#' # Check if a glycan has the given motifs using names.
#' # Note that the `motifs` vector will become names of the result vector.
#' glycan <- "Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' has_motifs(glycan, c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' # Using IUPAC-condensed structure strings.
#' # Note that names of `motifs` are preserved in the result vector.
#' motifs <- c(M1 = "Gal(b1-3)GalNAc(a1-", M2 = "GlcNAc(b1-4)GlcNAc(a1-")
#' has_motifs(glycan, motifs)
#'
#' # Or using 'glycan_graph' objects.
#' motifs <- list(
#'   M1 = parse_iupac_condensed("Gal(b1-3)GalNAc(a1-"),
#'   M2 = parse_iupac_condensed("GlcNAc(b1-4)GlcNAc(a1-")
#' )
#' has_motifs(glycan, motifs)
#'
#' # If motifs are missing, use all known motifs in the
#' # GlycoMotif GlyGen Collection.
#' has_motifs(glycan)[1:5]
#'
#' @export
has_motifs <- function(glycan, motifs = NULL, ..., alignments = "substructure", ignore_linkages = FALSE) {
  alignment_provided <- !missing(alignments)

  # Check input arguments
  valid_glycan_arg(glycan)
  valid_motifs_arg(motifs)
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  # Get motif type and default motifs
  if (missing(motifs)) {
    motifs <- available_motifs()
    motif_type <- "known"
  } else {
    motif_type <- get_motifs_type(motifs)
  }

  # Decide alignments
  if (motif_type == "known") {
    alignments <- decide_alignments(motifs, alignments, alignment_provided)
  }

  # Ensure motifs are graphs
  glycan <- ensure_glycan_is_graph(glycan)
  motifs <- ensure_motifs_are_graphs(motifs, motif_type)

  purrr::map2_lgl(
    motifs, alignments,
    ~ has_motif_(glycan, .x, alignment = .y, ignore_linkages = ignore_linkages)
  )
}


#' Check if these Glycans have the Given Motif
#'
#' @description
#' This is a vectorized version of [has_motif()].
#' It checks if a list of glycans (or a character vector of IUPAC-condensed structure strings)
#' have the given motif.
#'
#' A logical vector is returned.
#' Names of `glycans` are preserved as names of the output logical vectors.
#' If no names are provided, and the input is a character vector, the names
#' will be the vector itself.
#'
#' @param glycans A list of 'glycan_graph' objects, or a character vector of IUPAC-condensed structure strings.
#' @param motif A 'glycan_graph' object, an IUPAC-condensed structure string, or a known motif name.
#' @param ... Ignored.
#' @param alignment A character string of alignment name for the motif.
#' Could be "substructure", "core", "terminal", or "whole".
#' Default to "substructure".
#' If not provided and `motif` is a known motif name,
#' the alignment in the GlycoMotif GlyGen Collection will be used.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical vector, indicating if the glycans have the motif.
#'
#' @examples
#' glycan <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
#' have_motif(glycan, "Gal")
#'
#' @export
have_motif <- function(glycans, motif, ..., alignment = "substructure", ignore_linkages = FALSE) {
  alignment_provided <- !missing(alignment)

  # Check input arguments
  valid_glycans_arg(glycans)
  valid_motif_arg(motif)
  valid_alignment_arg(alignment)
  valid_ignore_linkages_arg(ignore_linkages)

  # Decide motif type and alignment
  motif_type <- get_motif_type(motif)
  if (motif_type == "known") {
    alignment <- decide_alignment(motif, alignment, alignment_provided)
  }

  # Ensure glycans and motif are graphs
  glycans <- ensure_glycans_are_graphs(glycans)
  motif <- ensure_motif_is_graph(motif, motif_type)

  purrr::map_lgl(glycans, has_motif_, motif, alignment = alignment, ignore_linkages = ignore_linkages)
}


valid_glycans_arg <- function(x) {
  error_msg <- "`glycans` must be either 'glycan_graph' objects or a character vector of IUPAC-condensed structure strings."
  if (is.character(x)) {
    NULL
  } else if (is.list(x)) {
    if (!purrr::every(x, glyrepr::is_glycan)) {
      rlang::abort(error_msg)
    }
  } else {
    rlang::abort(error_msg)
  }
}


valid_motifs_arg <- function(x) {
  error_msg <- "`motifs` must be either 'glycan_graph' objects, a character vector of IUPAC-condensed structure strings, or a character vector of known motif names."
  if (is.null(x) || is.character(x)) {
    NULL
  } else if (is.list(x)) {
    if (!purrr::every(x, glyrepr::is_glycan)) {
      rlang::abort(error_msg)
    }
  } else {
    rlang::abort(error_msg)
  }
}


valid_alignments_arg <- function(x, motifs) {
  if (!is.character(x)) {
    rlang::abort("`alignments` must be a character vector.")
  }
  if (length(x) != 1 && length(x) != length(motifs)) {
    cli::cli_abort(c(
      "`alignments` must be either a single character string or a character vector of the same length as `motifs`.",
      i = "`motif` length: {.val {length(motifs)}}, `alignments` length: {.val {length(x)}}"
    ))
  }
  if (!all(x %in% c("substructure", "core", "terminal", "whole"))) {
    rlang::abort("`alignments` must be one of 'substructure', 'core', 'terminal' or 'whole'.")
  }
}


get_motifs_type <- function(motifs) {
  if (purrr::every(motifs, glyrepr::is_glycan)) {
    return("glycan_graph")
  } else if (all(is.character(motifs))) {
    if (purrr::every(motifs, is_known_motif)) {
      return("known")
    } else if (purrr::some(motifs, is_known_motif)) {
      unknown_motifs <- motifs[!purrr::map_lgl(motifs, is_known_motif)]
      cli::cli_abort("Motifs {.val {unknown_motifs}} are not known motif names.")
    } else {
      # If all motifs are characters but not known motif names,
      # assume they are IUPAC-condensed structure strings.
      return("iupac")
    }
  } else {
    rlang::abort("`motifs` must be either 'glycan_graph' objects, a character vector of IUPAC-condensed structure strings, or a character vector of known motif names.")
  }
}


decide_alignments <- function(motif_names, alignments, alignments_provided) {
  db_alignments <- glygen_motifs$alignment[glygen_motifs$name %in% motif_names]
  if (alignments_provided) {
    rlang::warn("Use user-provided alignments, not the ones in the database.")
    alignments
  } else {
    db_alignments
  }
}


ensure_motifs_are_graphs <- function(motifs, motif_type) {
  # Convert to graphs
  if (motif_type == "glycan_graph") {
    graph_list <- motifs
  } else if (motif_type == "iupac") {
    graph_list <- purrr::map(motifs, try_parse_iupac_condensed)
    if (any(is.na(graph_list))) {
      invalid_indices <- which(is.na(graph_list))
      cli::cli_abort("Motifs at indices {.val {invalid_indices}} are neither known motif names or able to be parsed as IUPAC-condensed structure strings.")
    }
  } else if (motif_type == "known") {
    graph_list <- purrr::map(motifs, ~ get_motif_graph(.x)$graph)
  }

  # Add names
  if (is.null(names(motifs)) && is.character(motifs)) {
    names(graph_list) <- motifs
  }
  graph_list
}


ensure_glycans_are_graphs <- function(glycans) {
  if(is.character(glycans)) {
    graph_list <- purrr::map(glycans, try_parse_iupac_condensed)
    if (any(is.na(graph_list))) {
      invalid_indices <- which(is.na(graph_list))
      cli::cli_abort("Glycans at indices {.val {invalid_indices}} are not able to be parsed as IUPAC-condensed structure strings.")
    }
    if (is.null(names(graph_list))) {
      names(graph_list) <- glycans
    }
  } else {
    graph_list <- glycans
  }
  graph_list
}


try_parse_iupac_condensed <- function(iupac) {
  tryCatch({
    glyparse::parse_iupac_condensed(iupac)
  }, error = function(e) NA)
}
