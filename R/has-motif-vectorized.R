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

  if (missing(motifs)) {
    motifs <- available_motifs()
    motif_type <- "known"
  } else {
    motif_type <- get_motifs_type(motifs)
  }

  if (motif_type == "known" && alignment_provided) {
    rlang::warn("Use user-provided alignments, not the ones in the database.")
  }

  suppressWarnings({
    if (alignment_provided) {
      res <- purrr::map2_lgl(motifs, alignments, ~ try_has_motif(
        glycan, .x, alignment = .y, ignore_linkages = ignore_linkages
      ))
    } else {
      res <- purrr::map_lgl(motifs, ~ try_has_motif(
        glycan, .x, ignore_linkages = ignore_linkages
      ))
    }
  }, classes = "warning_custom_alignment"
  )

  if (any(is.na(res))) {
    cli::cli_abort("Motifs at indices {.val {which(is.na(res))}} are neither known motif names or able to be parsed as IUPAC-condensed structure strings.")
  }
  if (is.null(names(res)) && is.character(motifs)) {
    names(res) <- motifs
  }
  res
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
  motif_type <- get_motif_type(motif)

  if (motif_type == "known" && alignment_provided) {
    db_alignment <- glygen_motifs$alignment[glygen_motifs$name == motif]
    if (alignment != db_alignment) {
      cli::cli_warn("The provided alignment type {.val {alignment}} is different from the motif's alignment type {.val {db_alignment}} in database.")
    }
  }

  suppressWarnings({
    if (alignment_provided) {
      res <- purrr::map_lgl(
        glycans, try_has_motif, motif, alignment = alignment, ignore_linkages = ignore_linkages
      )
    } else {
      res <- purrr::map_lgl(
        glycans, try_has_motif, motif, ignore_linkages = ignore_linkages
      )
    }
  }, classes = "warning_custom_alignment"
  )

  if (any(is.na(res))) {
    cli::cli_abort("Glycans at indices {.val {which(is.na(res))}} are not able to be parsed as IUPAC-condensed structure strings.")
  }
  if (is.null(names(res)) && is.character(glycans)) {
    names(res) <- glycans
  }
  res
}


get_motif_type <- function(motif) {
  if (glyrepr::is_glycan(motif)) {
    return("glycan_graph")
  } else if (is.character(motif)) {
    if (is_known_motif(motif)) {
      return("known")
    } else {
      return("iupac")
    }
  } else {
    rlang::abort("`motif` must be either a 'glycan_graph' object, an IUPAC-condensed structure string, or a known motif name.")
  }
}


try_has_motif <- function(glycan, motif, ...) {
  tryCatch(has_motif(glycan, motif, ...), error = function(e) NA)
}
