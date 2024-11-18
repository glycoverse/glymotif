#' Vectorized Motif Checking Functions
#'
#' @description
#' Vectorized versions of [has_motif()]:
#' - `has_motifs()`: one glycan, many motifs.
#' - `have_motif()`: many glycans, one motif.
#' - `have_motifs()`: many glycans, many motifs.
#'
#' Vectorized versions of [counts_motif()]:
#' - `counts_motifs()`: one glycan, many motifs.
#' - `count_motif()`: many glycans, one motif.
#' - `count_motifs()`: many glycans, many motifs.
#'
#' They can be told apart by English grammar:
#' - Functions starting with "has" and "counts" work with one glycan,
#'   "have" and "count" with multiple glycans.
#' - Functions ending with "motif" work with one motif, "motifs" with multiple motifs.
#'
#' @details
#' # Why not `purrr`?
#'
#' These functions have performance benefits over simply using `purrr` functions
#' on `has_motif()` and `counts_motif()`.
#' These two functions have many internal checks and conversions,
#' which can be redundant using `purrr`.
#' For example, when passing N motifs to `has_motifs()`,
#' the checks and conversions will be performed N + 1 times
#' (N for the motifs and 1 for the glycan).
#' However, using `purrr::map_lgl(motifs, has_motif, glycan)`
#' will perform the checks and conversions N + N times.
#' The glycan will be processed repeatedly for each motif,
#' even thought N - 1 of them are redundant.
#' The situation is even worse if you are dealing with multiple glycans
#' and multiple motifs.
#' Using `purrr` functions will perform the checks and conversions
#' N * M times (N for the motifs and M for the glycans),
#' while using `have_motifs()` only need N + M times.
#' A huge save!
#'
#' # Why so many?
#'
#' You may wonder why we can't just make `has_motif()` and `counts_motif()` vectorized.
#' The reason is the the same between base "apply" family functions and `purrr` functions.
#' These vectorized functions have consistent and predictable return types.
#' - `has_motif()` and `counts_motif()` always return a scalar.
#' - `has_motifs()`, `have_motif()`, `counts_motifs()`, and `count_motif()` always return a vector.
#' - `have_motifs()` and `count_motifs()` always return a matrix.
#'
#' # More about arguments
#'
#' When using known motifs in the GlycoMotif GlyGen Collection,
#' the best practice is to not provide the alignment arguments,
#' and let the function decide the alignment based on the motif name.
#' However, it is still possible to override the default alignments.
#' In this case, the user-provided alignments will be used,
#' but a warning will be issued.
#'
#' Monosaccharide types must be the same within `glycans` or `motifs`.
#' However, an holistic obscurer monosaccharide type in `motifs` than in `glycans` is allowed.
#'
#' @param glycan A 'glycan_graph' object, or an IUPAC-condensed structure string.
#' @param motif A 'glycan_graph' object, an IUPAC-condensed structure string, or a known motif name.
#' @param glycans A list of 'glycan_graph' objects, or a character vector of IUPAC-condensed structure strings.
#' @param motifs A list of 'glycan_graph' objects, a character vector of
#' IUPAC-condensed structure strings, or a character vector of known motif names.
#' If missing, all known motifs in the GlycoMotif GlyGen Collection will be used.
#' @param alignment A single character, one of 'substructure', 'core', 'terminal' or 'whole'.
#' If not provided and `motif` is a known motif name,
#' the alignment in the GlycoMotif GlyGen Collection will be used.
#' @param alignments A character vector of alignment names for the motifs,
#' either of the same length as `motifs`, or a single character string.
#' If not provided and `motifs` are known motif names or missing (using all known motifs),
#' the alignments in the GlycoMotif GlyGen Collection will be used.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#' @param simplify A logical value. Only used in `have_motifs()`.
#' If `TRUE`, drop columns (motifs) with all FALSE values. Default is `FALSE`.
#'
#' @return `has_motifs()`, `have_motif()`, `counts_motifs()`, and `count_motif()`
#' return a logical vector of the same length as the input motifs or glycans, respectively.
#' `have_motifs()` and `count_motifs()` return a matrix.
#' Rows are glycans and columns are motifs.
#'
#' The dimension names are determined by the input arguments.
#' If `glycans` or `motifs` have names, they will be preserved in the result.
#' If not, but `glycans` or `motifs` are character vectors, they will be used as names.
#'
#' @seealso [has_motif()], [counts_motif()]
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
#' # Check if a list of glycans have the given motif.
#' glycans <- c(G1 = "Gal(b1-3)GalNAc(a1-", G2 = "GlcNAc(b1-4)GlcNAc(a1-")
#' have_motif(glycans, "O-Glycan core 1")
#'
#' # Check if a list of glycans have the given motifs.
#' glycans <- c(G1 = "Gal(b1-3)GalNAc(a1-", G2 = "GlcNAc(b1-4)GlcNAc(a1-")
#' motifs <- c(M1 = "Gal(b1-", M2 = "GlcNAc(b1-")
#' have_motifs(glycans, motifs)
#'
#' @export
has_motifs <- function(glycan, motifs = NULL, alignments = NULL, ignore_linkages = FALSE) {
  one_glycan_many_motifs_wrapper(glycan, motifs, alignments, ignore_linkages, simple_has_motifs)
}


#' @rdname has_motifs
#' @export
have_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  .f <- function(glycans, motif, alignment, ignore_linkages) {
    purrr::map_lgl(glycans, has_motif_, motif, alignment = alignment, ignore_linkages = ignore_linkages)
  }
  many_glycans_one_motif_wrapper(glycans, motif, alignment, ignore_linkages, .f)
}


#' @rdname has_motifs
#' @export
have_motifs <- function(glycans, motifs = NULL, alignments = NULL, ignore_linkages = FALSE, simplify = FALSE) {
  .f <- function(glycans, motifs, alignments, ignore_linkages) {
    lgl_list <- purrr::map(
      glycans, simple_has_motifs, motifs,
      alignments = alignments, ignore_linkages = ignore_linkages
    )
    do.call(rbind, lgl_list)
  }
  result <- many_glycans_many_motifs_wrapper(glycans, motifs, alignments, ignore_linkages, .f)
  if (simplify) {
    # Drop columns will all FALSE
    result <- result[, colSums(result) > 0, drop = FALSE]
  }
  result
}


#' @rdname has_motifs
#' @export
counts_motifs <- function(glycan, motifs = NULL, alignments = NULL, ignore_linkages = FALSE) {
  one_glycan_many_motifs_wrapper(glycan, motifs, alignments, ignore_linkages, simple_counts_motifs)
}


#' @rdname has_motifs
#' @export
count_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  .f <- function(glycans, motif, alignment, ignore_linkages) {
    purrr::map_int(glycans, counts_motif_, motif, alignment = alignment, ignore_linkages = ignore_linkages)
  }
  many_glycans_one_motif_wrapper(glycans, motif, alignment, ignore_linkages, .f)
}


#' @rdname has_motifs
#' @export
count_motifs <- function(glycans, motifs = NULL, alignments = NULL, ignore_linkages = FALSE, simplify = FALSE) {
  .f <- function(glycans, motifs, alignments, ignore_linkages) {
    int_list <- purrr::map(
      glycans, simple_counts_motifs, motifs,
      alignments = alignments, ignore_linkages = ignore_linkages
    )
    do.call(rbind, int_list)
  }
  result <- many_glycans_many_motifs_wrapper(glycans, motifs, alignments, ignore_linkages, .f)
  if (simplify) {
    # Drop columns will all FALSE
    result <- result[, colSums(result) > 0, drop = FALSE]
  }
  result
}


# ----- Argument Processing -----
one_glycan_many_motifs_wrapper <- function(glycan, motifs, alignments, ignore_linkages, func) {
  # Check input arguments
  valid_glycan_arg(glycan)
  valid_motifs_arg(motifs)
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  if (!is.null(motifs) && length(motifs) == 0) {
    return(rlang::set_names(logical(0L), character(0L)))
  }

  # Get motif type and default motifs
  if (is.null(motifs)) {
    motifs <- available_motifs()
    motif_type <- "known"
  } else {
    motif_type <- get_motifs_type(motifs)
  }

  # Decide alignments
  if (motif_type == "known") {
    alignments <- decide_alignments(motifs, alignments)
  } else if (is.null(alignments)) {
    alignments <- "substructure"
  }

  # Ensure motifs are graphs
  glycan <- ensure_glycan_is_graph(glycan)
  motifs <- ensure_motifs_are_graphs(motifs, motif_type)

  # Ensure NE graphs
  glycan <- ensure_ne_graph(glycan)
  motifs <- ensure_ne_graphs(motifs)

  # Ensure mono types are the same
  if (motif_type != "known" && !same_mono_types(motifs)) {
    rlang::abort("All motifs must have the same monosaccharide type.")
  }
  glycan <- ensure_glycan_mono_type(glycan, motifs[[1]])

  func(glycan, motifs, alignments, ignore_linkages)
}


many_glycans_one_motif_wrapper <- function(glycans, motif, alignment, ignore_linkages, func) {
  # Check input arguments
  valid_glycans_arg(glycans)
  valid_motif_arg(motif)
  valid_alignment_arg(alignment)
  valid_ignore_linkages_arg(ignore_linkages)

  if (length(glycans) == 0) {
    return(rlang::set_names(logical(0L), character(0L)))
  }

  # Decide motif type and alignment
  motif_type <- get_motif_type(motif)
  if (motif_type == "known") {
    alignment <- decide_alignment(motif, alignment)
  } else if (is.null(alignment)) {
    alignment <- "substructure"
  }

  # Ensure glycans and motif are graphs
  glycans <- ensure_glycans_are_graphs(glycans)
  motif <- ensure_motif_is_graph(motif, motif_type)

  # Ensure NE graphs
  glycans <- ensure_ne_graphs(glycans)
  motif <- ensure_ne_graph(motif)

  # Ensure mono types are the same
  if (!same_mono_types(glycans)) {
    rlang::abort("All glycans must have the same monosaccharide type.")
  }
  glycans <- ensure_glycans_mono_type(glycans, motif)

  func(glycans, motif, alignment, ignore_linkages)
}


many_glycans_many_motifs_wrapper <- function(glycans, motifs, alignments, ignore_linkages, func) {
  # Check input arguments
  valid_glycans_arg(glycans)
  valid_motifs_arg(motifs)
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  # Get motif type and default motifs
  if (is.null(motifs)) {
    motifs <- available_motifs()
    motif_type <- "known"
  } else {
    motif_type <- get_motifs_type(motifs)
  }

  # Decide alignments
  if (motif_type == "known") {
    alignments <- decide_alignments(motifs, alignments)
  } else if (is.null(alignments)) {
    alignments <- "substructure"
  }

  # Ensure glycans and motifs are graphs
  glycans <- ensure_glycans_are_graphs(glycans)
  motifs <- ensure_motifs_are_graphs(motifs, motif_type)

  # Ensure NE graphs
  glycans <- ensure_ne_graphs(glycans)
  motifs <- ensure_ne_graphs(motifs)

  # Ensure mono types are the same
  if (motif_type != "known" && !same_mono_types(motifs)) {
    rlang::abort("All motifs must have the same monosaccharide type.")
  }
  if (!same_mono_types(glycans)) {
    rlang::abort("All glycans must have the same monosaccharide type.")
  }
  glycans <- ensure_glycans_mono_type(glycans, motifs[[1]])

  func(glycans, motifs, alignments, ignore_linkages)
}


# ----- Utility Functions -----
simple_has_motifs <- function(glycan, motifs, alignments, ignore_linkages) {
  purrr::map2_lgl(
    motifs, alignments,
    ~ has_motif_(glycan, .x, alignment = .y, ignore_linkages = ignore_linkages)
  )
}


simple_counts_motifs <- function(glycan, motifs, alignments, ignore_linkages) {
  purrr::map2_int(
    motifs, alignments,
    ~ counts_motif_(glycan, .x, alignment = .y, ignore_linkages = ignore_linkages)
  )
}


try_parse_iupac_condensed <- function(iupac) {
  purrr::possibly(glyparse::parse_iupac_condensed, NA)(iupac)
}
