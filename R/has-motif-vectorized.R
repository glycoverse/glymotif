#' Check if the Glycan(s) have the Given Motif(s)
#'
#' @description
#' These are vectorized versions of [has_motif()].
#' - `has_motif()` checks if a motif is present in a glycan.
#' - `has_motifs()` checks if a list of motifs are present in a glycans.
#' - `have_motif()` checks if a motif is present in a list of glycans.
#' - `have_motifs()` checks if a list of motifs are present in a list of glycans.
#'
#' They can be told apart by English grammar:
#' - Functions starting with "has" work with one glycan, "have" with multiple glycans.
#' - Functions ending with "motif" work with one motif, "motifs" with multiple motifs.
#'
#' @details
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
#' Internally, all glycan and motif inputs are converted to 'glycan_graph' objects.
#' Therefore, if performance matters and you need to repeatedly use the same glycans or motifs,
#' it is recommended to convert them to 'glycan_graph' objects beforehand.
#' [get_motif_graph()] can be used to get the 'glycan_graph' object of a known motif.
#' [glyparse::parse_iupac_condensed()] can be used to convert IUPAC-condensed
#' structure strings to 'glycan_graph' objects.
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
#' @return #' `has_motifs()` and `have_motif()` return a logical vector
#' of the same length as the input motifs or glycans, respectively.
#' `have_motifs()` returns a matrix of logical values,
#' indicating if each glycan has each motif.
#' Rows are glycans and columns are motifs.
#'
#' The dimension names are determined by the input arguments.
#' If `glycans` or `motifs` have names, they will be preserved in the result.
#' If not, but `glycans` or `motifs` are character vectors, they will be used as names.
#'
#' @seealso [has_motif()]
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

  simple_has_motifs(glycan, motifs, alignments, ignore_linkages)
}


#' @rdname has_motifs
#' @export
have_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
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

  purrr::map_lgl(glycans, has_motif_, motif, alignment = alignment, ignore_linkages = ignore_linkages)
}


#' @rdname has_motifs
#' @export
have_motifs <- function(glycans, motifs = NULL, alignments = NULL, ignore_linkages = FALSE, simplify = FALSE) {
  # Check input arguments
  valid_glycans_arg(glycans)
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

  lgl_list <- purrr::map(glycans, simple_has_motifs, motifs, alignments = alignments, ignore_linkages = ignore_linkages)
  result <- do.call(rbind, lgl_list)

  if (simplify) {
    # Drop columns will all FALSE
    result <- result[, colSums(result) > 0]
  }

  result
}


simple_has_motifs <- function(glycan, motifs, alignments, ignore_linkages) {
  purrr::map2_lgl(
    motifs, alignments,
    ~ has_motif_(glycan, .x, alignment = .y, ignore_linkages = ignore_linkages)
  )
}


try_parse_iupac_condensed <- function(iupac) {
  purrr::possibly(glyparse::parse_iupac_condensed, NA)(iupac)
}
