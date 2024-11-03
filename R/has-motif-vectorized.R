#' Check if Glycan(s) have the Given Motif(s)
#'
#' @description
#' These are vectorized versions of [has_motif()], vectorizing over
#' motifs (`has_motifs()`) and glycans (`have_motif()`), respectively.
#' They can be told apart by English grammar :).
#'
#' Names of glycans or motifs
#' are preserved as names of the output logical vectors.
#' If no names are provided, and the input is a character vector, the names
#' will be the vector itself.
#' It is always recommended to provide names for the input vectors.
#'
#' If `motifs` are missing in `has_motifs()`, all known motifs from
#' GlycoMotif GlyGen Collection will be used.
#'
#' @param glycan A 'glycan_graph' object, or an IUPAC-condensed structure string.
#' @param motif A 'glycan_graph' object, an IUPAC-condensed structure string,
#' or a known motif name (use [available_motifs()] to see all available motifs).
#' @param glycans A list of 'glycan_graph' objects, or a character vector of
#' IUPAC-condensed structure strings.
#' @param motifs A list of 'glycan_graph' objects, a character vector of
#' IUPAC-condensed structure strings, or a character vector of known motif names.
#' @param ... Additional arguments passed to [has_motif()].
#'
#' @return A logical vector, indicating if the glycan(s) have the motif(s).
#'
#' @examples
#' library(glyparse)
#'
#' # Check if a glycan has the given motifs using names
#' glycan <- "Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' has_motifs(glycan, c("O-Glycan core 1", "O-Glycan core 2"))
#'
#' # using IUPAC-condensed structure strings
#' has_motifs(glycan, c("Gal(b1-3)GalNAc(a1-", "GlcNAc(b1-4)GlcNAc(a1-"))
#'
#' # or using 'glycan_graph' objects
#' motifs <- list(
#'   M1 = parse_iupac_condensed("Gal(b1-3)GalNAc(a1-"),
#'   M2 = parse_iupac_condensed("GlcNAc(b1-4)GlcNAc(a1-")
#' )
#' has_motifs(glycan, motifs)
#'
#' # If motifs are missing, use all known motifs
#' has_motifs(glycan)[1:5]
#'
#' # Check if glycans have the given motif
#' glycans <- c(
#'   G1 = "Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)Gal",
#'   G2 = "Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' )
#' have_motif(glycans, "Gal(b1-3)Gal")
#'
#' @export
has_motifs <- function(glycan, motifs = NULL, ...) {
  if (missing(motifs)) {
    motifs <- available_motifs()
  }
  res <- purrr::map_lgl(motifs, ~ try_has_motif(glycan, .x, ...))
  if (any(is.na(res))) {
    cli::cli_abort("Motifs at indices {.val {which(is.na(res))}} are neither
                   known motif names or able to be parsed as IUPAC-condensed
                   structure strings.")
  }
  if (is.null(names(res)) && is.character(motifs)) {
    names(res) <- motifs
  }
  res
}


#' @rdname has_motifs
#' @export
have_motif <- function(glycans, motif, ...) {
  res <- purrr::map_lgl(glycans, try_has_motif, motif, ...)
  if (any(is.na(res))) {
    cli::cli_abort("Glycans at indices {.val {which(is.na(res))}} are not able
                   to be parsed as IUPAC-condensed structure strings.")
  }
  if (is.null(names(res)) && is.character(glycans)) {
    names(res) <- glycans
  }
  res
}


try_has_motif <- function(glycan, motif, ...) {
  tryCatch(has_motif(glycan, motif, ...), error = function(e) NA)
}
