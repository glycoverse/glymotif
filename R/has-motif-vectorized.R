#' Check if Glycan(s) have the Given Motif(s)
#'
#' These are vectorized versions of [has_motif()], vectorizing over
#' motifs (`has_motifs`) and glycans (`have_motif`), respectively.
#' They can be told apart by English grammar :).
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
#' @export
has_motifs <- function(glycan, motifs, ...) {
  res <- purrr::map_lgl(motifs, ~ try_has_motif(glycan, .x, ...))
  if (any(is.na(res))) {
    cli::cli_abort("Motifs at indices {.val {which(is.na(res))}} are neither
                   known motif names or able to be parsed as IUPAC-condensed
                   structure strings.")
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
  res
}


try_has_motif <- function(glycan, motif, ...) {
  tryCatch(has_motif(glycan, motif, ...), error = function(e) NA)
}
