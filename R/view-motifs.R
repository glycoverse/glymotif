#' View motif matches on a glycan
#'
#' Visualize where a motif matches a glycan structure.
#'
#' `view_motif()` matches one `motif` against one `glycan` with the same
#' matching rules used by [match_motif()], then draws the glycan with the
#' matched residues highlighted.
#'
#' @param glycan One of:
#'   - A scalar [glyrepr::glycan_structure()].
#'   - A glycan structure string. All formats supported by [glyparse::auto_parse()]
#'     are accepted, including IUPAC-condensed, WURCS, GlycoCT, and others.
#' @inheritParams have_motif
#'
#' @returns
#' A `ggplot` object returned by [glydraw::draw_cartoon()]. If no match is found,
#' the glycan is drawn without highlighted residues and a cli alert is emitted.
#'
#' @seealso [match_motif()], [glydraw::draw_cartoon()]
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' glycan <- n_glycan_core()
#' motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
#'
#' \dontrun{
#' view_motif(glycan, motif)
#' }
#' @export
view_motif <- function(
  glycan,
  motif,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  params <- prepare_motif_args(
    glycans = glycan,
    motifs = motif,
    alignments = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree,
    strict_sub = strict_sub
  )
  if (length(glycan) != 1) {
    cli::cli_abort(c(
      "Only one glycan can be visualized at a time.",
      "x" = "Got {.val {length(glycan)}} glycans."
    ))
  }
  if (length(motif) != 1) {
    cli::cli_abort(c(
      "Only one motif can be visualized at a time.",
      "x" = "Got {.val {length(motif)}} motifs."
    ))
  }

  match_res <- match_motif(
    params$glycans,
    params$motifs,
    params$alignments,
    ignore_linkages = params$ignore_linkages,
    strict_sub = params$strict_sub,
    match_degree = params$match_degree
  )
  match_res <- match_res[[1]]
  idx <- unique(unlist(match_res))

  if (length(idx) == 0) {
    cli::cli_alert_danger("No matches found between the glycan and motif.")
    glydraw::draw_cartoon(params$glycans, highlight = integer(0))
  } else {
    glydraw::draw_cartoon(params$glycans, highlight = idx)
  }
}
