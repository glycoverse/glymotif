#' Match Motif(s) in Glycans
#'
#' @description
#' These functions find all occurrences of the given `motif`(s) in the `glycans`.
#' Node-to-node mapping is returned for each match.
#' This function is NOT useful for most users if you are not interested in the concrete node mapping.
#' See [have_motif()] and [count_motif()] for more information about the matching rules.
#'
#' - `match_motif()` matches a single motif against multiple glycans
#' - `match_motifs()` matches multiple motifs against multiple glycans
#'
#' Different from [have_motif()] and [count_motif()],
#' these functions return detailed match information.
#' More specifically, for each glycan-motif pair,
#' a integer vector is returned,
#' indicating the node mapping from the motif to the glycan.
#' For example, if the vector is `c(2, 3, 6)`,
#' it means that the first node in the motif matches the 2nd node in the glycan,
#' the second node in the motif matches the 3rd node in the glycan,
#' and the third node in the motif matches the 6th node in the glycan.
#'
#' Node indices are only meaningful for [glyrepr::glycan_structure()],
#' so only [glyrepr::glycan_structure()] is supported for `glycans` and `motifs`.
#'
#' @details
#' # Vertex and Linkage Indices
#'
#' The indices of vertices and linkages in a glycan correspond directly to their
#' order in the IUPAC-condensed string, which is printed when you print a
#' [glyrepr::glycan_structure()].
#' For example, for the glycan `Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-)`,
#' the vertices are "Man", "Man", "Man", "GlcNAc", "GlcNAc",
#' and the linkages are "a1-3", "a1-6", "b1-4", "b1-4".
#'
#' Thus, matching the motif "Man(a1-3)Man(b1-4)" to this glycan yields `c(1, 3)`.
#' This indicates that the first motif vertex (the a1-3 Man) corresponds to
#' the first vertex in the glycan, and the second motif vertex (the b1-4 Man)
#' corresponds to the third vertex in the glycan.
#'
#' @param glycans A `glyrepr_structure` object.
#' @param motif A `glyrepr_structure` object with length 1.
#' @param motifs A `glyrepr_structure` object.
#' @inheritParams have_motif
#'
#' @returns
#' A nested list of integer vectors.
#' - `match_motif()`: Two levels of nesting.
#'   The outer list corresponds to glycans, and the inner list corresponds to matches.
#'   Use `purrr::pluck(result, glycan_index, match_index)` to access the match information.
#'   For example, `purrr::pluck(result, 1, 2)` means the 2nd match in the 1st glycan.
#' - `match_motifs()`: Three levels of nesting.
#'   The outermost list corresponds to motifs, the middle list corresponds to glycans,
#'   and the innermost list corresponds to matches.
#'   Use `purrr::pluck(result, motif_index, glycan_index, match_index)` to access the match information.
#'   For example, `purrr::pluck(result, 1, 2, 3)` means the 3rd match in the 2nd glycan for the 1st motif.
#'
#' @seealso [have_motif()], [count_motif()]
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' (glycan <- n_glycan_core())
#'
#' # Let's peek under the hood of the nodes in the glycan
#' glycan_graph <- get_structure_graphs(glycan)
#' igraph::V(glycan_graph)$mono  # 1, 2, 3, 4, 5
#'
#' # Match a single motif against a single glycan
#' motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
#' match_motif(glycan, motif)
#'
#' # Match multiple motifs against a single glycan
#' motifs <- c(
#'   "Man(a1-3)[Man(a1-6)]Man(b1-",
#'   "Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
#' )
#' motifs <- parse_iupac_condensed(motifs)
#' match_motifs(glycan, motifs)
#'
#' @export
match_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  .assert_glycan_structure(glycans, "glycans")
  .assert_glycan_structure(motif, "motif")
  params <- prepare_have_motif_args(glycans, motif, alignment, ignore_linkages)
  rlang::exec("match_motif_", !!!params)
}

#' @rdname match_motif
#' @export
match_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  .assert_glycan_structure(glycans, "glycans")
  .assert_glycan_structure(motifs, "motifs")
  params <- prepare_have_motifs_args(glycans, motifs, alignments, ignore_linkages)
  rlang::exec("match_motifs_", !!!params)
}

.assert_glycan_structure <- function(x, arg_name) {
  if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "Argument {.arg {arg_name}} must be a 'glyrepr_structure' object.",
      "i" = "Use `glyparse::auto_parse()` to parse glycan structure strings into `glyrepr_structure`."
    ))
  }
}

#' Internal verison of `match_motif()`
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` object.
#' @param motif A `glyrepr_structure` object with length 1.
#' @param alignment A character scalar.
#' @param ignore_linkages A logical value.
#'
#' @noRd
match_motif_ <- function(glycans, motif, alignment, ignore_linkages = FALSE) {
  apply_single_motif_to_glycans(
    glycans = glycans,
    motif = motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    single_glycan_func = .match_motif_single,
    smap_func = glyrepr::smap
  )
}

#' Internal verison of `match_motifs()`
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` object.
#' @param motifs A `glyrepr_structure` object.
#' @param alignments A character vector with the same length as `motifs`.
#' @param ignore_linkages A logical value.
#'
#' @noRd
match_motifs_ <- function(glycans, motifs, alignments, glycan_names, motif_names, ignore_linkages = FALSE) {
  purrr::map2(motifs, alignments, ~ match_motif_(glycans, .x, alignment = .y, ignore_linkages = ignore_linkages))
}

.match_motif_single <- function(glycan_graph, motif_graph, alignment, ignore_linkages = FALSE) {
  c_graphs <- colorize_graphs(glycan_graph, motif_graph)
  glycan_graph <- c_graphs$glycan
  motif_graph <- c_graphs$motif
  res <- perform_vf2(glycan_graph, motif_graph)
  valid_mask <- purrr::map_lgl(
    res, is_valid_result, glycan = glycan_graph, motif = motif_graph,
    alignment = alignment, ignore_linkages = ignore_linkages
  )
  valid_res <- res[valid_mask]
  unique_res <- unique_vf2_res(valid_res)
  purrr::map(unique_res, as.integer)
}
