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
#' # About Names
#'
#' `match_motif()` perserve names from the input `glycans` vector.`
#'
#' For `match_motifs()`, the outermost list is named by motifs,
#' and the inner lists are named by glycans,
#' following the same rules as in `have_motifs()` and `count_motifs()`.
#'
#' # Floating parts and substituents
#'
#' Glycans with unresolved floating parts or substituents are matched across
#' every conflict-free localization allowed by their candidate-parent domains.
#' These functions return the union of node mappings from every localization,
#' using node indices from the original unresolved structure.
#'
#' Motifs must be connected structures and therefore cannot themselves contain
#' unresolved floating parts or substituents.
#'
#' Matching supports up to 256 raw candidate-parent combinations per glycan.
#' Localize floating parts and substituents with
#' [glyrepr::localize_floating_parts()] first for larger domains.
#'
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
#'   The outermost list is named by `motifs` if they have names.
#'   The middle list is named by `glycans` if they have names.
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
#' igraph::V(glycan_graph)$mono # 1, 2, 3, 4, 5
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
match_motif <- function(
  glycans,
  motif,
  ...,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
) {
  rlang::check_dots_empty()

  # Store input names before processing
  glycan_names <- names(glycans)

  .assert_glycan_structure(glycans, "glycans")
  .assert_glycan_structure(motif, "motif")
  params <- prepare_motif_args(
    glycans = glycans,
    motifs = motif,
    alignments = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree,
    single_motif = TRUE,
    strict_sub = strict_sub,
    mode = mode
  )
  result <- rlang::exec("match_motif_", !!!params)

  # Apply names to result if input had names
  if (!is.null(glycan_names)) {
    names(result) <- glycan_names
  }

  result
}

#' @rdname match_motif
#' @export
match_motifs <- function(
  glycans,
  motifs,
  ...,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
) {
  rlang::check_dots_empty()

  # Validate glycans first (must be glycan_structure)
  .assert_glycan_structure(glycans, "glycans")

  # Store input names before processing
  glycan_names <- names(glycans)
  is_branch_spec <- inherits(motifs, "branch_motifs_spec")

  params <- prepare_motif_args(
    glycans = glycans,
    motifs = motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree,
    single_motif = FALSE,
    strict_sub = strict_sub,
    mode = mode
  )

  # Validate resolved motifs (must be glycan_structure, not raw strings)
  .assert_glycan_structure(params$motifs, "motifs")

  # Use names from resolved motifs if available (e.g., from dynamic_motifs/branch_motifs)
  # Otherwise fall back to prepare_motif_names on the original input
  if (!is.null(names(params$motifs))) {
    motif_names <- names(params$motifs)
  } else {
    motif_names <- prepare_motif_names(motifs)
  }

  results <- rlang::exec(
    "match_motifs_",
    !!!params,
    glycan_names = glycan_names,
    motif_names = motif_names
  )

  if (is_branch_spec && is.list(params$match_degree)) {
    results <- .trim_result_for_branch_motifs(results)
  }

  results
}

.assert_glycan_structure <- function(x, arg_name) {
  if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "Argument {.arg {arg_name}} must be a 'glyrepr_structure' object.",
      "i" = "Use `glyparse::auto_parse()` to parse glycan structure strings into `glyrepr_structure`."
    ))
  }
}

.trim_result_for_branch_motifs <- function(res) {
  purrr::map(res, function(motif_res) {
    purrr::map(motif_res, function(glycan_res) {
      purrr::map(glycan_res, ~ .x[1:(length(.x) - 4)])
    })
  })
}

#' Internal verison of `match_motif()`
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` object.
#' @param motif A `glyrepr_structure` object with length 1.
#' @param alignment A character scalar.
#'   Possible values are "substructure", "core", "terminal", and "whole".
#' @param ignore_linkages A logical value.
#'
#' @noRd
match_motif_ <- function(
  glycans,
  motif,
  alignment,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  strict_floating = TRUE
) {
  apply_single_motif_to_glycans(
    glycans = glycans,
    motif = motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = mode,
    strict_floating = strict_floating,
    single_glycan_func = .match_motif_single,
    smap_func = glyrepr::smap,
    result_type = "list"
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
match_motifs_ <- function(
  glycans,
  motifs,
  alignments,
  glycan_names,
  motif_names,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  strict_floating = TRUE
) {
  match_degree_list <- if (is.null(match_degree)) {
    rep(list(NULL), length(motifs))
  } else {
    match_degree
  }

  apply_motifs_to_glycans(
    glycans = glycans,
    motifs = motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    single_glycan_func = .match_motif_single,
    glycan_names = glycan_names,
    motif_names = motif_names,
    strict_sub = strict_sub,
    match_degree = match_degree_list,
    mode = mode,
    strict_floating = strict_floating,
    result_type = "list"
  )
}

.match_motif_single <- function(
  glycan_graph,
  motif_graph,
  motif_has_linkages,
  motif_composition_profile,
  alignment,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL,
  compatibility_cache = NULL
) {
  if (
    !whole_alignment_size_can_match(
      glycan_graph,
      motif_graph,
      alignment,
      glycan_batch_profile = glycan_batch_profile,
      motif_batch_profile = motif_batch_profile
    )
  ) {
    return(list())
  }
  if (
    !core_alignment_root_can_match(
      glycan_graph,
      motif_graph,
      alignment,
      strict_sub = strict_sub,
      mode = mode,
      glycan_batch_profile = glycan_batch_profile,
      motif_batch_profile = motif_batch_profile
    )
  ) {
    return(list())
  }

  linkage_match_mode <- resolve_linkage_match_mode(
    glycan_graph,
    motif_has_linkages,
    ignore_linkages,
    mode = mode,
    glycan_batch_profile = glycan_batch_profile
  )

  # "none" is an early no-match result: the motif requires linkage information,
  # but the glycan has none, so there are no mappings to return.
  if (linkage_match_mode == "none") {
    return(list())
  }

  if (
    !composition_can_match(
      glycan_graph,
      motif_composition_profile,
      glycan_batch_profile = glycan_batch_profile,
      motif_batch_profile = motif_batch_profile
    )
  ) {
    return(list())
  }

  perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = alignment,
    ignore_linkages = linkage_match_mode == "ignore",
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = mode,
    glycan_batch_profile = glycan_batch_profile,
    motif_batch_profile = motif_batch_profile,
    compatibility_cache = compatibility_cache
  )
}
