#' Low-Level Motif Matching on Graphs
#'
#' @description
#' These functions are low-level variants of [have_motif()], [count_motif()],
#' and [match_motif()] for package code that already has compatible igraph
#' objects from [glyrepr::get_structure_graphs()].
#'
#' @param glycan_graph An igraph glycan graph.
#' @param motif_graph An igraph motif graph.
#' @param alignment A character scalar: `"substructure"`, `"core"`,
#'   `"terminal"`, or `"whole"`.
#' @param ignore_linkages A logical scalar. If `TRUE`, linkages are ignored.
#' @param strict_sub A logical scalar. If `TRUE`, substituents are matched
#'   strictly.
#' @param match_degree A logical vector indicating which motif nodes must match
#'   the glycan's in- and out-degree exactly. A scalar is recycled to the
#'   number of motif nodes.
#'
#' @details
#' These functions do no validation, parsing, naming, or monosaccharide-type
#' conversion. Callers must provide valid, already-compatible graph objects.
#' In particular, when matching a generic motif graph, callers must pass a
#' glycan graph whose `mono` vertex attributes have already been converted to
#' the compatible generic representation.
#'
#' These functions never call [glyrepr::as_glycan_structure()].
#'
#' @returns
#' - `.g_have_motif()` returns a logical scalar.
#' - `.g_count_motif()` returns an integer scalar.
#' - `.g_match_motif()` returns a list of integer vectors.
#'
#' @seealso [have_motif()], [count_motif()], [match_motif()]
#'
#' @name dot-g_motif
NULL

#' @rdname dot-g_motif
#' @export
.g_have_motif <- function(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  apply_single_motif_to_graph(
    glycan_graph = glycan_graph,
    motif_graph = motif_graph,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    single_glycan_func = .have_motif_single
  )
}

#' @rdname dot-g_motif
#' @export
.g_count_motif <- function(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  apply_single_motif_to_graph(
    glycan_graph = glycan_graph,
    motif_graph = motif_graph,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    single_glycan_func = .count_motif_single
  )
}

#' @rdname dot-g_motif
#' @export
.g_match_motif <- function(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  apply_single_motif_to_graph(
    glycan_graph = glycan_graph,
    motif_graph = motif_graph,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    single_glycan_func = .match_motif_single
  )
}

#' Apply a Single Motif Graph to a Single Glycan Graph
#'
#' @param glycan_graph An igraph glycan graph.
#' @param motif_graph An igraph motif graph.
#' @param alignment A character scalar.
#' @param ignore_linkages A logical scalar.
#' @param strict_sub A logical scalar.
#' @param match_degree A logical vector or `NULL`.
#' @param single_glycan_func A graph-level motif matching function.
#'
#' @return The result from `single_glycan_func`.
#' @noRd
apply_single_motif_to_graph <- function(
  glycan_graph,
  motif_graph,
  alignment,
  ignore_linkages,
  strict_sub,
  match_degree,
  single_glycan_func
) {
  motif_has_linkages <- graph_has_linkages(motif_graph)
  motif_composition_profile <- new_motif_composition_profile(motif_graph)
  match_degree <- normalize_graph_match_degree(match_degree, motif_graph)

  single_glycan_func(
    glycan_graph = glycan_graph,
    motif_graph = motif_graph,
    motif_has_linkages = motif_has_linkages,
    motif_composition_profile = motif_composition_profile,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree
  )
}

#' Normalize Graph-Level Match Degree
#'
#' @param match_degree A logical vector or `NULL`.
#' @param motif_graph An igraph motif graph.
#'
#' @return A normalized logical vector or `NULL`.
#' @noRd
normalize_graph_match_degree <- function(match_degree, motif_graph) {
  if (is.null(match_degree) || length(match_degree) != 1) {
    return(match_degree)
  }

  rep(match_degree, igraph::vcount(motif_graph))
}
