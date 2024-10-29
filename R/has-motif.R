#' Check if a Glycan has the Given Motif
#'
#' @description
#' This function checks if a `glycan` has a given `motif`.
#' Technically speaking, it performs a subgraph isomorphism test to
#' determine if the `motif` is a subgraph of the `glycan`.
#' Both monosaccharides and linkages are considered in the comparison by default.
#' If `ignore_linkages` is set to `TRUE`, linkages will be ignored in the comparison.
#'
#' Both `glycan` and `motif` should be 'glycan_graph' objects
#' (see [glyrepr::as_glycan_graph()]).
#' They can be either "NE" or "DN" glycan graphs (can be different).
#'
#' Also, they can have different monosaccharide types
#' ("concrete", "generic" or "simple", see [glyrepr::decide_mono_type()]).
#' However, the monosaccharide type of `glycan` cannot be obscurer than that of `motif`,
#' which will raise an error.
#' For example, a "concrete" `glycan` can have a "generic" `motif`, but not vice versa.
#'
#' Obscure linkages (e.g. "??-?") are allowed in the `motif` graph
#' (see [glyrepr::possible_linkages()]).
#' "?" in a motif graph means "anything could be OK",
#' so it will match any linkage in the `glycan` graph.
#' However, "?" in a `glycan` graph will only match "?" in the `motif` graph.
#'
#' According to the [GlycoMotif](https://glycomotif.glyomics.org) database,
#' a motif can be classified into four alignment types:
#' - "substructure": The motif can be anywhere in the glycan. This is the default.
#' See [substructure](https://glycomotif.glyomics.org/glycomotif/Substructure_Alignment)
#' for details.
#' - "core": The motif must align with at least one connected substructure
#' (subtree) at the reducing end of the glycan.
#' See [glycan core](https://glycomotif.glyomics.org/glycomotif/Glycan_Core_Alignment)
#' for details.
#' - "terminal": The motif must align with at least one connected substructure
#' (subtree) at the nonreducing end of the glycan.
#' See [nonreducing end](https://glycomotif.glyomics.org/glycomotif/Nonreducing-End_Alignment)
#' for details.
#' - "whole": The motif must align with the entire glycan.
#' See [whole-glycan](https://glycomotif.glyomics.org/glycomotif/Whole-Glycan_Alignment)
#' for details.
#'
#' Please see the Examples section if you are confused.
#' And also see the documentation of key functions listed above.
#'
#' @details
#' Under the hood, if `alignment` is "whole", the function uses
#' [igraph::isomorphic()] with "vf2" method to perform the isomorphism test.
#' Otherwise, it uses [igraph::graph.get.subisomorphisms.vf2()] to get all possible
#' subgraph isomorphisms between the `glycan` and `motif` graphs.
#' Vextex attributes and edge attributes of the `glycan` and `motif` graphs are colorized
#' to add "color" attributes to the vertices and edges of the graphs.
#'
#' To allow obscure linkage matching, [glyrepr::possible_linkages()] is used to
#' generate all possible versions of the `motif` graph.
#' Then, the subgraph isomorphism test is performed for each version,
#' returning `TRUE` if any of them is isomorphic to the `glycan` graph.
#' This implementation could suffer from performance issues when the `motif` graph
#' has many obscure linkages.
#' However, it is the most straightforward way to handle obscure linkages.
#' Future improvements may fine-tune the "vf2" method or use other methods
#' to support wildcard matching directly.
#'
#' @param glycan A 'glycan_graph' object.
#' @param motif A 'glycan_graph' object.
#' @param ... Not used.
#' @param alignment A character string. Possible values are "substructure", "core", "terminal" and "whole".
#' See description for details. Default is "substructure".
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical value indicating if the `glycan` has the `motif`.
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' (glycan <- glyrepr::o_glycan_core_2(mode = "ne", mono_type = "concrete"))
#'
#' # The glycan has the motif "Gal(b1-3)GalNAc"
#' motif_1 <- parse_iupac_condensed("Gal(b1-3)GalNAc")
#' has_motif(glycan, motif_1)
#'
#' # But not "Gal(b1-4)GalNAc" (wrong linkage)
#' motif_2 <- parse_iupac_condensed("Gal(b1-4)GalNAc")
#' has_motif(glycan, motif_2)
#'
#' # Set `ignore_linkages` to `TRUE` to ignore linkages
#' has_motif(glycan, motif_2, ignore_linkages = TRUE)
#'
#' # Glycan and motif can have different graph modes
#' motif_1_dn <- convert_graph_mode(motif_1, to = "dn")
#' has_motif(glycan, motif_1_dn)
#'
#' # And different monosaccharide types
#' motif_1_generic <- convert_glycan_mono_type(motif_1, "generic")
#' has_motif(glycan, motif_1_generic)
#'
#' # However, the monosaccharide type of `glycan` cannot be obscurer than that of `motif`
#' glycan_simple <- convert_glycan_mono_type(glycan, "simple")
#' try(has_motif(glycan_simple, motif_1))
#'
#' # Obscure linkages in the `motif` graph are allowed
#' motif_3 <- parse_iupac_condensed("Gal(b1-?)GalNAc")
#' has_motif(glycan, motif_3)
#'
#' # However, obscure linkages in `glycan` will only match "?" in the `motif` graph
#' glycan_2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' has_motif(glycan_2, motif_1)
#' has_motif(glycan_2, motif_3)
#'
#' # Alignment types
#' glycan_3 <- parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_4 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_5 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal")
#' motif_6 <-  parse_iupac_condensed(         "Gal(a1-4)Gal(a1-6)Gal")
#' motif_7 <-  parse_iupac_condensed(         "Gal(a1-4)Gal")
#' motifs <- list(motif_4, motif_5, motif_6, motif_7)
#'
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "whole"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "core"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "terminal"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "substructure"))
#'
#' @export
has_motif <- function(glycan, motif, ..., alignment = "substructure", ignore_linkages = FALSE) {
  # Check input arguments
  if (!glyrepr::is_glycan(glycan) || !glyrepr::is_glycan(motif)) {
    rlang::abort("`glycan` and `motif` must be 'glycan_graph' objects.")
  }

  # Ensure that `glycan` and `motif` have the same mode
  motif <- ensure_motif_graph_mode(glycan, motif)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycan_mono_type()`.
  glycan <- ensure_glycan_mono_type(glycan, motif)

  # Check if the motif is a subgraph of the glycan
  if (ignore_linkages) {
    glycan <- glyrepr::remove_linkages(glycan)
    motif <- glyrepr::remove_linkages(motif)
    vf2_subgraph_isomorphic(glycan, motif, alignment)
  } else {
    # Impute the "?" in the linkage of the motif graph with every possibility
    motif_candidates <- impute_linkages(motif)
    any(purrr::map_lgl(
      motif_candidates,
      ~ vf2_subgraph_isomorphic(glycan, .x, alignment)
    ))
  }
}


ensure_motif_graph_mode <- function(glycan, motif) {
  # Make motif the same graph mode (NE or DN) as the glycan
  if (glyrepr::is_ne_glycan(glycan)) {
    glyrepr::convert_graph_mode(motif, to = "ne", strict = FALSE)
  } else {
    glyrepr::convert_graph_mode(motif, to = "dn", strict = FALSE)
  }
}


ensure_glycan_mono_type <- function(glycan, motif) {
  # Make the glycan the same monosaccharide type (concrete, generic or simple) as the motif.
  # If the glycan type is lower than the motif type, an error will be raised.
  glycan_type <- glyrepr::decide_glycan_mono_type(glycan)
  motif_type <- glyrepr::decide_glycan_mono_type(motif)

  conditions <- (
    (glycan_type == "concrete" && motif_type != "concrete") ||
    (glycan_type == "generic" && motif_type == "simple")
  )
  if (conditions) return(glyrepr::convert_glycan_mono_type(glycan, motif_type))

  if (glycan_type != motif_type) {
    cli::cli_abort(c(
      "The monosaccharide type of `glycan` cannot be obscurer than `motif`.",
      "x" = "{.val {glycan_type}} is obscurer than {.val {motif_type}}."
    ))
  }
  return(glycan)
}


colorize_glycan_graphs <- function(glycan, motif) {
  # Add "color" attributes to the vertices and edges of a glycan graph.
  # The "color" attrs wiil be usd by `igraph::subgraph_isomorphic()`
  # to compare glycan graphs (with method = "vf2").
  #
  # This function assumes that `glycan` and `input` are both "NE" or "DN" glycan graphs.
  if (inherits(glycan, "ne_glycan_graph")) {
    colorize_ne_glycan_graphs(glycan, motif)
  } else {
    colorize_dn_glycan_graphs(glycan, motif)
  }
}


colorize_ne_glycan_graphs <- function(glycan, motif) {
  node_colors <- colorize_labels(igraph::V(glycan)$mono, igraph::V(motif)$mono)
  igraph::V(glycan)$color <- node_colors[["glycan"]]
  igraph::V(motif)$color <- node_colors[["motif"]]

  edge_colors <- colorize_labels(igraph::E(glycan)$linkage, igraph::E(motif)$linkage)
  igraph::E(glycan)$color <- edge_colors[["glycan"]]
  igraph::E(motif)$color <- edge_colors[["motif"]]

  list(glycan = glycan, motif = motif)
}


colorize_dn_glycan_graphs <- function(glycan, motif) {
  get_labels <- function(graph) {
    dplyr::if_else(
      igraph::V(graph)$type == "mono",
      igraph::V(graph)$mono,
      igraph::V(graph)$linkage
    )
  }
  glycan_labels <- get_labels(glycan)
  motif_labels <- get_labels(motif)
  colors <- colorize_labels(glycan_labels, motif_labels)
  igraph::V(glycan)$color <- colors[["glycan"]]
  igraph::V(motif)$color <- colors[["motif"]]

  list(glycan = glycan, motif = motif)
}


colorize_labels <- function(glycan_labels, motif_labels) {
  # Helper function to colorize labels.
  unique_values <- unique(c(glycan_labels, motif_labels))
  color_map <- seq_along(unique_values)
  names(color_map) <- unique_values

  glycan_colors <- color_map[glycan_labels]
  names(glycan_colors) <- NULL
  motif_colors <- color_map[motif_labels]
  names(motif_colors) <- NULL

  list(glycan = glycan_colors, motif = motif_colors)
}


impute_linkages <- function(motif) {
  # Impute the "?" in the linkage of the motif graph with every possibility.
  # Return a list of motif graphs with different linkage imputations.
  if (glyrepr::is_ne_glycan(motif)) {
    impute_ne_linkages(motif)
  } else {
    impute_dn_linkages(motif)
  }
}


impute_ne_linkages <- function(motif) {
  linkages <- igraph::E(motif)$linkage
  candidates <- expand.grid(
    purrr::map(linkages, glyrepr::possible_linkages, include_unknown = TRUE),
    stringsAsFactors = FALSE
  )
  candidates <- purrr::list_transpose(unclass(candidates))
  candidates <- purrr::map(candidates, ~ unname(.x))
  purrr::map(candidates, ~ igraph::set_edge_attr(motif, "linkage", value = .x))
}


impute_dn_linkages <- function(motif) {
  ne_motif <- glyrepr::convert_graph_mode(motif, to = "ne")
  ne_result <- impute_ne_linkages(ne_motif)
  purrr::map(ne_result, glyrepr::convert_graph_mode, to = "dn")
}


vf2_subgraph_isomorphic <- function(glycan, motif, alignment) {
  # Colorize the glycan and motif graphs.
  # This is to add "color" attributes to the vertices and edges of a glycan graph,
  # which is demanded by the "vf2" method of `igraph::subgraph_isomorphic()`.
  colored_graphs <- colorize_glycan_graphs(glycan, motif)
  c_glycan <- colored_graphs[["glycan"]]
  c_motif <- colored_graphs[["motif"]]

  # Perform the subgraph isomorphism test with the "vf2" method.
  # Possible `alignment`: "substructure", "core", "terminal", "whole"
  if (alignment == "whole") {
    return(igraph::isomorphic(c_glycan, c_motif, method = "vf2"))
  }

  res <- igraph::graph.get.subisomorphisms.vf2(c_glycan, c_motif)
  if (alignment == "substructure") {
    return(length(res) > 0)
  }
  if (alignment == "core") {
    glycan_core_node <- core_node(c_glycan)
    return(any(purrr::map_lgl(res, ~ glycan_core_node %in% .x)))
  }
  if (alignment == "terminal") {
    glycan_terminal_nodes <- terminal_nodes(c_glycan)
    motif_terminal_nodes <- terminal_nodes(c_motif)
    return(any(purrr::map_lgl(
      res, ~ all(.x[motif_terminal_nodes] %in% glycan_terminal_nodes)))
    )
  }
}


terminal_nodes <- function(glycan) {
  out_degree <- igraph::degree(glycan, mode = "out")
  igraph::V(glycan)[out_degree == 0]
}


core_node <- function(glycan) {
  in_degree <- igraph::degree(glycan, mode = "in")
  igraph::V(glycan)[in_degree == 0]
}
