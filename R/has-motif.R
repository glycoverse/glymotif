#' Check if a Glycan has the Given Motif
#'
#' @description
#' This function checks if a `glycan` has a given `motif`.
#' Technically speaking, it performs a subgraph isomorphism test to
#' determine if the `motif` is a subgraph of the `glycan`.
#' Both monosaccharides and linkages are considered in the comparison by default.
#'
#' @details
#' # Graph mode and monosaccharide type
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
#' # Linkages
#'
#' Obscure linkages (e.g. "??-?") are allowed in the `motif` graph
#' (see [glyrepr::possible_linkages()]).
#' "?" in a motif graph means "anything could be OK",
#' so it will match any linkage in the `glycan` graph.
#' However, "?" in a `glycan` graph will only match "?" in the `motif` graph.
#' You can set `ignore_linkages = TRUE` to ignore linkages in the comparison.
#'
#' # Alignment
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
#' # Substituents
#'
#' Substituents (e.g. "Ac", "SO3") are matched in strict mode.
#' "Neu5Ac-9Ac" will only match "Neu5Ac-9Ac" but not "Neu5Ac",
#' and "Neu5Ac" will not match "Neu5Ac-9Ac".
#'
#' # Implementation
#'
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
#' # The anomer of the motif will be matched to linkages in the glycan
#' motif_4 <- parse_iupac_condensed("GlcNAc(b1-")
#' has_motif(glycan_2, motif_4)
#'
#' # Alignment types
#' # The default type is "substructure", which means the motif can be anywhere in the glycan.
#' # Other options include "core", "terminal" and "whole".
#' glycan_3 <- parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_5 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_6 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal")
#' motif_7 <-  parse_iupac_condensed(         "Gal(a1-4)Gal(a1-6)Gal")
#' motif_8 <-  parse_iupac_condensed(         "Gal(a1-4)Gal")
#' motifs <- list(motif_5, motif_6, motif_7, motif_8)
#'
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "whole"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "core"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "terminal"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "substructure"))
#'
#' # Substituents
#' glycan_4 <- parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc")
#' glycan_5 <- parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
#'
#' has_motif(glycan_4, glycan_5)
#' has_motif(glycan_5, glycan_4)
#' has_motif(glycan_4, glycan_4)
#' has_motif(glycan_5, glycan_5)
#'
#' @export
has_motif <- function(glycan, motif, ..., alignment = "substructure", ignore_linkages = FALSE) {
  # Check input arguments
  if (!glyrepr::is_glycan(glycan) || !glyrepr::is_glycan(motif)) {
    rlang::abort("`glycan` and `motif` must be 'glycan_graph' objects.")
  }
  if (!alignment %in% c("substructure", "core", "terminal", "whole")) {
    rlang::abort("`alignment` must be one of 'substructure', 'core', 'terminal' or 'whole'.")
  }

  # Ensure that `glycan` and `motif` are all "DN" type
  glycan <- glyrepr::convert_graph_mode(glycan, to = "dn", strict = FALSE)
  motif <- glyrepr::convert_graph_mode(motif, to = "dn", strict = FALSE)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycan_mono_type()`.
  glycan <- ensure_glycan_mono_type(glycan, motif)

  # Check if the motif is a subgraph of the glycan
  if (ignore_linkages) {
    glycan <- glyrepr::remove_linkages(glycan)
    motif <- glyrepr::remove_linkages(motif)
    vf2_subgraph_isomorphic(glycan, motif, alignment, ignore_linkages)
  } else {
    # Impute the "?" in the linkage of the motif graph with every possibility
    motif_candidates <- impute_linkages(motif)
    any(purrr::map_lgl(
      motif_candidates,
      ~ vf2_subgraph_isomorphic(glycan, .x, alignment, ignore_linkages)
    ))
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
  get_labels <- function(graph) {
    dplyr::if_else(
      igraph::V(graph)$type == "mono",
      stringr::str_c(igraph::V(graph)$mono, igraph::V(graph)$sub),
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
  if (igraph::vcount(motif) == 1) {
    return(list(motif))
  }
  ne_motif <- glyrepr::convert_graph_mode(motif, to = "ne")
  ne_result <- impute_ne_linkages(ne_motif)
  return(purrr::map(ne_result, glyrepr::convert_graph_mode, to = "dn"))
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


vf2_subgraph_isomorphic <- function(glycan, motif, alignment, ignore_linkage) {
  # Colorize the glycan and motif graphs.
  # This is to add "color" attributes to the vertices and edges of a glycan graph,
  # which is demanded by the "vf2" method of `igraph::subgraph_isomorphic()`.
  colored_graphs <- colorize_glycan_graphs(glycan, motif)
  c_glycan <- colored_graphs[["glycan"]]
  c_motif <- colored_graphs[["motif"]]

  # Perform the subgraph isomorphism test with the "vf2" method.
  # Possible `alignment`: "substructure", "core", "terminal", "whole"
  if (alignment == "whole") {
    return(has_whole_motif(c_glycan, c_motif, ignore_linkage))
  }

  # `res` below is a list of integer vectors, each of which is a mapping
  # from the vertices of the motif graph to the vertices of the glycan graph.
  # For example, for a two-nodes motif, if `res[[1]]` is c(2, 3),
  # that means node 2 in the glycan graph matches node 1 in the motif graph,
  # and node 3 in the glycan graph matches node 2 in the motif graph.
  res <- igraph::graph.get.subisomorphisms.vf2(c_glycan, c_motif)
  if (length(res) == 0) return(FALSE)

  if (alignment == "substructure") {
    return(process_substructure_res(res, glycan, motif, ignore_linkage))
  }
  if (alignment == "core") {
    return(process_core_res(res, glycan, motif, ignore_linkage))
  }
  if (alignment == "terminal") {
    return(process_terminal_res(res, glycan, motif, ignore_linkage))
  }
}


has_whole_motif <- function(c_glycan, c_motif, ignore_linkage) {
  is_isomo <- igraph::isomorphic(c_glycan, c_motif, method = "vf2")
  if (ignore_linkage) {
    is_isomo
  } else {
    is_isomo && match_anomer(c_glycan$anomer, c_motif$anomer)
  }
}


process_substructure_res <- function(res, glycan, motif, ignore_linkage) {
  if (ignore_linkage) {
    return(TRUE)
  } else {
    for (r in res) {
      if (check_motif_anomer(r, glycan, motif)) return(TRUE)
    }
    return(FALSE)
  }
}


check_motif_anomer <- function(r, glycan, motif) {
  # Check if the anomer of the motif is right.
  motif_core_node <- core_node(motif)
  glycan_core_node <- core_node(glycan)
  if (r[[motif_core_node]] == glycan_core_node) {
    # 1. When two core nodes overlap, check if the anomer matches.
    return(match_anomer(glycan$anomer, motif$anomer))
  } else {
    # 2. When motif is within glycan, check if anomer of motif
    # matches the linkage.
    glycan_node_to_motif_core <- r[[motif_core_node]]
    linkage_node <- igraph::neighbors(glycan, glycan_node_to_motif_core, mode = "in")
    linkage <- igraph::V(glycan)[linkage_node]$linkage
    return(match_anomer(stringr::str_sub(linkage, 1, 2), motif$anomer))
  }
}


process_core_res <- function(res, glycan, motif, ignore_linkage) {
  glycan_core_node <- core_node(glycan)
  motif_core_node <- core_node(motif)
  if (!ignore_linkage) {
    if (!match_anomer(glycan$anomer, motif$anomer)) return(FALSE)
  }
  return(any(purrr::map_lgl(res, ~ .x[[motif_core_node]] == glycan_core_node)))
}


process_terminal_res <- function(res, glycan, motif, ignore_linkage) {
  glycan_terminal_nodes <- terminal_nodes(glycan)
  motif_terminal_nodes <- terminal_nodes(motif)
  if (ignore_linkage) {
    return(any(purrr::map_lgl(
      res, ~ all(.x[motif_terminal_nodes] %in% glycan_terminal_nodes)))
    )
  } else {
    return(any(purrr::map_lgl(
      res,
      ~ all(.x[motif_terminal_nodes] %in% glycan_terminal_nodes) &&
        check_motif_anomer(.x, glycan, motif)
    )))
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


match_anomer <- function(glycan_anomer, motif_anomer) {
  # Check if the anomer of the glycan and motif are matched.
  # - "??" in motif will match any anomer and position in glycan.
  # - "?1" in motif will match any anomer at position 1 in glycan.
  # - "a?" in motif will match anomer "a" at any position in glycan.
  # - "a1" in motif will only match anomer "a" at position 1 in glycan.
  anomer_ok <- (
    stringr::str_sub(motif_anomer, 1, 1) == "?" ||
    stringr::str_sub(glycan_anomer, 1, 1) == stringr::str_sub(motif_anomer, 1, 1)
  )
  position_ok <- (
    stringr::str_sub(motif_anomer, 2) == "?" ||
    stringr::str_sub(glycan_anomer, 2) == stringr::str_sub(motif_anomer, 2)
  )
  anomer_ok && position_ok
}
