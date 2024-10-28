has_motif <- function(glycan, motif, ignore_linkages = FALSE) {
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
    vf2_subgraph_isomorphic(glycan, motif)
  } else {
    # Impute the "?" in the linkage of the motif graph with every possibility
    motif_candidates <- impute_linkages(motif)
    any(purrr::map_lgl(
      motif_candidates,
      ~ vf2_subgraph_isomorphic(glycan, .x)
    ))
  }
}


ensure_motif_graph_mode <- function(glycan, motif) {
  # Make motif the same graph mode (NE or DN) as the glycan
  if (glyrepr::is_ne_glycan(glycan) && glyrepr::is_dn_glycan(motif)) {
    return(glyrepr::convert_dn_to_ne(motif))
  }
  if (glyrepr::is_dn_glycan(glycan) && glyrepr::is_ne_glycan(motif)) {
    return(glyrepr::convert_ne_to_dn(motif))
  }
  return(motif)
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
  candidates <- expand.grid(purrr::map(linkages, glyrepr::possible_linkages), stringsAsFactors = FALSE)
  candidates <- purrr::list_transpose(unclass(candidates))
  candidates <- purrr::map(candidates, ~ unname(.x))
  purrr::map(candidates, ~ igraph::set_edge_attr(motif, "linkage", value = .x))
}


impute_dn_linkages <- function(motif) {
  ne_motif <- glyrepr::convert_dn_to_ne(motif)
  ne_result <- impute_ne_linkages(ne_motif)
  purrr::map(ne_result, ~ glyrepr::convert_ne_to_dn(.x))
}


vf2_subgraph_isomorphic <- function(glycan, motif) {
  # Colorize the glycan and motif graphs.
  # This is to add "color" attributes to the vertices and edges of a glycan graph,
  # which is demanded by the "vf2" method of `igraph::subgraph_isomorphic()`.
  colored_graphs <- colorize_glycan_graphs(glycan, motif)
  c_glycan <- colored_graphs[["glycan"]]
  c_motif <- colored_graphs[["motif"]]
  # Perform the subgraph isomorphism test with the "vf2" method.
  igraph::subgraph_isomorphic(c_motif, c_glycan, method = "vf2")
}
