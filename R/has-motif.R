has_motif <- function(glycan, motif, ignore_linkages = FALSE) {
  # Check input arguments
  if (!glyrepr::is_glycan(glycan) || !glyrepr::is_glycan(motif)) {
    rlang::abort("`glycan` and `motif` must be 'glycan_graph' objects.")
  }

  # Ensure that `glycan` and `motif` have the same mode
  motif <- ensure_same_mode(glycan, motif)

  # Colorize the glycan and motif graphs.
  # This is to add "color" attributes to the vertices and edges of a glycan graph,
  # which is demanded by the "vf2" method of `igraph::subgraph_isomorphic()`.
  colored_graphs <- colorize_glycan_graphs(glycan, motif, ignore_linkages)

  # Check if the motif is a subgraph of the glycan
  vf2_subgraph_isomorphic(colored_graphs[["glycan"]], colored_graphs[["motif"]])
}


ensure_same_mode <- function(glycan, motif) {
  # Make motif the same mode as the glycan
  if (glyrepr::is_ne_glycan(glycan) && glyrepr::is_dn_glycan(motif)) {
    return(glyrepr::convert_dn_to_ne(motif))
  }
  if (glyrepr::is_dn_glycan(glycan) && glyrepr::is_ne_glycan(motif)) {
    return(glyrepr::convert_ne_to_dn(motif))
  }
  return(motif)
}


colorize_glycan_graphs <- function(glycan, motif, ignore_linkages) {
  # Add "color" attributes to the vertices and edges of a glycan graph.
  # The "color" attrs wiil be usd by `igraph::subgraph_isomorphic()`
  # to compare glycan graphs (with method = "vf2").
  #
  # This function assumes that `glycan` and `input` are both "NE" or "DN" glycan graphs.
  if (inherits(glycan, "ne_glycan_graph")) {
    colorize_ne_glycan_graphs(glycan, motif, ignore_linkages)
  } else {
    colorize_dn_glycan_graphs(glycan, motif, ignore_linkages)
  }
}


colorize_ne_glycan_graphs <- function(glycan, motif, ignore_linkages) {
  node_colors <- colorize_labels(igraph::V(glycan)$mono, igraph::V(motif)$mono)
  igraph::V(glycan)$color <- node_colors[["glycan"]]
  igraph::V(motif)$color <- node_colors[["motif"]]

  if (!ignore_linkages) {
    edge_colors <- colorize_labels(igraph::E(glycan)$linkage, igraph::E(motif)$linkage)
    igraph::E(glycan)$color <- edge_colors[["glycan"]]
    igraph::E(motif)$color <- edge_colors[["motif"]]
  }

  list(glycan = glycan, motif = motif)
}


colorize_dn_glycan_graphs <- function(glycan, motif, ignore_linkages) {
  get_labels <- function(graph, ignore_linkages) {
    dplyr::if_else(
      igraph::V(graph)$type == "mono",
      igraph::V(graph)$mono,
      if (ignore_linkages) "??-?" else igraph::V(graph)$linkage
    )
  }
  glycan_labels <- get_labels(glycan, ignore_linkages)
  motif_labels <- get_labels(motif, ignore_linkages)
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


vf2_subgraph_isomorphic <- function(glycan, motif) {
  # `glycan` and `motif` should be colored.
  igraph::subgraph_isomorphic(motif, glycan, method = "vf2")
}
