has_motif <- function(glycan, motif) {
  if (!glyrepr::is_glycan(glycan)) {
    rlang::abort("`glycan` must be a 'glycan_graph' object.")
  }
  if (!glyrepr::is_glycan(motif)) {
    rlang::abort("`motif` must be a 'glycan_graph' object.")
  }
  if (inherits(glycan, "ne_glycan_graph") && inherits(motif, "dn_glycan_graph")) {
      motif <- glyrepr::convert_dn_to_ne(motif)
  }
  if (inherits(glycan, "dn_glycan_graph") && inherits(motif, "ne_glycan_graph")) {
      motif <- glyrepr::convert_ne_to_dn(motif)
  }
  colored_graphs <- colorize_glycan_graphs(glycan, motif)
  igraph::subgraph_isomorphic(colored_graphs[["motif"]], colored_graphs[["glycan"]], method = "vf2")
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
  colors <- colorize_labels(get_labels(glycan), get_labels(motif))
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
