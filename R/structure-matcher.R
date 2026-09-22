# Use public glyrepr/igraph accessors at the representation boundary. Native
# matching owns these compact profiles and never inspects igraph internals.
native_structure_input <- function(structures) {
  index <- index_unique_structures(structures)
  index$graphs <- lapply(index$graphs, function(graph) {
    list(
      n = as.integer(igraph::vcount(graph)),
      edges = graph_edge_matrix(graph),
      attributes = igraph::graph_attr(graph),
      vertices = igraph::vertex_attr(graph),
      edge_attributes = igraph::edge_attr(graph)
    )
  })
  index$codes <- unname(as.character(structures))
  index
}

match_structures_native <- function(
  glycans,
  motifs,
  alignments,
  ignore_linkages,
  strict_sub,
  match_degree,
  mode,
  strict_floating,
  result_type
) {
  if (is.null(match_degree)) {
    match_degree <- rep(list(NULL), length(motifs))
  }
  cpp_match_structures(
    native_structure_input(glycans),
    native_structure_input(motifs),
    native_monosaccharide_dictionary(),
    alignments,
    ignore_linkages,
    strict_sub,
    mode == "lenient",
    match_degree,
    switch(result_type, logical = "have", integer = "count", list = "match"),
    strict_floating,
    .max_floating_localizations
  )
}

# Build the compatibility dictionary through glyrepr's public residue API.
native_monosaccharide_dictionary <- function() {
  monos <- glyrepr::available_monosaccharides()
  data.frame(concrete = monos, generic = glyrepr::convert_to_generic(monos))
}
