# Vertex and edge sequences retain weak references to their graph until the
# graph itself is collected (https://github.com/igraph/rigraph/issues/288).
# These helpers pass plain integer indices so repeated work on a long-lived
# glycan graph does not accumulate those references.
graph_vertex_ids <- function(graph) {
  seq_len(igraph::vcount(graph))
}

graph_edge_ids <- function(graph) {
  seq_len(igraph::ecount(graph))
}

graph_vertex_attr <- function(graph, name, index = NULL) {
  if (is.null(index)) {
    index <- graph_vertex_ids(graph)
  }
  igraph::vertex_attr(graph, name, index = as.integer(index))
}

graph_edge_attr <- function(graph, name, index = NULL) {
  if (is.null(index)) {
    index <- graph_edge_ids(graph)
  }
  igraph::edge_attr(graph, name, index = as.integer(index))
}

graph_degree <- function(graph, mode = "all") {
  mode <- rlang::arg_match(mode, c("all", "out", "in", "total"))
  edges <- igraph::as_edgelist(graph, names = FALSE)
  vertex_count <- igraph::vcount(graph)
  if (nrow(edges) == 0L) {
    return(numeric(vertex_count))
  }

  out <- tabulate(as.integer(edges[, 1]), nbins = vertex_count)
  incoming <- tabulate(as.integer(edges[, 2]), nbins = vertex_count)
  if (!igraph::is_directed(graph) || mode %in% c("all", "total")) {
    return(as.numeric(out + incoming))
  }
  if (mode == "out") {
    return(as.numeric(out))
  }
  as.numeric(incoming)
}

graph_incident_edge_ids <- function(graph, vertex, mode = "all") {
  mode <- rlang::arg_match(mode, c("all", "out", "in"))
  edges <- igraph::as_edgelist(graph, names = FALSE)
  if (nrow(edges) == 0L) {
    return(integer())
  }

  vertex <- as.integer(vertex)
  if (!igraph::is_directed(graph) || mode == "all") {
    return(which(edges[, 1] %in% vertex | edges[, 2] %in% vertex))
  }
  if (mode == "out") {
    return(which(edges[, 1] %in% vertex))
  }
  which(edges[, 2] %in% vertex)
}

graph_neighbor_ids <- function(graph, vertex, mode = "all") {
  mode <- rlang::arg_match(mode, c("all", "out", "in"))
  edges <- igraph::as_edgelist(graph, names = FALSE)
  if (nrow(edges) == 0L) {
    return(integer())
  }

  vertex <- as.integer(vertex)
  if (!igraph::is_directed(graph) || mode == "all") {
    from <- edges[edges[, 2] %in% vertex, 1]
    to <- edges[edges[, 1] %in% vertex, 2]
    return(sort(unique(as.integer(c(from, to)))))
  }
  if (mode == "out") {
    return(sort(unique(as.integer(edges[edges[, 1] %in% vertex, 2]))))
  }
  sort(unique(as.integer(edges[edges[, 2] %in% vertex, 1])))
}

graph_subcomponent_ids <- function(graph, vertex, mode = "all") {
  mode <- rlang::arg_match(mode, c("all", "out", "in"))
  queue <- as.integer(vertex)
  visited <- rep(FALSE, igraph::vcount(graph))
  result <- integer()

  while (length(queue) > 0L) {
    current <- queue[[1]]
    queue <- queue[-1]
    if (visited[[current]]) {
      next
    }

    visited[[current]] <- TRUE
    result <- c(result, current)
    neighbors <- graph_neighbor_ids(graph, current, mode = mode)
    queue <- c(queue, neighbors[!visited[neighbors]])
  }

  result
}

graph_mono_type <- function(graph) {
  monos <- graph_vertex_attr(graph, "mono")
  types <- tryCatch(
    glyrepr::get_mono_type(monos),
    error = function(cnd) "unknown"
  )
  unique_types <- unique(types)
  if (length(unique_types) > 1L) {
    return("mixed")
  }
  unique_types[[1]]
}

structure_mono_type <- function(structures) {
  if (length(structures) == 0L) {
    return(character())
  }
  graphs <- attr(structures, "graphs")
  if (length(graphs) == 0L) {
    return(NA_character_)
  }
  graph_mono_type(graphs[[1]])
}

graph_has_linkages <- function(graph) {
  any(graph_edge_attr(graph, "linkage") != "??-?") ||
    igraph::graph_attr(graph, "anomer") != "??"
}

graph_has_strict_linkages <- function(graph) {
  linkages <- graph_edge_attr(graph, "linkage")
  anomer <- igraph::graph_attr(graph, "anomer")
  all(!stringr::str_detect(c(linkages, anomer), stringr::fixed("?"))) &&
    all(!stringr::str_detect(linkages, stringr::fixed("/")))
}

structure_level <- function(structures) {
  codes <- as.character(structures)
  valid_codes <- unique(codes[!is.na(codes)])
  if (length(valid_codes) == 0L) {
    return(NA_character_)
  }

  graphs <- attr(structures, "graphs")[valid_codes]
  if (graph_mono_type(graphs[[1]]) == "generic") {
    return("basic")
  }

  strict <- vapply(graphs, graph_has_strict_linkages, logical(1L))
  if (all(strict)) {
    return("intact")
  }
  lenient <- vapply(graphs, graph_has_linkages, logical(1L))
  if (any(lenient)) {
    return("partial")
  }
  "topological"
}

graph_subisomorphisms_vf2 <- function(
  graph,
  pattern,
  graph_colors = NULL,
  pattern_colors = NULL
) {
  if (is.null(graph_colors) && is.null(pattern_colors)) {
    if ("color" %in% igraph::vertex_attr_names(graph)) {
      graph_colors <- graph_vertex_attr(graph, "color")
    }
    if ("color" %in% igraph::vertex_attr_names(pattern)) {
      pattern_colors <- graph_vertex_attr(pattern, "color")
    }
  }

  graph_edge_colors <- if ("color" %in% igraph::edge_attr_names(graph)) {
    graph_edge_attr(graph, "color")
  } else {
    NULL
  }
  pattern_edge_colors <- if ("color" %in% igraph::edge_attr_names(pattern)) {
    graph_edge_attr(pattern, "color")
  } else {
    NULL
  }

  # The public wrapper converts every numeric mapping to an igraph vertex
  # sequence. The implementation returns the same zero-based mappings without
  # creating graph weak references.
  impl <- get(
    "get_subisomorphisms_vf2_impl",
    envir = asNamespace("igraph"),
    inherits = FALSE
  )
  mappings <- impl(
    graph1 = graph,
    graph2 = pattern,
    vertex_color1 = graph_colors,
    vertex_color2 = pattern_colors,
    edge_color1 = graph_edge_colors,
    edge_color2 = pattern_edge_colors
  )
  lapply(mappings, \(mapping) as.integer(mapping) + 1L)
}
