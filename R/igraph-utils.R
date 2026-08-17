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
  restore_structure_graph_values(structures, graph_mono_type)
}

restore_structure_graph_values <- function(structures, summarize_graph) {
  index <- index_unique_structures(structures)
  values <- vapply(index$graphs, summarize_graph, character(1L))
  result <- values[index$restore]
  names(result) <- names(structures)
  result
}

graph_linkages <- function(graph) {
  linkages <- graph_edge_attr(graph, "linkage")
  floating_parts <- igraph::graph_attr(graph, "floating_parts")
  if (length(floating_parts) > 0L) {
    linkages <- c(
      linkages,
      vapply(floating_parts, \(part) part$linkage, character(1L))
    )
  }
  linkages
}

graph_has_linkages <- function(graph) {
  any(graph_linkages(graph) != "??-?") ||
    igraph::graph_attr(graph, "anomer") != "??"
}

graph_has_strict_linkages <- function(graph) {
  linkages <- graph_linkages(graph)
  anomer <- igraph::graph_attr(graph, "anomer")
  all(!stringr::str_detect(c(linkages, anomer), stringr::fixed("?"))) &&
    all(!stringr::str_detect(linkages, stringr::fixed("/")))
}

structure_level <- function(structures) {
  restore_structure_graph_values(structures, graph_structure_level)
}

graph_structure_level <- function(graph) {
  if (graph_has_strict_linkages(graph)) {
    return("intact")
  }
  if (graph_has_linkages(graph)) {
    return("partial")
  }
  "topological"
}
