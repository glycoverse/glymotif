source("benchmarks/rcpp-floating-matcher/setup.R")
error <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = identity
  )
}
check_rejection <- function(graph, label) {
  a <- error(cpp_localizations(graph))
  b <- error(enumerate_floating_graph_localizations(graph))
  stopifnot(inherits(a, "error"), inherits(b, "error"))
  data.frame(
    case = label,
    cpp_error = conditionMessage(a),
    r_error = conditionMessage(b),
    both_rejected = TRUE
  )
}
x <- as_glycan_structure("{Man(a1-?)}{Man(a1-?)}Man(a1-")
g <- attr(x, "graphs")[[1]]
p <- igraph::graph_attr(g, "floating_parts")
p[[1]]$parents <- p[[2]]$root
p[[2]]$parents <- p[[1]]$root
g <- igraph::set_graph_attr(g, "floating_parts", p)
rows <- list(check_rejection(g, "component-cycle"))
x <- as_glycan_structure("{6S|2,3}{Fuc(a1-6)|2,3}Gal(a1-3)Glc(a1-")
g <- attr(x, "graphs")[[1]]
p <- igraph::graph_attr(g, "floating_parts")
s <- igraph::graph_attr(g, "floating_substituents")
p[[1]]$parents <- 2L
s[[1]]$parents <- 2L
g <- igraph::set_graph_attr(
  igraph::set_graph_attr(g, "floating_parts", p),
  "floating_substituents",
  s
)
rows[[2]] <- check_rejection(g, "part-substituent-slot-collision")
# Graph-level fixtures avoid canonicalization collapsing symmetric assignments.
make_limit_graph <- function(nparts) {
  n <- nparts + 2L
  graph <- igraph::make_empty_graph(n, directed = TRUE)
  graph <- igraph::add_edges(graph, c(n, n - 1L), linkage = "??-?")
  graph <- igraph::set_vertex_attr(
    graph,
    "name",
    value = as.character(seq_len(n))
  )
  graph <- igraph::set_vertex_attr(graph, "mono", value = rep("Man", n))
  graph <- igraph::set_vertex_attr(graph, "sub", value = rep("", n))
  graph <- igraph::set_graph_attr(graph, "anomer", "??")
  graph <- igraph::set_graph_attr(graph, "alditol", FALSE)
  parts <- lapply(seq_len(nparts), function(i) {
    list(
      root = as.integer(i),
      nodes = as.integer(i),
      linkage = "??-?",
      parents = as.integer(c(n - 1L, n))
    )
  })
  igraph::set_graph_attr(graph, "floating_parts", parts)
}
g <- make_limit_graph(8L)
a <- cpp_localizations(g)
b <- enumerate_floating_graph_localizations(g)
stopifnot(a$raw == 256, length(a$variants) == 256L, nrow(b) == 256L)
for (i in seq_len(256L)) {
  stopifnot(identical(a$variants[[i]]$parents, b$assignments[[i]]$parent_node))
}
rows[[3]] <- data.frame(
  case = "exact-256-accepted",
  cpp_error = "",
  r_error = "",
  both_rejected = FALSE
)
rows[[4]] <- check_rejection(make_limit_graph(9L), "512-rejected")
write.csv(
  do.call(rbind, rows),
  file.path(out, "guard-audit.csv"),
  row.names = FALSE
)
cat("GUARDS DONE\n")
