has_any_floating_metadata <- function(structures) {
  has_floating <- glyrepr::has_floating_parts(structures) |
    glyrepr::has_floating_substituents(structures)
  any(has_floating, na.rm = TRUE)
}

.max_floating_localizations <- 256L

graph_has_floating_metadata <- function(graph) {
  any(
    c("floating_parts", "floating_substituents") %in%
      igraph::graph_attr_names(graph)
  )
}

floating_localization_graphs <- function(graph) {
  if (!graph_has_floating_metadata(graph)) {
    return(list(graph))
  }

  glyrepr::enumerate_floating_graph_localizations(
    graph,
    max_variants = .max_floating_localizations
  )$graph
}

aggregate_floating_graph_results <- function(
  graph,
  .f,
  result_type,
  strict_floating
) {
  if (!graph_has_floating_metadata(graph)) {
    return(.f(graph))
  }

  results <- lapply(
    floating_localization_graphs(graph),
    .f
  )
  aggregate_floating_results(
    results,
    result_type = result_type,
    strict_floating = strict_floating
  )
}

map_floating_structures <- function(
  structures,
  .f,
  ...,
  result_type,
  strict_floating,
  localization_index = NULL
) {
  if (is.null(localization_index)) {
    localization_index <- prepare_floating_graph_index(structures)
  }

  apply_one <- function(i) {
    graphs <- localization_index$graphs[[i]]
    variant_results <- lapply(graphs, .f, ...)
    aggregate_floating_results(
      variant_results,
      result_type = result_type,
      strict_floating = strict_floating
    )
  }

  unique_results <- switch(
    result_type,
    logical = vapply(
      seq_along(localization_index$graphs),
      apply_one,
      logical(1L)
    ),
    integer = vapply(
      seq_along(localization_index$graphs),
      apply_one,
      integer(1L)
    ),
    list = lapply(seq_along(localization_index$graphs), apply_one)
  )

  restore_batch_results(
    unique_results,
    localization_index$restore,
    result_type
  )
}

prepare_floating_graph_index <- function(structures) {
  index <- index_unique_structures(structures)
  index$graphs <- purrr::map(index$graphs, floating_localization_graphs)
  index
}

aggregate_floating_results <- function(
  results,
  result_type,
  strict_floating
) {
  switch(
    result_type,
    logical = if (strict_floating) {
      all(unlist(results))
    } else {
      any(unlist(results))
    },
    integer = if (strict_floating) {
      min(unlist(results))
    } else {
      max(unlist(results))
    },
    list = unique_floating_matches(results)
  )
}

unique_floating_matches <- function(results) {
  matches <- unlist(results, recursive = FALSE)
  if (length(matches) <= 1L) {
    return(matches)
  }

  keys <- vapply(
    matches,
    function(x) paste(x, collapse = ","),
    character(1L),
    USE.NAMES = FALSE
  )
  unname(matches[!duplicated(keys)])
}
