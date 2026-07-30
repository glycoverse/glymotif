has_any_floating_parts <- function(structures) {
  any(glyrepr::has_floating_parts(structures), na.rm = TRUE)
}

.max_floating_localizations <- 256L

floating_localization_graphs <- function(structure, graph) {
  if (!isTRUE(glyrepr::has_floating_parts(structure))) {
    return(list(graph))
  }

  localizations <- glyrepr::enumerate_floating_localizations(
    structure,
    max_variants = .max_floating_localizations,
    deduplicate = FALSE
  )
  candidate_edges <- glyrepr::structure_candidate_edges(structure)

  purrr::map(
    localizations$assignments,
    add_floating_assignment_edges,
    graph = graph,
    candidate_edges = candidate_edges
  )
}

add_floating_assignment_edges <- function(
  assignments,
  graph,
  candidate_edges
) {
  edge_rows <- purrr::map_int(
    seq_len(nrow(assignments)),
    function(i) {
      matches <- which(
        candidate_edges$part_id == assignments$part_id[[i]] &
          candidate_edges$from_node == assignments$parent_node[[i]]
      )
      if (length(matches) != 1L) {
        cli::cli_abort(
          "Unable to resolve a floating-part localization edge."
        )
      }
      matches
    }
  )
  edges <- candidate_edges[edge_rows, , drop = FALSE]

  igraph::add_edges(
    graph,
    as.integer(t(cbind(edges$from_node, edges$to_node))),
    attr = list(linkage = edges$linkage)
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
  codes <- unname(as.character(structures))
  first <- match(unique(codes[!is.na(codes)]), codes)
  index$graphs <- purrr::map2(
    structures[first],
    index$graphs,
    floating_localization_graphs
  )
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
