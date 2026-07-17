prepare_match_batch <- function(glycans, motifs, mode = "strict") {
  glycan_index <- index_unique_structures(glycans)
  motif_index <- index_unique_structures(motifs)

  glycan_graphs <- glycan_index$graphs
  motif_graphs <- motif_index$graphs

  if (identical(glyrepr::get_mono_type(motifs), "generic")) {
    glycan_graphs <- purrr::map(glycan_graphs, convert_graph_to_generic)
  }

  glycan_monos <- purrr::map(
    glycan_graphs,
    igraph::vertex_attr,
    name = "mono"
  )
  motif_monos <- purrr::map(
    motif_graphs,
    igraph::vertex_attr,
    name = "mono"
  )
  exact_keys <- unique(c(
    unlist(glycan_monos, use.names = FALSE),
    unlist(motif_monos, use.names = FALSE)
  ))

  motif_key_modes <- purrr::map_chr(
    motif_graphs,
    resolve_residue_key_mode,
    mode = mode
  )
  base_keys <- if (any(motif_key_modes == "base")) {
    unique(residue_color_keys(exact_keys))
  } else {
    NULL
  }
  generic_keys <- if (any(motif_key_modes == "generic")) {
    generic_key_vectors <- c(
      purrr::map(glycan_monos, residue_match_keys, key_mode = "generic"),
      purrr::map(motif_monos, residue_match_keys, key_mode = "generic")
    )
    unique(unlist(generic_key_vectors, use.names = FALSE))
  } else {
    NULL
  }

  glycan_profiles <- purrr::map2(
    glycan_graphs,
    glycan_monos,
    ~ new_batch_graph_profile(
      .x,
      .y,
      exact_keys = exact_keys,
      base_keys = base_keys,
      generic_keys = generic_keys
    )
  )
  motif_profiles <- purrr::map2(
    motif_graphs,
    seq_along(motif_graphs),
    ~ new_batch_graph_profile(
      .x,
      motif_monos[[.y]],
      exact_keys = exact_keys,
      base_keys = base_keys,
      generic_keys = generic_keys,
      key_mode = motif_key_modes[[.y]],
      include_composition_profile = TRUE,
      mode = mode
    )
  )

  list(
    glycan_profiles = glycan_profiles,
    motif_profiles = motif_profiles,
    glycan_restore = glycan_index$restore,
    motif_restore = motif_index$restore
  )
}

index_unique_structures <- function(structures) {
  codes <- unname(as.character(structures))
  unique_codes <- unique(codes[!is.na(codes)])
  first <- match(unique_codes, codes)
  graphs <- if (length(first) == 0L) {
    list()
  } else {
    glyrepr::get_structure_graphs(structures[first], return_list = TRUE)
  }

  list(
    graphs = graphs,
    restore = match(codes, unique_codes)
  )
}

convert_graph_to_generic <- function(graph) {
  monos <- igraph::vertex_attr(graph, "mono")
  igraph::set_vertex_attr(
    graph,
    "mono",
    value = glyrepr::convert_to_generic(monos)
  )
}

new_batch_graph_profile <- function(
  graph,
  monos,
  exact_keys,
  base_keys = NULL,
  generic_keys = NULL,
  key_mode = "exact",
  include_composition_profile = FALSE,
  mode = "strict"
) {
  exact_colors <- match(monos, exact_keys)
  base_colors <- if (is.null(base_keys)) {
    NULL
  } else {
    match(residue_color_keys(monos), base_keys)
  }
  generic_colors <- if (is.null(generic_keys)) {
    NULL
  } else {
    match(residue_match_keys(monos, "generic"), generic_keys)
  }

  list(
    graph = graph,
    monos = monos,
    subs = igraph::vertex_attr(graph, "sub"),
    vcount = igraph::vcount(graph),
    ecount = igraph::ecount(graph),
    core = core_node(graph),
    has_linkages = graph_has_linkages(graph),
    key_mode = key_mode,
    exact_colors = exact_colors,
    base_colors = base_colors,
    generic_colors = generic_colors,
    exact_counts = tabulate(exact_colors, nbins = length(exact_keys)),
    base_counts = if (is.null(base_colors)) {
      NULL
    } else {
      tabulate(base_colors, nbins = length(base_keys))
    },
    generic_counts = if (is.null(generic_colors)) {
      NULL
    } else {
      tabulate(generic_colors, nbins = length(generic_keys))
    },
    composition_profile = if (include_composition_profile) {
      new_motif_composition_profile(graph, mode = mode)
    } else {
      NULL
    }
  )
}

apply_batch_motif <- function(
  batch,
  motif_position,
  alignment,
  ignore_linkages,
  strict_sub,
  match_degree,
  mode,
  single_glycan_func,
  result_type
) {
  motif_profile <- batch$motif_profiles[[
    batch$motif_restore[[motif_position]]
  ]]

  apply_one <- function(glycan_profile) {
    single_glycan_func(
      glycan_graph = glycan_profile$graph,
      motif_graph = motif_profile$graph,
      motif_has_linkages = motif_profile$has_linkages,
      motif_composition_profile = motif_profile$composition_profile,
      alignment = alignment,
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub,
      match_degree = match_degree,
      mode = mode,
      glycan_batch_profile = glycan_profile,
      motif_batch_profile = motif_profile
    )
  }

  unique_results <- switch(
    result_type,
    logical = vapply(batch$glycan_profiles, apply_one, logical(1L)),
    integer = vapply(batch$glycan_profiles, apply_one, integer(1L)),
    list = lapply(batch$glycan_profiles, apply_one)
  )

  restore_batch_results(
    unique_results,
    batch$glycan_restore,
    result_type
  )
}

restore_batch_results <- function(results, restore, result_type) {
  if (result_type != "list") {
    return(results[restore])
  }

  restored <- vector("list", length(restore))
  valid <- !is.na(restore)
  restored[valid] <- results[restore[valid]]
  restored
}
