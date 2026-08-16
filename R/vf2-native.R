graph_edge_matrix <- function(graph) {
  edges <- igraph::as_edgelist(graph, names = FALSE)
  matrix(as.integer(edges), ncol = 2L)
}


new_native_graph_profile <- function(graph, monos = NULL, subs = NULL) {
  if (is.null(monos)) {
    monos <- graph_vertex_attr(graph, "mono")
  }
  if (is.null(subs)) {
    subs <- graph_vertex_attr(graph, "sub")
  }

  edges <- graph_edge_matrix(graph)
  linkages <- graph_edge_attr(graph, "linkage")
  in_degree <- graph_degree(graph, mode = "in")
  out_degree <- graph_degree(graph, mode = "out")
  incoming_anomers <- rep(NA_character_, igraph::vcount(graph))
  if (nrow(edges) > 0L) {
    incoming_anomers[edges[, 2]] <- stringr::str_extract(
      linkages,
      "^[^-]+"
    )
  }
  core <- which(in_degree == 0)
  incoming_anomers[core] <- igraph::graph_attr(graph, "anomer")

  list(
    edges = edges,
    residue_keys = paste(monos, subs, sep = "\r"),
    linkages = linkages,
    node_anomers = incoming_anomers,
    core = core,
    terminals = which(out_degree == 0),
    in_degree = in_degree,
    out_degree = out_degree
  )
}


compatibility_matrix_by_key <- function(
  glycan_keys,
  motif_keys,
  predicate
) {
  glycan_count <- length(glycan_keys)
  motif_count <- length(motif_keys)
  if (glycan_count == 0L || motif_count == 0L) {
    return(matrix(FALSE, nrow = motif_count, ncol = glycan_count))
  }

  glycan_unique <- unique(glycan_keys)
  motif_unique <- unique(motif_keys)
  glycan_first <- match(glycan_unique, glycan_keys)
  motif_first <- match(motif_unique, motif_keys)
  class_compatibility <- matrix(
    FALSE,
    nrow = length(motif_unique),
    ncol = length(glycan_unique)
  )

  for (motif_class in seq_along(motif_unique)) {
    for (glycan_class in seq_along(glycan_unique)) {
      class_compatibility[motif_class, glycan_class] <- predicate(
        glycan_first[[glycan_class]],
        motif_first[[motif_class]]
      )
    }
  }

  class_compatibility[
    match(motif_keys, motif_unique),
    match(glycan_keys, glycan_unique),
    drop = FALSE
  ]
}


new_native_compatibility_cache <- function(
  glycan_profiles,
  motif_profiles,
  strict_sub,
  ignore_linkages,
  mode
) {
  glycan_monos <- unlist(
    lapply(glycan_profiles, `[[`, "monos"),
    use.names = FALSE
  )
  motif_monos <- unlist(
    lapply(motif_profiles, `[[`, "monos"),
    use.names = FALSE
  )
  glycan_subs <- unlist(
    lapply(glycan_profiles, `[[`, "subs"),
    use.names = FALSE
  )
  motif_subs <- unlist(
    lapply(motif_profiles, `[[`, "subs"),
    use.names = FALSE
  )
  glycan_residue_keys <- paste(glycan_monos, glycan_subs, sep = "\r")
  motif_residue_keys <- paste(motif_monos, motif_subs, sep = "\r")
  unique_glycan_residue_keys <- unique(glycan_residue_keys)
  unique_motif_residue_keys <- unique(motif_residue_keys)
  glycan_first <- match(unique_glycan_residue_keys, glycan_residue_keys)
  motif_first <- match(unique_motif_residue_keys, motif_residue_keys)
  unique_glycan_monos <- glycan_monos[glycan_first]
  unique_motif_monos <- motif_monos[motif_first]
  unique_glycan_subs <- glycan_subs[glycan_first]
  unique_motif_subs <- motif_subs[motif_first]
  residue_compatibility <- compatibility_matrix_by_key(
    unique_glycan_residue_keys,
    unique_motif_residue_keys,
    function(glycan_index, motif_index) {
      match_residue(
        unique_glycan_monos[[glycan_index]],
        unique_glycan_subs[[glycan_index]],
        unique_motif_monos[[motif_index]],
        unique_motif_subs[[motif_index]],
        strict_sub = strict_sub,
        mode = mode
      )
    }
  )

  glycan_linkages <- unique(unlist(
    lapply(glycan_profiles, \(profile) profile$native$linkages),
    use.names = FALSE
  ))
  motif_linkages <- unique(unlist(
    lapply(motif_profiles, \(profile) profile$native$linkages),
    use.names = FALSE
  ))
  linkage_compatibility <- if (ignore_linkages) {
    NULL
  } else {
    compatibility_matrix_by_key(
      glycan_linkages,
      motif_linkages,
      function(glycan_index, motif_index) {
        match_linkage(
          glycan_linkages[[glycan_index]],
          motif_linkages[[motif_index]],
          mode = mode
        )
      }
    )
  }

  glycan_anomers <- unique(unlist(
    lapply(glycan_profiles, \(profile) profile$native$node_anomers),
    use.names = FALSE
  ))
  motif_anomers <- unique(vapply(
    motif_profiles,
    function(profile) {
      profile$native$node_anomers[[profile$native$core]]
    },
    character(1L)
  ))
  anomer_compatibility <- if (ignore_linkages) {
    NULL
  } else {
    compatibility_matrix_by_key(
      glycan_anomers,
      motif_anomers,
      function(glycan_index, motif_index) {
        match_anomer(
          glycan_anomers[[glycan_index]],
          motif_anomers[[motif_index]],
          mode = mode
        )
      }
    )
  }

  list(
    glycan_residue_keys = unique_glycan_residue_keys,
    motif_residue_keys = unique_motif_residue_keys,
    residue_compatibility = residue_compatibility,
    glycan_linkages = glycan_linkages,
    motif_linkages = motif_linkages,
    linkage_compatibility = linkage_compatibility,
    glycan_anomers = glycan_anomers,
    motif_anomers = motif_anomers,
    anomer_compatibility = anomer_compatibility
  )
}


residue_compatibility_matrix <- function(
  context,
  strict_sub,
  mode
) {
  glycan_keys <- paste(context$glycan_mono, context$glycan_sub, sep = "\r")
  motif_keys <- paste(context$motif_mono, context$motif_sub, sep = "\r")
  compatibility_matrix_by_key(
    glycan_keys,
    motif_keys,
    function(glycan_index, motif_index) {
      match_residue(
        context$glycan_mono[[glycan_index]],
        context$glycan_sub[[glycan_index]],
        context$motif_mono[[motif_index]],
        context$motif_sub[[motif_index]],
        strict_sub = strict_sub,
        mode = mode
      )
    }
  )
}


add_alignment_compatibility <- function(
  compatibility,
  context,
  alignment
) {
  if (alignment == "core") {
    incompatible <- setdiff(
      seq_len(ncol(compatibility)),
      context$glycan_core
    )
    compatibility[context$motif_core, incompatible] <- FALSE
  }

  if (alignment == "terminal") {
    incompatible <- setdiff(
      seq_len(ncol(compatibility)),
      context$glycan_terminals
    )
    compatibility[context$motif_terminals, incompatible] <- FALSE
  }

  compatibility
}


add_degree_compatibility <- function(
  compatibility,
  context,
  match_degree
) {
  if (is.null(match_degree) || !any(match_degree)) {
    return(compatibility)
  }

  for (motif_index in which(match_degree)) {
    compatible <-
      context$motif_in[[motif_index]] == context$glycan_in &
      context$motif_out[[motif_index]] == context$glycan_out
    compatibility[motif_index, !compatible] <- FALSE
  }
  compatibility
}


glycan_node_anomers <- function(context) {
  vapply(
    seq_along(context$glycan_mono),
    function(glycan_index) {
      if (glycan_index %in% context$glycan_core) {
        return(context$glycan_anomer)
      }
      stringr::str_split_1(
        context$glycan_incoming_linkages[[glycan_index]],
        stringr::fixed("-")
      )[[1]]
    },
    character(1L)
  )
}


add_anomer_compatibility <- function(
  compatibility,
  context,
  mode
) {
  anomers <- glycan_node_anomers(context)
  compatible <- vapply(
    anomers,
    match_anomer,
    motif_anomer = context$motif_anomer,
    mode = mode,
    logical(1L)
  )
  compatibility[context$motif_core, !compatible] <- FALSE
  compatibility
}


edge_compatibility_matrix <- function(
  glycan,
  motif,
  ignore_linkages,
  mode
) {
  glycan_count <- igraph::ecount(glycan)
  motif_count <- igraph::ecount(motif)
  if (ignore_linkages) {
    return(matrix(TRUE, nrow = motif_count, ncol = glycan_count))
  }

  glycan_linkages <- graph_edge_attr(glycan, "linkage")
  motif_linkages <- graph_edge_attr(motif, "linkage")
  compatibility_matrix_by_key(
    glycan_linkages,
    motif_linkages,
    function(glycan_index, motif_index) {
      match_linkage(
        glycan_linkages[[glycan_index]],
        motif_linkages[[motif_index]],
        mode = mode
      )
    }
  )
}


prepare_vf2_compatibility_from_profiles <- function(
  glycan_profile,
  motif_profile,
  compatibility_cache,
  alignment,
  ignore_linkages,
  match_degree
) {
  glycan_native <- glycan_profile$native
  motif_native <- motif_profile$native
  glycan_residue_ids <- match(
    glycan_native$residue_keys,
    compatibility_cache$glycan_residue_keys
  )
  motif_residue_ids <- match(
    motif_native$residue_keys,
    compatibility_cache$motif_residue_keys
  )
  vertex <- compatibility_cache$residue_compatibility[
    motif_residue_ids,
    glycan_residue_ids,
    drop = FALSE
  ]

  context <- list(
    glycan_core = glycan_native$core,
    motif_core = motif_native$core,
    glycan_terminals = glycan_native$terminals,
    motif_terminals = motif_native$terminals,
    glycan_in = glycan_native$in_degree,
    motif_in = motif_native$in_degree,
    glycan_out = glycan_native$out_degree,
    motif_out = motif_native$out_degree
  )
  if (is.null(match_degree)) {
    vertex <- add_alignment_compatibility(vertex, context, alignment)
  }
  vertex <- add_degree_compatibility(vertex, context, match_degree)

  if (!ignore_linkages) {
    glycan_anomer_ids <- match(
      glycan_native$node_anomers,
      compatibility_cache$glycan_anomers
    )
    motif_anomer_id <- match(
      motif_native$node_anomers[[motif_native$core]],
      compatibility_cache$motif_anomers
    )
    compatible_anomers <- compatibility_cache$anomer_compatibility[
      motif_anomer_id,
      glycan_anomer_ids
    ]
    vertex[motif_native$core, !compatible_anomers] <- FALSE
  }

  edge <- if (ignore_linkages) {
    matrix(
      TRUE,
      nrow = length(motif_native$linkages),
      ncol = length(glycan_native$linkages)
    )
  } else {
    compatibility_cache$linkage_compatibility[
      match(motif_native$linkages, compatibility_cache$motif_linkages),
      match(glycan_native$linkages, compatibility_cache$glycan_linkages),
      drop = FALSE
    ]
  }

  list(vertex = vertex, edge = edge)
}


prepare_vf2_compatibility <- function(
  glycan,
  motif,
  alignment,
  ignore_linkages,
  strict_sub,
  match_degree,
  mode
) {
  context <- prepare_validation_context(
    glycan,
    motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree
  )
  vertex <- residue_compatibility_matrix(
    context,
    strict_sub = strict_sub,
    mode = mode
  )
  if (is.null(match_degree)) {
    vertex <- add_alignment_compatibility(vertex, context, alignment)
  }
  vertex <- add_degree_compatibility(vertex, context, match_degree)
  if (!ignore_linkages) {
    vertex <- add_anomer_compatibility(vertex, context, mode)
  }

  list(
    vertex = vertex,
    edge = edge_compatibility_matrix(
      glycan,
      motif,
      ignore_linkages = ignore_linkages,
      mode = mode
    )
  )
}


perform_compatible_vf2 <- function(
  glycan,
  motif,
  alignment,
  ignore_linkages,
  strict_sub,
  match_degree,
  mode,
  first_only = FALSE,
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL,
  compatibility_cache = NULL
) {
  use_batch_profiles <-
    !is.null(glycan_batch_profile) &&
    !is.null(motif_batch_profile) &&
    !is.null(compatibility_cache)
  compatibility <- if (use_batch_profiles) {
    prepare_vf2_compatibility_from_profiles(
      glycan_batch_profile,
      motif_batch_profile,
      compatibility_cache,
      alignment = alignment,
      ignore_linkages = ignore_linkages,
      match_degree = match_degree
    )
  } else {
    prepare_vf2_compatibility(
      glycan,
      motif,
      alignment = alignment,
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub,
      match_degree = match_degree,
      mode = mode
    )
  }
  glycan_edges <- if (use_batch_profiles) {
    glycan_batch_profile$native$edges
  } else {
    graph_edge_matrix(glycan)
  }
  motif_edges <- if (use_batch_profiles) {
    motif_batch_profile$native$edges
  } else {
    graph_edge_matrix(motif)
  }

  cpp_vf2_subgraph_mono(
    glycan_vertex_count = igraph::vcount(glycan),
    glycan_edges = glycan_edges,
    motif_vertex_count = igraph::vcount(motif),
    motif_edges = motif_edges,
    vertex_compatibility = compatibility$vertex,
    edge_compatibility = compatibility$edge,
    first_only = first_only,
    unique_vertex_sets = !first_only
  )
}
