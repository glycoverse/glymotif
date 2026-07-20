colorize_graphs <- function(
  glycan,
  motif,
  mode = "strict",
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL
) {
  if (!is.null(glycan_batch_profile) && !is.null(motif_batch_profile)) {
    key_mode <- motif_batch_profile$key_mode
    if (key_mode == "none") {
      return(list(glycan_colors = NULL, motif_colors = NULL))
    }
    if (key_mode == "exact") {
      return(list(
        glycan_colors = glycan_batch_profile$exact_colors,
        motif_colors = motif_batch_profile$exact_colors
      ))
    }
    if (key_mode == "base") {
      return(list(
        glycan_colors = glycan_batch_profile$base_colors,
        motif_colors = motif_batch_profile$base_colors
      ))
    }
    return(list(
      glycan_colors = glycan_batch_profile$generic_colors,
      motif_colors = motif_batch_profile$generic_colors
    ))
  }

  # Prepare VF2 color vectors from the "mono" vertex attributes without
  # mutating the input graphs.
  glycan_monos <- graph_vertex_attr(glycan, "mono")
  motif_monos <- graph_vertex_attr(motif, "mono")
  key_mode <- resolve_residue_key_mode(motif, mode)
  if (key_mode == "none") {
    return(list(glycan_colors = NULL, motif_colors = NULL))
  }
  glycan_monos <- residue_match_keys(glycan_monos, key_mode)
  motif_monos <- residue_match_keys(motif_monos, key_mode)

  unique_monos <- unique(c(glycan_monos, motif_monos))
  color_map <- seq_along(unique_monos)
  names(color_map) <- unique_monos
  glycan_colors <- color_map[glycan_monos]
  motif_colors <- color_map[motif_monos]
  names(glycan_colors) <- NULL
  names(motif_colors) <- NULL
  list(glycan_colors = glycan_colors, motif_colors = motif_colors)
}

#' Create a Motif Composition Profile
#'
#' @param motif A motif graph.
#' @param mode Matching mode.
#'
#' @return A list containing residue keys, required counts, and keying mode.
#' @noRd
new_motif_composition_profile <- function(motif, mode = "strict") {
  key_mode <- resolve_residue_key_mode(motif, mode)
  if (key_mode == "none") {
    return(list(
      keys = character(),
      counts = integer(),
      key_mode = key_mode
    ))
  }

  motif_monos <- residue_match_keys(
    graph_vertex_attr(motif, "mono"),
    key_mode
  )
  keys <- unique(motif_monos)
  counts <- tabulate(match(motif_monos, keys), nbins = length(keys))
  list(
    keys = keys,
    counts = counts,
    key_mode = key_mode
  )
}

#' Check Whether Glycan Composition Can Contain a Motif
#'
#' This is a conservative negative filter. `TRUE` means VF2 is still required;
#' `FALSE` means the glycan lacks enough residues to contain the motif.
#'
#' @param glycan A glycan graph.
#' @param motif_profile A motif composition profile.
#'
#' @return A logical scalar.
#' @noRd
composition_can_match <- function(
  glycan,
  motif_profile,
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL
) {
  key_mode <- motif_profile$key_mode
  if (key_mode == "none") {
    return(TRUE)
  }

  if (!is.null(glycan_batch_profile) && !is.null(motif_batch_profile)) {
    if (key_mode == "exact") {
      return(all(
        glycan_batch_profile$exact_counts >= motif_batch_profile$exact_counts
      ))
    }
    if (key_mode == "base") {
      return(all(
        glycan_batch_profile$base_counts >= motif_batch_profile$base_counts
      ))
    }
    return(all(
      glycan_batch_profile$generic_counts >= motif_batch_profile$generic_counts
    ))
  }

  glycan_monos <- residue_match_keys(
    graph_vertex_attr(glycan, "mono"),
    key_mode
  )

  glycan_counts <- tabulate(
    match(glycan_monos, motif_profile$keys),
    nbins = length(motif_profile$keys)
  )
  all(glycan_counts >= motif_profile$counts)
}


#' Resolve Residue Keys for Conservative VF2 Pruning
#'
#' @param motif A motif graph.
#' @param mode Matching mode.
#'
#' @return One of `"exact"`, `"base"`, `"generic"`, or `"none"`.
#' @noRd
resolve_residue_key_mode <- function(motif, mode = "strict") {
  fuzzy <- has_fuzzy_modification(motif)
  if (mode == "lenient") {
    if (fuzzy) {
      return("none")
    }
    return("generic")
  }

  if (fuzzy) {
    "base"
  } else {
    "exact"
  }
}


#' Create Residue Keys for Conservative VF2 Pruning
#'
#' @param monos A character vector of monosaccharide names from one graph.
#' @param key_mode The residue keying mode.
#'
#' @return A character vector of residue keys.
#' @noRd
residue_match_keys <- function(monos, key_mode) {
  switch(
    key_mode,
    exact = monos,
    base = residue_color_keys(monos),
    generic = glyrepr::convert_to_generic(monos),
    none = character()
  )
}


#' Check Whether Whole-Alignment Graph Sizes Can Match
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param alignment Alignment mode.
#'
#' @return A logical scalar.
#' @noRd
whole_alignment_size_can_match <- function(
  glycan,
  motif,
  alignment,
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL
) {
  if (alignment != "whole") {
    return(TRUE)
  }

  if (!is.null(glycan_batch_profile) && !is.null(motif_batch_profile)) {
    return(
      glycan_batch_profile$vcount == motif_batch_profile$vcount &&
        glycan_batch_profile$ecount == motif_batch_profile$ecount
    )
  }

  igraph::vcount(glycan) == igraph::vcount(motif) &&
    igraph::ecount(glycan) == igraph::ecount(motif)
}


#' Check Whether Core-Alignment Root Residues Can Match
#'
#' This is a conservative negative filter. `TRUE` means VF2 is still required;
#' `FALSE` means the glycan and motif core residues are incompatible.
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param alignment Alignment mode.
#' @param strict_sub Whether substituent matching should be strict.
#'
#' @return A logical scalar.
#' @noRd
core_alignment_root_can_match <- function(
  glycan,
  motif,
  alignment,
  strict_sub,
  mode = "strict",
  glycan_batch_profile = NULL,
  motif_batch_profile = NULL
) {
  if (alignment != "core") {
    return(TRUE)
  }

  if (!is.null(glycan_batch_profile) && !is.null(motif_batch_profile)) {
    glycan_core <- glycan_batch_profile$core
    motif_core <- motif_batch_profile$core
    glycan_mono <- glycan_batch_profile$monos[[glycan_core]]
    glycan_sub <- glycan_batch_profile$subs[[glycan_core]]
    motif_mono <- motif_batch_profile$monos[[motif_core]]
    motif_sub <- motif_batch_profile$subs[[motif_core]]
  } else {
    glycan_core <- core_node(glycan)
    motif_core <- core_node(motif)
    glycan_mono <- igraph::vertex_attr(glycan, "mono", index = glycan_core)
    glycan_sub <- igraph::vertex_attr(glycan, "sub", index = glycan_core)
    motif_mono <- igraph::vertex_attr(motif, "mono", index = motif_core)
    motif_sub <- igraph::vertex_attr(motif, "sub", index = motif_core)
  }

  match_residue(
    glycan_mono,
    glycan_sub,
    motif_mono,
    motif_sub,
    strict_sub = strict_sub,
    mode = mode
  )
}


#' Check Whether a Motif Has Fuzzy Built-In Modifications
#'
#' @param motif A motif graph.
#'
#' @return A logical scalar.
#' @noRd
has_fuzzy_modification <- function(motif) {
  any(purrr::map_lgl(graph_vertex_attr(motif, "sub"), is_fuzzy_sub))
}


#' Create Residue Color Keys for Fuzzy Modification Matching
#'
#' @param monos A character vector of monosaccharide names.
#'
#' @return A character vector of color keys.
#' @noRd
residue_color_keys <- function(monos) {
  purrr::map_chr(monos, function(mono) {
    built_in <- decompose_builtin_modification(mono)
    if (is.null(built_in)) {
      mono
    } else {
      built_in$mono
    }
  })
}


#' Create a Base Validation Context
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#'
#' @return A base validation context.
#' @noRd
new_validation_context <- function(glycan, motif) {
  list(
    glycan = glycan,
    motif = motif,
    glycan_mono = graph_vertex_attr(glycan, "mono"),
    motif_mono = graph_vertex_attr(motif, "mono"),
    glycan_sub = graph_vertex_attr(glycan, "sub"),
    motif_sub = graph_vertex_attr(motif, "sub"),
    glycan_in = NULL,
    glycan_out = NULL,
    motif_in = NULL,
    motif_out = NULL,
    glycan_core = NULL,
    motif_core = NULL,
    glycan_terminals = NULL,
    motif_terminals = NULL,
    glycan_v = NULL,
    motif_v = NULL,
    glycan_e = NULL,
    motif_e = NULL,
    motif_edge_list = NULL,
    motif_linkages = NULL,
    glycan_edge_linkages = NULL,
    glycan_anomer = NULL,
    motif_anomer = NULL,
    glycan_incoming_linkages = NULL
  )
}


#' Add Alignment Data to a Validation Context
#'
#' @param context A validation context.
#' @param alignment Alignment mode.
#'
#' @return A validation context with alignment-specific data.
#' @noRd
add_alignment_context <- function(context, alignment) {
  switch(
    alignment,
    "substructure" = context,
    "core" = add_core_context(context),
    "terminal" = {
      context$glycan_terminals <- which(
        graph_degree(context$glycan, mode = "out") == 0
      )
      context$motif_terminals <- which(
        graph_degree(context$motif, mode = "out") == 0
      )
      context
    },
    "whole" = {
      context$glycan_v <- igraph::vcount(context$glycan)
      context$motif_v <- igraph::vcount(context$motif)
      context$glycan_e <- igraph::ecount(context$glycan)
      context$motif_e <- igraph::ecount(context$motif)
      context
    }
  )
}


#' Add Degree Data to a Validation Context
#'
#' @param context A validation context.
#'
#' @return A validation context with degree vectors.
#' @noRd
add_degree_context <- function(context) {
  context$motif_in <- graph_degree(context$motif, mode = "in")
  context$motif_out <- graph_degree(context$motif, mode = "out")
  context$glycan_in <- graph_degree(context$glycan, mode = "in")
  context$glycan_out <- graph_degree(context$glycan, mode = "out")
  context
}


#' Add Core Node Data to a Validation Context
#'
#' @param context A validation context.
#'
#' @return A validation context with core node indices.
#' @noRd
add_core_context <- function(context) {
  if (is.null(context$glycan_core)) {
    context$glycan_core <- which(
      graph_degree(context$glycan, mode = "in") == 0
    )
    context$motif_core <- which(graph_degree(context$motif, mode = "in") == 0)
  }
  context
}


#' Add Linkage and Anomer Data to a Validation Context
#'
#' @param context A validation context.
#'
#' @return A validation context with linkage and anomer data.
#' @noRd
add_linkage_context <- function(context) {
  context <- add_core_context(context)

  context$glycan_anomer <- igraph::graph_attr(context$glycan, "anomer")
  context$motif_anomer <- igraph::graph_attr(context$motif, "anomer")

  glycan_edge_list <- igraph::as_edgelist(context$glycan, names = FALSE)
  glycan_linkages <- graph_edge_attr(context$glycan, "linkage")
  context$glycan_incoming_linkages <- rep(
    NA_character_,
    igraph::vcount(context$glycan)
  )
  if (nrow(glycan_edge_list) > 0) {
    context$glycan_incoming_linkages[glycan_edge_list[, 2]] <- glycan_linkages
  }

  motif_edge_list <- igraph::as_edgelist(context$motif, names = FALSE)
  context$motif_edge_list <- motif_edge_list
  context$motif_linkages <- graph_edge_attr(context$motif, "linkage")
  if (nrow(motif_edge_list) > 0) {
    context$glycan_edge_linkages <- matrix(
      NA_character_,
      nrow = igraph::vcount(context$glycan),
      ncol = igraph::vcount(context$glycan)
    )
    context$glycan_edge_linkages[glycan_edge_list] <- glycan_linkages
  }

  context
}


#' Prepare Reusable Graph Data for Candidate Validation
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param alignment Alignment mode.
#' @param ignore_linkages Whether linkage/anomer checks are skipped.
#' @param match_degree Degree matching vector.
#'
#' @return A validation context containing only data needed for this validation mode.
#' @noRd
prepare_validation_context <- function(
  glycan,
  motif,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL
) {
  context <- new_validation_context(glycan, motif)
  context <- add_alignment_context(context, alignment)
  if (!is.null(match_degree) && any(match_degree)) {
    context <- add_degree_context(context)
  }
  if (!ignore_linkages) {
    context <- add_linkage_context(context)
  }
  context
}


perform_vf2 <- function(
  glycan,
  motif,
  glycan_colors = NULL,
  motif_colors = NULL
) {
  # Perform "VF2" algorithm
  # `res` is a list of all possible matches between `glycan` and `motif`.
  # Each match is an integer vector, with the same length as `vcount(motif)`.
  # The `i`-th element of the vector is the index of the vertex in `glycan`
  # that matches the `i`-th vertex in `motif`.
  # e.g. If `res[[1]] = c(2, 3)`,
  # it means the 1st vertex in `motif` matches the 2nd vertex in `glycan`,
  # and the 2nd vertex in `motif` matches the 3rd vertex in `glycan`.
  graph_subisomorphisms_vf2(
    graph = glycan,
    pattern = motif,
    graph_colors = glycan_colors,
    pattern_colors = motif_colors
  )
}


unique_vf2_res <- function(res) {
  if (length(res) <= 1) {
    return(res)
  }
  keys <- vapply(
    res,
    function(x) paste(sort(x), collapse = ","),
    character(1L),
    USE.NAMES = FALSE
  )
  res[!duplicated(keys)]
}


is_valid_result <- function(
  r,
  glycan,
  motif,
  alignment,
  ignore_linkages,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  context = NULL
) {
  # Optimized early exit using most selective checks first
  # Alignment check is often the most selective and fastest
  if (is.null(match_degree)) {
    if (
      !alignment_check(
        r,
        glycan,
        motif,
        context = context,
        alignment = alignment
      )
    ) {
      return(FALSE)
    }
  }

  # Residue check combines monosaccharide and substituent validation because
  # built-in modifications such as NAc and 5Ac may be represented either in
  # the mono name or in the sub attribute.
  if (
    !residue_check(
      r,
      glycan,
      motif,
      context = context,
      strict_sub = strict_sub,
      mode = mode
    )
  ) {
    return(FALSE)
  }

  if (
    !degree_check(
      r,
      glycan,
      motif,
      context = context,
      match_degree = match_degree
    )
  ) {
    return(FALSE)
  }

  # Only check linkages and anomer if linkages are not ignored
  # These are the most expensive checks, so do them last
  if (!ignore_linkages) {
    if (!linkage_check(r, glycan, motif, context = context, mode = mode)) {
      return(FALSE)
    }
    if (!anomer_check(r, glycan, motif, context = context, mode = mode)) {
      return(FALSE)
    }
  }

  return(TRUE)
}


#' Validate Residue Matching
#'
#' @param r A vector mapping motif vertices to glycan vertices.
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param strict_sub Whether substituent matching should be strict.
#' @param context A validation context.
#'
#' @return `TRUE` if all matched residues are compatible.
#' @noRd
residue_check <- function(
  r,
  glycan = NULL,
  motif = NULL,
  strict_sub,
  mode = "strict",
  context = NULL
) {
  if (!is.null(context)) {
    glycan_monos <- context$glycan_mono[r]
    motif_monos <- context$motif_mono
    glycan_subs <- context$glycan_sub[r]
    motif_subs <- context$motif_sub
  } else {
    glycan_monos <- graph_vertex_attr(glycan, "mono", r)
    motif_monos <- graph_vertex_attr(motif, "mono")
    glycan_subs <- graph_vertex_attr(glycan, "sub", r)
    motif_subs <- graph_vertex_attr(motif, "sub")
  }

  purrr::every(seq_along(glycan_monos), function(i) {
    match_residue(
      glycan_monos[[i]],
      glycan_subs[[i]],
      motif_monos[[i]],
      motif_subs[[i]],
      strict_sub,
      mode
    )
  })
}

#' Resolve Linkage Matching Mode
#'
#' @param glycan A glycan graph.
#' @param motif_has_linkages Whether the motif has informative linkages.
#' @param ignore_linkages Whether linkage matching should be ignored.
#'
#' @return One of three linkage modes:
#'   * `"ignore"`: skip `linkage_check()` and `anomer_check()` because the
#'     user requested `ignore_linkages = TRUE`, or because the motif has no
#'     informative linkage constraints to enforce.
#'   * `"none"`: return no matches before VF2 because the motif has informative
#'     linkage constraints but the glycan has no informative linkages to match.
#'   * `"check"`: run the normal VF2 candidate validation with linkage and
#'     anomer checks enabled.
#' @noRd
resolve_linkage_match_mode <- function(
  glycan,
  motif_has_linkages,
  ignore_linkages,
  mode = "strict",
  glycan_batch_profile = NULL
) {
  # Unlinked motifs are linkage-agnostic: once mono/sub/alignment checks pass,
  # wildcard linkage checks cannot reject additional candidates.
  if (ignore_linkages || !motif_has_linkages) {
    return("ignore")
  }

  # A linked motif cannot match an unlinked glycan. Returning "none" lets callers
  # bypass VF2 entirely and return the empty result for their output type.
  glycan_has_linkages <- if (is.null(glycan_batch_profile)) {
    graph_has_linkages(glycan)
  } else {
    glycan_batch_profile$has_linkages
  }
  if (mode == "strict" && !glycan_has_linkages) {
    return("none")
  }

  "check"
}

alignment_check <- function(
  r,
  glycan = NULL,
  motif = NULL,
  alignment,
  context = NULL
) {
  switch(
    alignment,
    "substructure" = TRUE,
    "core" = {
      if (is.null(context)) {
        glycan_core <- core_node(glycan)
        motif_core <- core_node(motif)
        r[[motif_core]] == glycan_core
      } else {
        r[[context$motif_core]] == context$glycan_core
      }
    },
    "terminal" = {
      if (is.null(context)) {
        glycan_terminals <- terminal_nodes(glycan)
        motif_terminals <- terminal_nodes(motif)
        all(r[motif_terminals] %in% glycan_terminals)
      } else {
        all(r[context$motif_terminals] %in% context$glycan_terminals)
      }
    },
    "whole" = {
      if (is.null(context)) {
        glycan_v <- igraph::vcount(glycan)
        motif_v <- igraph::vcount(motif)
        glycan_e <- igraph::ecount(glycan)
        motif_e <- igraph::ecount(motif)
        motif_v == glycan_v &&
          length(unique(r)) == glycan_v &&
          motif_e == glycan_e
      } else {
        context$motif_v == context$glycan_v &&
          length(unique(r)) == context$glycan_v &&
          context$motif_e == context$glycan_e
      }
    }
  )
}


#' Validate Degree Matching
#'
#' @param r A vector mapping motif vertices to glycan vertices.
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param match_degree A logical vector indicating which motif nodes to enforce degree matching.
#'
#' @return `TRUE` if all selected nodes have matching in- and out-degrees.
#' @noRd
degree_check <- function(
  r,
  glycan = NULL,
  motif = NULL,
  match_degree,
  context = NULL
) {
  if (is.null(match_degree) || !any(match_degree)) {
    return(TRUE)
  }

  if (!is.null(context)) {
    motif_in <- context$motif_in
    motif_out <- context$motif_out
    glycan_in <- context$glycan_in
    glycan_out <- context$glycan_out
  } else {
    motif_in <- graph_degree(motif, mode = "in")
    motif_out <- graph_degree(motif, mode = "out")
    glycan_in <- graph_degree(glycan, mode = "in")
    glycan_out <- graph_degree(glycan, mode = "out")
  }

  motif_idx <- which(match_degree)
  all(purrr::map_lgl(motif_idx, function(i) {
    g_idx <- r[[i]]
    motif_in[[i]] == glycan_in[[g_idx]] && motif_out[[i]] == glycan_out[[g_idx]]
  }))
}


match_sub <- function(glycan_sub, motif_sub, strict_sub, mode = "strict") {
  # Handle unstrict matching:
  if (!strict_sub && motif_sub == "") {
    return(TRUE)
  }

  # Handle empty substituents - both must be empty to match
  if (motif_sub == "" || glycan_sub == "") {
    return(motif_sub == "" && glycan_sub == "")
  }

  # Split substituents by comma
  glycan_subs <- stringr::str_split(glycan_sub, ",")[[1]]
  motif_subs <- stringr::str_split(motif_sub, ",")[[1]]

  # Remove empty strings
  glycan_subs <- glycan_subs[glycan_subs != ""]
  motif_subs <- motif_subs[motif_subs != ""]

  has_one_to_one_sub_match(glycan_subs, motif_subs, mode = mode)
}

#' Check Whether Substituents Have a One-to-One Assignment
#'
#' @param glycan_subs Glycan substituent tokens.
#' @param motif_subs Motif substituent tokens.
#' @param mode Matching mode.
#'
#' @return A logical scalar.
#' @noRd
has_one_to_one_sub_match <- function(glycan_subs, motif_subs, mode = "strict") {
  if (length(glycan_subs) != length(motif_subs)) {
    return(FALSE)
  }

  if (length(glycan_subs) == 0L) {
    return(TRUE)
  }

  candidates <- purrr::map(glycan_subs, function(g_sub) {
    which(purrr::map_lgl(motif_subs, function(m_sub) {
      match_single_sub(g_sub, m_sub, mode = mode)
    }))
  })

  if (any(lengths(candidates) == 0L)) {
    return(FALSE)
  }

  candidates <- candidates[order(lengths(candidates))]

  assign_next_sub <- function(i, used) {
    if (i > length(candidates)) {
      return(TRUE)
    }

    available <- setdiff(candidates[[i]], used)
    any(purrr::map_lgl(available, function(j) {
      assign_next_sub(i + 1L, c(used, j))
    }))
  }

  assign_next_sub(1L, integer())
}


#' Match a Glycan Residue to a Motif Residue
#'
#' @param glycan_mono Glycan monosaccharide name.
#' @param glycan_sub Glycan substituent string.
#' @param motif_mono Motif monosaccharide name.
#' @param motif_sub Motif substituent string.
#' @param strict_sub Whether substituent matching should be strict.
#'
#' @return A logical scalar.
#' @noRd
match_residue <- function(
  glycan_mono,
  glycan_sub,
  motif_mono,
  motif_sub,
  strict_sub,
  mode = "strict"
) {
  if (
    match_mono(glycan_mono, motif_mono, mode) &&
      match_sub(glycan_sub, motif_sub, strict_sub, mode)
  ) {
    return(TRUE)
  }

  if (!is_fuzzy_sub(motif_sub)) {
    return(FALSE)
  }

  built_in <- decompose_builtin_modification(glycan_mono)
  if (is.null(built_in) || !match_mono(built_in$mono, motif_mono, mode)) {
    return(FALSE)
  }

  match_sub(combine_subs(built_in$sub, glycan_sub), motif_sub, strict_sub, mode)
}

#' Match a Glycan Monosaccharide to a Motif Monosaccharide
#'
#' @param glycan_mono Glycan monosaccharide name.
#' @param motif_mono Motif monosaccharide name.
#' @param mode Matching mode.
#'
#' @return A logical scalar.
#' @noRd
match_mono <- function(glycan_mono, motif_mono, mode = "strict") {
  if (glycan_mono == motif_mono) {
    return(TRUE)
  }

  if (mode != "lenient") {
    return(FALSE)
  }

  glycan_type <- glyrepr::get_mono_type(glycan_mono)
  motif_type <- glyrepr::get_mono_type(motif_mono)

  if (glycan_type == "generic" && motif_type == "concrete") {
    return(glycan_mono == glyrepr::convert_to_generic(motif_mono))
  }

  if (glycan_type == "concrete" && motif_type == "generic") {
    return(glyrepr::convert_to_generic(glycan_mono) == motif_mono)
  }

  FALSE
}


#' Test Whether a Motif Substituent Contains a Fuzzy Position
#'
#' @param sub A substituent string.
#'
#' @return A logical scalar.
#' @noRd
is_fuzzy_sub <- function(sub) {
  isTRUE(stringr::str_detect(sub, "(^|,)\\?"))
}


#' Decompose Built-In Residue Modifications
#'
#' @param mono A monosaccharide name.
#'
#' @return `NULL` or a list with `mono` and `sub` fields.
#' @noRd
decompose_builtin_modification <- function(mono) {
  patterns <- c(
    NAc = "NAc",
    N = "N"
  )

  if (mono == "Neu5Ac") {
    return(list(mono = "Neu", sub = "5Ac"))
  }

  for (sub in names(patterns)) {
    if (stringr::str_ends(mono, stringr::fixed(patterns[[sub]]))) {
      base <- stringr::str_remove(mono, stringr::fixed(patterns[[sub]]))
      if (base != "" && is_known_residue_base(base)) {
        return(list(mono = base, sub = stringr::str_c("?", sub)))
      }
    }
  }

  NULL
}


#' Check Whether a Decomposed Residue Base Is Known
#'
#' @param mono A monosaccharide name.
#'
#' @return A logical scalar.
#' @noRd
is_known_residue_base <- function(mono) {
  mono %in% glyrepr::available_monosaccharides()
}


#' Combine Implicit and Explicit Substituents
#'
#' @param implicit_sub A substituent from built-in residue decomposition.
#' @param explicit_sub A substituent already stored on the glycan graph.
#'
#' @return A comma-separated substituent string.
#' @noRd
combine_subs <- function(implicit_sub, explicit_sub) {
  subs <- c(implicit_sub, explicit_sub)
  subs <- subs[!is.na(subs) & subs != ""]
  stringr::str_c(subs, collapse = ",")
}

# Helper function to match a single substituent (handles obscure linkages)
match_single_sub <- function(glycan_sub, motif_sub, mode = "strict") {
  # Extract position and substituent parts
  # Format: "3Me", "6S", "?Me", "?S", etc.

  # Extract position (first character) and substituent (rest)
  motif_pos <- stringr::str_sub(motif_sub, 1, 1)
  motif_rest <- stringr::str_sub(motif_sub, 2)

  glycan_pos <- stringr::str_sub(glycan_sub, 1, 1)
  glycan_rest <- stringr::str_sub(glycan_sub, 2)

  pos_match <- function(motif_pos, glycan_pos) {
    motif_pos == "?" ||
      (mode == "lenient" && glycan_pos == "?") ||
      motif_pos == glycan_pos
  }
  sub_match <- function(motif_rest, glycan_rest) {
    motif_rest == glycan_rest
  }

  pos_match(motif_pos, glycan_pos) && sub_match(motif_rest, glycan_rest)
}


linkage_check <- function(
  r,
  glycan = NULL,
  motif = NULL,
  context = NULL,
  mode = "strict"
) {
  if (is.null(context)) {
    edges <- get_corresponding_edges(r, glycan, motif)
    glycan_linkages <- edges$glycan$linkage
    motif_linkages <- edges$motif$linkage
  } else {
    if (nrow(context$motif_edge_list) == 0) {
      return(TRUE)
    }

    glycan_edge_list <- matrix(
      r[as.vector(context$motif_edge_list)],
      ncol = 2L
    )
    glycan_linkages <- context$glycan_edge_linkages[glycan_edge_list]
    motif_linkages <- context$motif_linkages
  }
  for (i in seq_along(glycan_linkages)) {
    if (
      !match_linkage(glycan_linkages[[i]], motif_linkages[[i]], mode = mode)
    ) {
      return(FALSE)
    }
  }
  return(TRUE)
}


get_corresponding_edges <- function(r, glycan, motif) {
  motif_edge_list <- igraph::as_edgelist(motif, names = FALSE)

  glycan_edge_ids <- purrr::map_int(
    seq_len(nrow(motif_edge_list)),
    function(i) {
      motif_edge <- motif_edge_list[i, ] # c(node_id_1, node_id_2)
      glycan_edge <- r[motif_edge]
      # Convert igraph.vs to vector as get_edge_ids expects a simple vector
      # This fix a bug introduced by igraph v2.2.0
      igraph::get_edge_ids(glycan, as.vector(glycan_edge))
    }
  )

  list(
    glycan = list(
      linkage = graph_edge_attr(glycan, "linkage", glycan_edge_ids)
    ),
    motif = list(linkage = graph_edge_attr(motif, "linkage"))
  )
}


match_linkage <- function(glycan_linkage, motif_linkage, mode = "strict") {
  gl <- parse_linkage(glycan_linkage)
  ml <- parse_linkage(motif_linkage)

  anomer_ok <- function(gl, ml) {
    ml[["anomer"]] == "?" ||
      (mode == "lenient" && gl[["anomer"]] == "?") ||
      ml[["anomer"]] == gl[["anomer"]]
  }
  pos1_ok <- function(gl, ml) {
    ml[["pos1"]] == "?" ||
      (mode == "lenient" && gl[["pos1"]] == "?") ||
      ml[["pos1"]] == gl[["pos1"]]
  }
  pos2_ok <- function(gl, ml) {
    glycan_pos2 <- parse_pos2(gl[["pos2"]])
    motif_pos2 <- parse_pos2(ml[["pos2"]])

    if (any(motif_pos2 == "?")) {
      return(TRUE)
    }
    if (mode == "lenient" && any(glycan_pos2 == "?")) {
      return(TRUE)
    }
    if (mode == "lenient") {
      return(any(glycan_pos2 %in% motif_pos2))
    }

    all(glycan_pos2 %in% motif_pos2)
  }

  anomer_ok(gl, ml) && pos1_ok(gl, ml) && pos2_ok(gl, ml)
}


parse_linkage <- function(linkage) {
  c(
    anomer = stringr::str_sub(linkage, 1, 1),
    pos1 = stringr::str_sub(linkage, 2, 2),
    pos2 = stringr::str_sub(linkage, 4, -1)
  )
}


parse_pos2 <- function(pos2) {
  if (stringr::str_detect(pos2, "/")) {
    stringr::str_split(pos2, "/")[[1]]
  } else {
    pos2
  }
}


anomer_check <- function(
  r,
  glycan = NULL,
  motif = NULL,
  context = NULL,
  mode = "strict"
) {
  if (is.null(context)) {
    glycan_core <- core_node(glycan)
    motif_core <- core_node(motif)
    matched_g_node <- r[[motif_core]]

    if (matched_g_node == glycan_core) {
      # This means two cores are matched.
      match_anomer(glycan$anomer, motif$anomer, mode = mode)
    } else {
      # The motif anomer should match a linkage in the glycan.
      in_edge <- graph_incident_edge_ids(glycan, matched_g_node, mode = "in")
      linkage <- graph_edge_attr(glycan, "linkage", in_edge)
      linkage_anomer <- stringr::str_split_1(linkage, "-")[[1]]
      match_anomer(linkage_anomer, motif$anomer, mode = mode)
    }
  } else {
    matched_g_node <- r[[context$motif_core]]

    if (matched_g_node == context$glycan_core) {
      # This means two cores are matched.
      match_anomer(context$glycan_anomer, context$motif_anomer, mode = mode)
    } else {
      # The motif anomer should match a linkage in the glycan.
      linkage <- context$glycan_incoming_linkages[[matched_g_node]]
      linkage_anomer <- stringr::str_split_1(linkage, "-")[[1]]
      match_anomer(linkage_anomer, context$motif_anomer, mode = mode)
    }
  }
}


match_anomer <- function(glycan_anomer, motif_anomer, mode = "strict") {
  # Check if the anomer of the glycan and motif are matched.
  # - "??" in motif will match any anomer and position in glycan.
  # - "?1" in motif will match any anomer at position 1 in glycan.
  # - "a?" in motif will match anomer "a" at any position in glycan.
  # - "a1" in motif will only match anomer "a" at position 1 in glycan.
  ga <- parse_anomer(glycan_anomer)
  ma <- parse_anomer(motif_anomer)

  anomer_ok <- function(ga, ma) {
    ma[["anomer"]] == "?" ||
      (mode == "lenient" && ga[["anomer"]] == "?") ||
      ma[["anomer"]] == ga[["anomer"]]
  }
  position_ok <- function(ga, ma) {
    ma[["pos"]] == "?" ||
      (mode == "lenient" && ga[["pos"]] == "?") ||
      ma[["pos"]] == ga[["pos"]]
  }

  anomer_ok(ga, ma) && position_ok(ga, ma)
}


parse_anomer <- function(anomer) {
  c(
    anomer = stringr::str_sub(anomer, 1, 1),
    pos = stringr::str_sub(anomer, 2)
  )
}


terminal_nodes <- function(glycan) {
  out_degree <- graph_degree(glycan, mode = "out")
  which(out_degree == 0)
}


core_node <- function(glycan) {
  # This is a hack based on one important property of glyrepr's glycan structure vectors:
  # the order of the nodes in the graph is consistent with the order they appear
  # in the IUPAC-condensed sequences.
  # Therefore, the last node in the graph is guaranteed to be the root node.
  igraph::vcount(glycan)
}
