colorize_graphs <- function(glycan, motif) {
  # Prepare VF2 color vectors from the "mono" vertex attributes without
  # mutating the input graphs.
  glycan_monos <- igraph::V(glycan)$mono
  motif_monos <- igraph::V(motif)$mono
  if (has_fuzzy_modification(motif)) {
    glycan_monos <- residue_color_keys(glycan_monos)
    motif_monos <- residue_color_keys(motif_monos)
  }

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
#'
#' @return A list containing residue keys, required counts, and keying mode.
#' @noRd
new_motif_composition_profile <- function(motif) {
  use_base_keys <- has_fuzzy_modification(motif)
  motif_monos <- igraph::vertex_attr(motif, "mono")
  if (use_base_keys) {
    motif_monos <- residue_color_keys(motif_monos)
  }

  keys <- unique(motif_monos)
  counts <- tabulate(match(motif_monos, keys), nbins = length(keys))
  list(
    keys = keys,
    counts = counts,
    use_base_keys = use_base_keys
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
composition_can_match <- function(glycan, motif_profile) {
  glycan_monos <- igraph::vertex_attr(glycan, "mono")
  if (motif_profile$use_base_keys) {
    glycan_monos <- residue_color_keys(glycan_monos)
  }

  glycan_counts <- tabulate(
    match(glycan_monos, motif_profile$keys),
    nbins = length(motif_profile$keys)
  )
  all(glycan_counts >= motif_profile$counts)
}


#' Check Whether Whole-Alignment Graph Sizes Can Match
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param alignment Alignment mode.
#'
#' @return A logical scalar.
#' @noRd
whole_alignment_size_can_match <- function(glycan, motif, alignment) {
  if (alignment != "whole") {
    return(TRUE)
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
  strict_sub
) {
  if (alignment != "core") {
    return(TRUE)
  }

  glycan_core <- core_node(glycan)
  motif_core <- core_node(motif)

  match_residue(
    igraph::vertex_attr(glycan, "mono", index = glycan_core),
    igraph::vertex_attr(glycan, "sub", index = glycan_core),
    igraph::vertex_attr(motif, "mono", index = motif_core),
    igraph::vertex_attr(motif, "sub", index = motif_core),
    strict_sub = strict_sub
  )
}


#' Check Whether a Motif Has Fuzzy Built-In Modifications
#'
#' @param motif A motif graph.
#'
#' @return A logical scalar.
#' @noRd
has_fuzzy_modification <- function(motif) {
  any(purrr::map_lgl(igraph::V(motif)$sub, is_fuzzy_sub))
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
    glycan_mono = igraph::vertex_attr(glycan, "mono"),
    motif_mono = igraph::vertex_attr(motif, "mono"),
    glycan_sub = igraph::vertex_attr(glycan, "sub"),
    motif_sub = igraph::vertex_attr(motif, "sub"),
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
        igraph::degree(context$glycan, mode = "out") == 0
      )
      context$motif_terminals <- which(
        igraph::degree(context$motif, mode = "out") == 0
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
  context$motif_in <- igraph::degree(context$motif, mode = "in")
  context$motif_out <- igraph::degree(context$motif, mode = "out")
  context$glycan_in <- igraph::degree(context$glycan, mode = "in")
  context$glycan_out <- igraph::degree(context$glycan, mode = "out")
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
      igraph::degree(context$glycan, mode = "in") == 0
    )
    context$motif_core <- which(igraph::degree(context$motif, mode = "in") == 0)
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
  glycan_linkages <- igraph::edge_attr(context$glycan, "linkage")
  context$glycan_incoming_linkages <- rep(
    NA_character_,
    igraph::vcount(context$glycan)
  )
  if (nrow(glycan_edge_list) > 0) {
    context$glycan_incoming_linkages[glycan_edge_list[, 2]] <- glycan_linkages
  }

  motif_edge_list <- igraph::as_edgelist(context$motif, names = FALSE)
  context$motif_edge_list <- motif_edge_list
  context$motif_linkages <- igraph::edge_attr(context$motif, "linkage")
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
  if (is.null(glycan_colors) && is.null(motif_colors)) {
    igraph::graph.get.subisomorphisms.vf2(glycan, motif)
  } else {
    igraph::graph.get.subisomorphisms.vf2(
      glycan,
      motif,
      vertex.color1 = glycan_colors,
      vertex.color2 = motif_colors
    )
  }
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
      strict_sub = strict_sub
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
    if (!linkage_check(r, glycan, motif, context = context)) {
      return(FALSE)
    }
    if (!anomer_check(r, glycan, motif, context = context)) {
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
  context = NULL
) {
  if (!is.null(context)) {
    glycan_monos <- context$glycan_mono[r]
    motif_monos <- context$motif_mono
    glycan_subs <- context$glycan_sub[r]
    motif_subs <- context$motif_sub
  } else {
    glycan_monos <- igraph::V(glycan)$mono[r]
    motif_monos <- igraph::V(motif)$mono
    glycan_subs <- igraph::V(glycan)$sub[r]
    motif_subs <- igraph::V(motif)$sub
  }

  purrr::every(seq_along(glycan_monos), function(i) {
    match_residue(
      glycan_monos[[i]],
      glycan_subs[[i]],
      motif_monos[[i]],
      motif_subs[[i]],
      strict_sub
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
  ignore_linkages
) {
  # Unlinked motifs are linkage-agnostic: once mono/sub/alignment checks pass,
  # wildcard linkage checks cannot reject additional candidates.
  if (ignore_linkages || !motif_has_linkages) {
    return("ignore")
  }

  # A linked motif cannot match an unlinked glycan. Returning "none" lets callers
  # bypass VF2 entirely and return the empty result for their output type.
  if (!graph_has_linkages(glycan)) {
    return("none")
  }

  "check"
}

#' Determine Whether a Graph Has Informative Linkages
#'
#' @param glycan A glycan graph.
#'
#' @return A logical scalar.
#' @noRd
graph_has_linkages <- function(glycan) {
  any(igraph::edge_attr(glycan, "linkage") != "??-?") ||
    igraph::graph_attr(glycan, "anomer") != "??"
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
    motif_in <- igraph::degree(motif, mode = "in")
    motif_out <- igraph::degree(motif, mode = "out")
    glycan_in <- igraph::degree(glycan, mode = "in")
    glycan_out <- igraph::degree(glycan, mode = "out")
  }

  motif_idx <- which(match_degree)
  all(purrr::map_lgl(motif_idx, function(i) {
    g_idx <- r[[i]]
    motif_in[[i]] == glycan_in[[g_idx]] && motif_out[[i]] == glycan_out[[g_idx]]
  }))
}


match_sub <- function(glycan_sub, motif_sub, strict_sub) {
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

  # For strict matching: every motif substituent must match a glycan substituent
  # and every glycan substituent must be matched by some motif substituent

  # Check if all motif substituents are matched
  motif_matched <- purrr::map_lgl(motif_subs, function(m_sub) {
    any(purrr::map_lgl(glycan_subs, function(g_sub) {
      match_single_sub(g_sub, m_sub)
    }))
  })

  # Check if all glycan substituents are matched (unless motif has wildcards)
  glycan_matched <- purrr::map_lgl(glycan_subs, function(g_sub) {
    any(purrr::map_lgl(motif_subs, function(m_sub) {
      match_single_sub(g_sub, m_sub)
    }))
  })

  all(motif_matched) && all(glycan_matched)
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
  strict_sub
) {
  if (
    glycan_mono == motif_mono &&
      match_sub(glycan_sub, motif_sub, strict_sub)
  ) {
    return(TRUE)
  }

  if (!is_fuzzy_sub(motif_sub)) {
    return(FALSE)
  }

  built_in <- decompose_builtin_modification(glycan_mono)
  if (is.null(built_in) || built_in$mono != motif_mono) {
    return(FALSE)
  }

  match_sub(combine_subs(built_in$sub, glycan_sub), motif_sub, strict_sub)
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
match_single_sub <- function(glycan_sub, motif_sub) {
  # Extract position and substituent parts
  # Format: "3Me", "6S", "?Me", "?S", etc.

  # Extract position (first character) and substituent (rest)
  motif_pos <- stringr::str_sub(motif_sub, 1, 1)
  motif_rest <- stringr::str_sub(motif_sub, 2)

  glycan_pos <- stringr::str_sub(glycan_sub, 1, 1)
  glycan_rest <- stringr::str_sub(glycan_sub, 2)

  pos_match <- function(motif_pos, glycan_pos) {
    motif_pos == "?" || motif_pos == glycan_pos
  }
  sub_match <- function(motif_rest, glycan_rest) {
    motif_rest == glycan_rest
  }

  pos_match(motif_pos, glycan_pos) && sub_match(motif_rest, glycan_rest)
}


linkage_check <- function(r, glycan = NULL, motif = NULL, context = NULL) {
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
    if (!match_linkage(glycan_linkages[[i]], motif_linkages[[i]])) {
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

  motif_edges <- igraph::E(motif)
  glycan_edges <- igraph::E(glycan)[glycan_edge_ids]
  list(glycan = glycan_edges, motif = motif_edges)
}


match_linkage <- function(glycan_linkage, motif_linkage) {
  gl <- parse_linkage(glycan_linkage)
  ml <- parse_linkage(motif_linkage)

  anomer_ok <- function(gl, ml) {
    ml[["anomer"]] == "?" || ml[["anomer"]] == gl[["anomer"]]
  }
  pos1_ok <- function(gl, ml) {
    ml[["pos1"]] == "?" || ml[["pos1"]] == gl[["pos1"]]
  }
  pos2_ok <- function(gl, ml) {
    ml[["pos2"]] == "?" ||
      all(parse_pos2(gl[["pos2"]]) %in% parse_pos2(ml[["pos2"]]))
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


anomer_check <- function(r, glycan = NULL, motif = NULL, context = NULL) {
  if (is.null(context)) {
    glycan_core <- core_node(glycan)
    motif_core <- core_node(motif)
    matched_g_node <- r[[motif_core]]

    if (matched_g_node == glycan_core) {
      # This means two cores are matched.
      match_anomer(glycan$anomer, motif$anomer)
    } else {
      # The motif anomer should match a linkage in the glycan.
      linkage <- igraph::incident(glycan, matched_g_node, mode = "in")$linkage
      linkage_anomer <- stringr::str_split_1(linkage, "-")[[1]]
      match_anomer(linkage_anomer, motif$anomer)
    }
  } else {
    matched_g_node <- r[[context$motif_core]]

    if (matched_g_node == context$glycan_core) {
      # This means two cores are matched.
      match_anomer(context$glycan_anomer, context$motif_anomer)
    } else {
      # The motif anomer should match a linkage in the glycan.
      linkage <- context$glycan_incoming_linkages[[matched_g_node]]
      linkage_anomer <- stringr::str_split_1(linkage, "-")[[1]]
      match_anomer(linkage_anomer, context$motif_anomer)
    }
  }
}


match_anomer <- function(glycan_anomer, motif_anomer) {
  # Check if the anomer of the glycan and motif are matched.
  # - "??" in motif will match any anomer and position in glycan.
  # - "?1" in motif will match any anomer at position 1 in glycan.
  # - "a?" in motif will match anomer "a" at any position in glycan.
  # - "a1" in motif will only match anomer "a" at position 1 in glycan.
  ga <- parse_anomer(glycan_anomer)
  ma <- parse_anomer(motif_anomer)

  anomer_ok <- function(ga, ma) {
    ma[["anomer"]] == "?" || ma[["anomer"]] == ga[["anomer"]]
  }
  position_ok <- function(ga, ma) {
    ma[["pos"]] == "?" || ma[["pos"]] == ga[["pos"]]
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
  out_degree <- igraph::degree(glycan, mode = "out")
  igraph::V(glycan)[out_degree == 0]
}


core_node <- function(glycan) {
  # This is a hack based on one important property of glyrepr's glycan structure vectors:
  # the order of the nodes in the graph is consistent with the order they appear
  # in the IUPAC-condensed sequences.
  # Therefore, the last node in the graph is guaranteed to be the root node.
  igraph::vcount(glycan)
}
