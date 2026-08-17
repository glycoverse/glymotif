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

  motif_type <- graph_mono_type(motif)
  if (motif_type == "mixed") {
    return("none")
  }

  if (fuzzy) {
    if (motif_type == "generic") {
      return("none")
    }
    return("base")
  }

  if (motif_type == "generic") {
    return("generic")
  }
  "exact"
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
    base = residue_base_keys(monos),
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


#' Create Base Residue Keys for Fuzzy Modification Matching
#'
#' @param monos A character vector of monosaccharide names.
#'
#' @return A character vector of base residue keys.
#' @noRd
residue_base_keys <- function(monos) {
  purrr::map_chr(monos, function(mono) {
    built_in <- decompose_builtin_modification(mono)
    if (is.null(built_in)) {
      mono
    } else {
      built_in$mono
    }
  })
}


#' Create a Base VF2 Compatibility Context
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#'
#' @return A base compatibility context.
#' @noRd
new_vf2_context <- function(glycan, motif) {
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
    glycan_anomer = NULL,
    motif_anomer = NULL,
    glycan_incoming_linkages = NULL
  )
}


#' Add Alignment Data to a VF2 Compatibility Context
#'
#' @param context A VF2 compatibility context.
#' @param alignment Alignment mode.
#'
#' @return A compatibility context with alignment-specific data.
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
    "whole" = context
  )
}


#' Add Degree Data to a VF2 Compatibility Context
#'
#' @param context A VF2 compatibility context.
#'
#' @return A compatibility context with degree vectors.
#' @noRd
add_degree_context <- function(context) {
  context$motif_in <- graph_degree(context$motif, mode = "in")
  context$motif_out <- graph_degree(context$motif, mode = "out")
  context$glycan_in <- graph_degree(context$glycan, mode = "in")
  context$glycan_out <- graph_degree(context$glycan, mode = "out")
  context
}


#' Add Core Node Data to a VF2 Compatibility Context
#'
#' @param context A VF2 compatibility context.
#'
#' @return A compatibility context with core node indices.
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


#' Add Linkage and Anomer Data to a VF2 Compatibility Context
#'
#' @param context A VF2 compatibility context.
#'
#' @return A compatibility context with linkage and anomer data.
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

  context
}


#' Prepare Reusable Graph Data for VF2 Compatibility
#'
#' @param glycan A glycan graph.
#' @param motif A motif graph.
#' @param alignment Alignment mode.
#' @param ignore_linkages Whether linkage/anomer checks are skipped.
#' @param match_degree Degree matching vector.
#'
#' @return A context containing only data needed for compatibility predicates.
#' @noRd
prepare_vf2_context <- function(
  glycan,
  motif,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL
) {
  context <- new_vf2_context(glycan, motif)
  context <- add_alignment_context(context, alignment)
  if (!is.null(match_degree) && any(match_degree)) {
    context <- add_degree_context(context)
  }
  if (!ignore_linkages) {
    context <- add_linkage_context(context)
  }
  context
}


#' Resolve Linkage Matching Mode
#'
#' @param glycan A glycan graph.
#' @param motif_has_linkages Whether the motif has informative linkages.
#' @param ignore_linkages Whether linkage matching should be ignored.
#'
#' @return One of three linkage modes:
#'   * `"ignore"`: skip linkage and anomer compatibility because the user
#'     requested `ignore_linkages = TRUE`, or because the motif has no
#'     informative constraints to enforce.
#'   * `"none"`: return no matches before VF2 because the motif has informative
#'     linkage constraints but the glycan has no informative linkages to match.
#'   * `"check"`: enable linkage and anomer compatibility predicates in VF2.
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

  glycan_type <- glyrepr::get_mono_type(glycan_mono)
  motif_type <- glyrepr::get_mono_type(motif_mono)

  if (glycan_type == "concrete" && motif_type == "generic") {
    return(glyrepr::convert_to_generic(glycan_mono) == motif_mono)
  }

  if (
    mode == "lenient" &&
      glycan_type == "generic" &&
      motif_type == "concrete"
  ) {
    return(glycan_mono == glyrepr::convert_to_generic(motif_mono))
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


core_node <- function(glycan) {
  # This is a hack based on one important property of glyrepr's glycan structure vectors:
  # the order of the nodes in the graph is consistent with the order they appear
  # in the IUPAC-condensed sequences.
  # Therefore, the last node in the graph is guaranteed to be the root node.
  igraph::vcount(glycan)
}
