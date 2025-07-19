colorize_graphs <- function(glycan, motif) {
  # Add "color" vertex attributes to the graph.
  # The colors are converted from the "mono" vertex attributes.
  unique_monos <- unique(c(igraph::V(glycan)$mono, igraph::V(motif)$mono))
  color_map <- seq_along(unique_monos)
  names(color_map) <- unique_monos
  glycan_colors <- color_map[igraph::V(glycan)$mono]
  motif_colors <- color_map[igraph::V(motif)$mono]
  names(glycan_colors) <- NULL
  names(motif_colors) <- NULL
  igraph::V(glycan)$color <- glycan_colors
  igraph::V(motif)$color <- motif_colors
  list(glycan = glycan, motif = motif)
}


perform_vf2 <- function(glycan, motif) {
  # Perform "VF2" algorithm
  # `res` is a list of all possible matches between `glycan` and `motif`.
  # Each match is an integer vector, with the same length as `vcount(motif)`.
  # The `i`-th element of the vector is the index of the vertex in `glycan`
  # that matches the `i`-th vertex in `motif`.
  # e.g. If `res[[1]] = c(2, 3)`,
  # it means the 1st vertex in `motif` matches the 2nd vertex in `glycan`,
  # and the 2nd vertex in `motif` matches the 3rd vertex in `glycan`.
  igraph::graph.get.subisomorphisms.vf2(glycan, motif)
}


is_valid_result <- function(r, glycan, motif, alignment, ignore_linkages) {
  # Optimized early exit using most selective checks first
  # Alignment check is often the most selective and fastest
  if (!alignment_check(r, glycan, motif, alignment = alignment)) {
    return(FALSE)
  }
  
  # Substituent check is relatively fast and often selective
  if (!substituent_check(r, glycan, motif, alignment = alignment)) {
    return(FALSE)
  }
  
  # Only check linkages and anomer if linkages are not ignored
  # These are the most expensive checks, so do them last
  if (!ignore_linkages) {
    if (!linkage_check(r, glycan, motif, alignment = alignment)) {
      return(FALSE)
    }
    if (!anomer_check(r, glycan, motif, alignment = alignment)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}


alignment_check <- function(r, glycan, motif, alignment) {
  # Fast path for most common case
  if (alignment == "substructure") return(TRUE)

  # Core alignment check (most common for N-glycan motifs)
  if (alignment == "core") {
    glycan_core <- core_node(glycan)
    motif_core <- core_node(motif)
    return(r[[motif_core]] == glycan_core)
  }
  
  # Terminal alignment check
  if (alignment == "terminal") {
    glycan_terminals <- terminal_nodes(glycan)
    motif_terminals <- terminal_nodes(motif)
    return(all(r[motif_terminals] %in% glycan_terminals))
  }

  # Whole glycan alignment check (most expensive, do last)
  if (alignment == "whole") {
    return(igraph::isomorphic(glycan, motif, method = "vf2"))
  }
  
  # Default return for unknown alignment types
  return(FALSE)
}


substituent_check <- function(r, glycan, motif, ...) {
  glycan_subs <- igraph::V(glycan)$sub[r]
  motif_subs <- igraph::V(motif)$sub
  all(purrr::map2_lgl(glycan_subs, motif_subs, match_sub))
}


match_sub <- function(glycan_sub, motif_sub) {
  # Handle empty substituents
  if (motif_sub == "" && glycan_sub == "") {
    return(TRUE)
  }
  if (motif_sub == "" && glycan_sub != "") {
    return(FALSE)
  }
  if (motif_sub != "" && glycan_sub == "") {
    return(FALSE)
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

# Helper function to match a single substituent (handles obscure linkages)
match_single_sub <- function(glycan_sub, motif_sub) {
  # Extract position and substituent parts
  # Format: "3Me", "6S", "?Me", "?S", etc.

  # Extract position (first character) and substituent (rest)
  motif_pos <- stringr::str_sub(motif_sub, 1, 1)
  motif_rest <- stringr::str_sub(motif_sub, 2)

  glycan_pos <- stringr::str_sub(glycan_sub, 1, 1)
  glycan_rest <- stringr::str_sub(glycan_sub, 2)

  # Check if positions match (? in motif matches any position)
  pos_match <- (motif_pos == "?" || motif_pos == glycan_pos)

  # Check if substituent parts match
  sub_match <- (motif_rest == glycan_rest)

  pos_match && sub_match
}


linkage_check <- function(r, glycan, motif, ...) {
  edges <- get_corresponding_edges(r, glycan, motif)
  glycan_linkages <- edges$glycan$linkage
  motif_linkages <- edges$motif$linkage
  all(purrr::map2_lgl(glycan_linkages, motif_linkages, match_linkage))
}


get_corresponding_edges <- function(r, glycan, motif) {
  motif_edge_list <- igraph::as_edgelist(motif, names = FALSE)

  glycan_edge_ids <- purrr::map_int(seq_len(nrow(motif_edge_list)), function(i) {
    motif_edge <- motif_edge_list[i, ]  # c(node_id_1, node_id_2)
    glycan_edge <- r[motif_edge]
    igraph::get_edge_ids(glycan, glycan_edge)
  })

  motif_edges <- igraph::E(motif)
  glycan_edges <- igraph::E(glycan)[glycan_edge_ids]
  list(glycan = glycan_edges, motif = motif_edges)
}


match_linkage <- function(glycan_linkage, motif_linkage) {
  gl <- parse_linkage(glycan_linkage)
  ml <- parse_linkage(motif_linkage)

  anomer_ok <- (ml[["anomer"]] == "?" || ml[["anomer"]] == gl[["anomer"]])
  pos1_ok <- (ml[["pos1"]] == "?" || ml[["pos1"]] == gl[["pos1"]])
  pos2_ok <- (ml[["pos2"]] == "?" || all(parse_pos2(gl[["pos2"]]) %in% parse_pos2(ml[["pos2"]])))

  anomer_ok && pos1_ok && pos2_ok
}


parse_linkage <- function(linkage) {
  c(
    anomer = stringr::str_sub(linkage, 1, 1),
    pos1   = stringr::str_sub(linkage, 2, 2),
    pos2   = stringr::str_sub(linkage, 4, -1)
  )
}


parse_pos2 <- function(pos2) {
  if (stringr::str_detect(pos2, "/")) {
    stringr::str_split(pos2, "/")[[1]]
  } else {
    pos2
  }
}


anomer_check <- function(r, glycan, motif, ...) {
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
}


match_anomer <- function(glycan_anomer, motif_anomer) {
  # Check if the anomer of the glycan and motif are matched.
  # - "??" in motif will match any anomer and position in glycan.
  # - "?1" in motif will match any anomer at position 1 in glycan.
  # - "a?" in motif will match anomer "a" at any position in glycan.
  # - "a1" in motif will only match anomer "a" at position 1 in glycan.
  ga <- parse_anomer(glycan_anomer)
  ma <- parse_anomer(motif_anomer)

  anomer_ok <- ma[["anomer"]] == "?" || ma[["anomer"]] == ga[["anomer"]]
  position_ok <- ma[["pos"]] == "?" || ma[["pos"]] == ga[["pos"]]

  anomer_ok && position_ok
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
  in_degree <- igraph::degree(glycan, mode = "in")
  igraph::V(glycan)[in_degree == 0]
}
