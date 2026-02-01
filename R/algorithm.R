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


unique_vf2_res <- function(res) {
  if (length(res) <= 1) return(res)
  keys <- vapply(
    res,
    function(x) paste(sort(x), collapse = ","),
    character(1L),
    USE.NAMES = FALSE
  )
  res[!duplicated(keys)]
}


is_valid_result <- function(r, glycan, motif, alignment, ignore_linkages, strict_sub = TRUE, match_degree = NULL) {
  # Optimized early exit using most selective checks first
  # Alignment check is often the most selective and fastest
  if (is.null(match_degree)) {
    if (!alignment_check(r, glycan, motif, alignment = alignment)) {
      return(FALSE)
    }
  }

  # Substituent check is relatively fast and often selective
  if (!substituent_check(r, glycan, motif, strict_sub = strict_sub)) {
    return(FALSE)
  }

  if (!degree_check(r, glycan, motif, match_degree)) {
    return(FALSE)
  }

  # Only check linkages and anomer if linkages are not ignored
  # These are the most expensive checks, so do them last
  if (!ignore_linkages) {
    if (!linkage_check(r, glycan, motif)) {
      return(FALSE)
    }
    if (!anomer_check(r, glycan, motif)) {
      return(FALSE)
    }
  }

  return(TRUE)
}


alignment_check <- function(r, glycan, motif, alignment) {
  switch(alignment,
    "substructure" = TRUE,
    "core" = {
      glycan_core <- core_node(glycan)
      motif_core <- core_node(motif)
      r[[motif_core]] == glycan_core
    },
    "terminal" = {
      glycan_terminals <- terminal_nodes(glycan)
      motif_terminals <- terminal_nodes(motif)
      all(r[motif_terminals] %in% glycan_terminals)
    },
    "whole" = {
      glycan_v <- igraph::vcount(glycan)
      motif_v <- igraph::vcount(motif)
      glycan_e <- igraph::ecount(glycan)
      motif_e <- igraph::ecount(motif)
      motif_v == glycan_v && length(unique(r)) == glycan_v && motif_e == glycan_e
    }
  )
}


substituent_check <- function(r, glycan, motif, strict_sub) {
  glycan_subs <- igraph::V(glycan)$sub[r]
  motif_subs <- igraph::V(motif)$sub
  for (i in seq_along(glycan_subs)) {
    if (!match_sub(glycan_subs[[i]], motif_subs[[i]], strict_sub)) {
      return(FALSE)
    }
  }
  return(TRUE)
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
degree_check <- function(r, glycan, motif, match_degree) {
  if (is.null(match_degree) || !any(match_degree)) {
    return(TRUE)
  }

  motif_in <- igraph::degree(motif, mode = "in")
  motif_out <- igraph::degree(motif, mode = "out")
  glycan_in <- igraph::degree(glycan, mode = "in")
  glycan_out <- igraph::degree(glycan, mode = "out")

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


linkage_check <- function(r, glycan, motif) {
  edges <- get_corresponding_edges(r, glycan, motif)
  glycan_linkages <- edges$glycan$linkage
  motif_linkages <- edges$motif$linkage
  for (i in seq_along(glycan_linkages)) {
    if (!match_linkage(glycan_linkages[[i]], motif_linkages[[i]])) {
      return(FALSE)
    }
  }
  return(TRUE)
}


get_corresponding_edges <- function(r, glycan, motif) {
  motif_edge_list <- igraph::as_edgelist(motif, names = FALSE)

  glycan_edge_ids <- purrr::map_int(seq_len(nrow(motif_edge_list)), function(i) {
    motif_edge <- motif_edge_list[i, ]  # c(node_id_1, node_id_2)
    glycan_edge <- r[motif_edge]
    # Convert igraph.vs to vector as get_edge_ids expects a simple vector
    # This fix a bug introduced by igraph v2.2.0
    igraph::get_edge_ids(glycan, as.vector(glycan_edge))
  })

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
    ml[["pos2"]] == "?" || all(parse_pos2(gl[["pos2"]]) %in% parse_pos2(ml[["pos2"]]))
  }

  anomer_ok(gl, ml) && pos1_ok(gl, ml) && pos2_ok(gl, ml)
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


anomer_check <- function(r, glycan, motif) {
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
  in_degree <- igraph::degree(glycan, mode = "in")
  igraph::V(glycan)[in_degree == 0]
}
