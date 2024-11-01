#' Check if a Glycan has the Given Motif
#'
#' @description
#' This function checks if a `glycan` has a given `motif`.
#' Technically speaking, it performs a subgraph isomorphism test to
#' determine if the `motif` is a subgraph of the `glycan`.
#' Monosaccharides, linkages , and substituents are all considered.
#'
#' @details
#' # Graph mode and monosaccharide type
#'
#' Both `glycan` and `motif` should be 'glycan_graph' objects
#' (see [glyrepr::as_glycan_graph()]).
#' They can be either "NE" or "DN" glycan graphs (can be different).
#'
#' Also, they can have different monosaccharide types
#' ("concrete", "generic" or "simple", see [glyrepr::decide_mono_type()]).
#' However, the monosaccharide type of `glycan` cannot be obscurer than that of `motif`,
#' which will raise an error.
#' For example, a "concrete" `glycan` can have a "generic" `motif`, but not vice versa.
#'
#' # Linkages
#'
#' Obscure linkages (e.g. "??-?") are allowed in the `motif` graph
#' (see [glyrepr::possible_linkages()]).
#' "?" in a motif graph means "anything could be OK",
#' so it will match any linkage in the `glycan` graph.
#' However, "?" in a `glycan` graph will only match "?" in the `motif` graph.
#' You can set `ignore_linkages = TRUE` to ignore linkages in the comparison.
#'
#' Both motifs and glycans can have a "half-linkage" at the reducing end,
#' e.g. "GlcNAc(b1-".
#' The half linkage in the motif will be matched to any linkage in the glycan,
#' or the half linkage of the glycan.
#' e.g. Glycan "GlcNAc(b1-4)Gal(a1-" will have both "GlcNAc(b1-" and "Gal(a1-" motifs.
#'
#' # Alignment
#'
#' According to the [GlycoMotif](https://glycomotif.glyomics.org) database,
#' a motif can be classified into four alignment types:
#' - "substructure": The motif can be anywhere in the glycan. This is the default.
#' See [substructure](https://glycomotif.glyomics.org/glycomotif/Substructure_Alignment)
#' for details.
#' - "core": The motif must align with at least one connected substructure
#' (subtree) at the reducing end of the glycan.
#' See [glycan core](https://glycomotif.glyomics.org/glycomotif/Glycan_Core_Alignment)
#' for details.
#' - "terminal": The motif must align with at least one connected substructure
#' (subtree) at the nonreducing end of the glycan.
#' See [nonreducing end](https://glycomotif.glyomics.org/glycomotif/Nonreducing-End_Alignment)
#' for details.
#' - "whole": The motif must align with the entire glycan.
#' See [whole-glycan](https://glycomotif.glyomics.org/glycomotif/Whole-Glycan_Alignment)
#' for details.
#'
#' # Substituents
#'
#' Substituents (e.g. "Ac", "SO3") are matched in strict mode.
#' "Neu5Ac-9Ac" will only match "Neu5Ac-9Ac" but not "Neu5Ac",
#' and "Neu5Ac" will not match "Neu5Ac-9Ac".
#' Obscure linkage in the motif will match any linkage in the glycan.
#' e.g. Motif "Neu5Ac-?Ac" will match "Neu5Ac-9Ac" in the glycan.
#'
#' # Implementation
#'
#' Under the hood, the function uses [igraph::graph.get.subisomorphisms.vf2()]
#' to get all possible subgraph isomorphisms between `glycan` and `motif`.
#' `color` vertex attributes are added to the graphs to distinguish monosaccharides.
#' For all possible matches, the function checks the following:
#' - Alignment: using `alignment_check()`
#' - Substituents: using `substituent_check()`
#' - Linkages: using `linkage_check()`
#' - Anomer: using `anomer_check()`
#' The function returns `TRUE` if any of the matches pass all checks.
#'
#' @param glycan A 'glycan_graph' object.
#' @param motif A 'glycan_graph' object.
#' @param ... Not used.
#' @param alignment A character string. Possible values are "substructure", "core", "terminal" and "whole".
#' See description for details. Default is "substructure".
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical value indicating if the `glycan` has the `motif`.
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' (glycan <- glyrepr::o_glycan_core_2(mode = "ne", mono_type = "concrete"))
#'
#' # The glycan has the motif "Gal(b1-3)GalNAc"
#' motif_1 <- parse_iupac_condensed("Gal(b1-3)GalNAc")
#' has_motif(glycan, motif_1)
#'
#' # But not "Gal(b1-4)GalNAc" (wrong linkage)
#' motif_2 <- parse_iupac_condensed("Gal(b1-4)GalNAc")
#' has_motif(glycan, motif_2)
#'
#' # Set `ignore_linkages` to `TRUE` to ignore linkages
#' has_motif(glycan, motif_2, ignore_linkages = TRUE)
#'
#' # Glycan and motif can have different graph modes
#' motif_1_dn <- convert_graph_mode(motif_1, to = "dn")
#' has_motif(glycan, motif_1_dn)
#'
#' # And different monosaccharide types
#' motif_1_generic <- convert_glycan_mono_type(motif_1, "generic")
#' has_motif(glycan, motif_1_generic)
#'
#' # However, the monosaccharide type of `glycan` cannot be obscurer than that of `motif`
#' glycan_simple <- convert_glycan_mono_type(glycan, "simple")
#' try(has_motif(glycan_simple, motif_1))
#'
#' # Obscure linkages in the `motif` graph are allowed
#' motif_3 <- parse_iupac_condensed("Gal(b1-?)GalNAc")
#' has_motif(glycan, motif_3)
#'
#' # However, obscure linkages in `glycan` will only match "?" in the `motif` graph
#' glycan_2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' has_motif(glycan_2, motif_1)
#' has_motif(glycan_2, motif_3)
#'
#' # The anomer of the motif will be matched to linkages in the glycan
#' motif_4 <- parse_iupac_condensed("GlcNAc(b1-")
#' has_motif(glycan_2, motif_4)
#'
#' # Alignment types
#' # The default type is "substructure", which means the motif can be anywhere in the glycan.
#' # Other options include "core", "terminal" and "whole".
#' glycan_3 <- parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_5 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motif_6 <-  parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal")
#' motif_7 <-  parse_iupac_condensed(         "Gal(a1-4)Gal(a1-6)Gal")
#' motif_8 <-  parse_iupac_condensed(         "Gal(a1-4)Gal")
#' motifs <- list(motif_5, motif_6, motif_7, motif_8)
#'
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "whole"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "core"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "terminal"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "substructure"))
#'
#' # Substituents
#' glycan_4 <- parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc")
#' glycan_5 <- parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
#'
#' has_motif(glycan_4, glycan_5)
#' has_motif(glycan_5, glycan_4)
#' has_motif(glycan_4, glycan_4)
#' has_motif(glycan_5, glycan_5)
#'
#' @export
has_motif <- function(glycan, motif, ..., alignment = "substructure", ignore_linkages = FALSE) {
  # Check input arguments
  if (!glyrepr::is_glycan(glycan) || !glyrepr::is_glycan(motif)) {
    rlang::abort("`glycan` and `motif` must be 'glycan_graph' objects.")
  }
  if (!alignment %in% c("substructure", "core", "terminal", "whole")) {
    rlang::abort("`alignment` must be one of 'substructure', 'core', 'terminal' or 'whole'.")
  }

  # Ensure that `glycan` and `motif` are all "NE" type
  glycan <- glyrepr::convert_graph_mode(glycan, to = "ne", strict = FALSE)
  motif <- glyrepr::convert_graph_mode(motif, to = "ne", strict = FALSE)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycan_mono_type()`.
  glycan <- ensure_glycan_mono_type(glycan, motif)

  # Colorize the graphs
  c_graphs <- colorize_graphs(glycan, motif)
  glycan <- c_graphs$glycan
  motif <- c_graphs$motif

  # Perform "VF2" algorithm
  # `res` is a list of all possible matches between `glycan` and `motif`.
  # Each match is an integer vector, with the same length as `vcount(motif)`.
  # The `i`-th element of the vector is the index of the vertex in `glycan`
  # that matches the `i`-th vertex in `motif`.
  # e.g. If `res[[1]] = c(2, 3)`,
  # it means the 1st vertex in `motif` matches the 2nd vertex in `glycan`,
  # and the 2nd vertex in `motif` matches the 3rd vertex in `glycan`.
  res <- igraph::graph.get.subisomorphisms.vf2(glycan, motif)

  # Check each result in `res` for the following:
  # alignment, linkages, anomer, substituents
  any(purrr::map_lgl(
    res, is_vaild_result, glycan = glycan, motif = motif,
    alignment = alignment, ignore_linkages = ignore_linkages
  ))
}


ensure_glycan_mono_type <- function(glycan, motif) {
  # Make the glycan the same monosaccharide type (concrete, generic or simple) as the motif.
  # If the glycan type is lower than the motif type, an error will be raised.
  glycan_type <- glyrepr::decide_glycan_mono_type(glycan)
  motif_type <- glyrepr::decide_glycan_mono_type(motif)

  conditions <- (
    (glycan_type == "concrete" && motif_type != "concrete") ||
    (glycan_type == "generic" && motif_type == "simple")
  )
  if (conditions) return(glyrepr::convert_glycan_mono_type(glycan, motif_type))

  if (glycan_type != motif_type) {
    cli::cli_abort(c(
      "The monosaccharide type of `glycan` cannot be obscurer than `motif`.",
      "x" = "{.val {glycan_type}} is obscurer than {.val {motif_type}}."
    ))
  }
  return(glycan)
}


colorize_graphs <- function(glycan, motif) {
  # Add "color" vevtex attributes to the graph.
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


is_vaild_result <- function(r, glycan, motif, alignment, ignore_linkages) {
  # This function actually allocate each check to different functions,
  # including: `alignment_check`, `substituent_check`, `linkage_check`, `anomer_check`.
  # If any of them returns `FALSE`, the result is invalid.
  if (ignore_linkages) {
    check_funs <- list(alignment_check, substituent_check)
  } else {
    check_funs <- list(alignment_check, substituent_check, linkage_check, anomer_check)
  }
  all(purrr::map_lgl(check_funs, ~ .x(r, glycan, motif, alignment = alignment)))
}


alignment_check <- function(r, glycan, motif, alignment) {
  if (alignment == "substructure") return(TRUE)

  if (alignment == "whole") {
    return(igraph::isomorphic(glycan, motif, method = "vf2"))
  }

  if (alignment == "core") {
    glycan_core <- core_node(glycan)
    motif_core <- core_node(motif)
    return(r[[motif_core]] == glycan_core)
  }

  if (alignment == "terminal") {
    glycan_terminals <- terminal_nodes(glycan)
    motif_terminals <- terminal_nodes(motif)
    return(all(r[motif_terminals] %in% glycan_terminals))
  }
}


substituent_check <- function(r, glycan, motif, ...) {
  glycan_subs <- igraph::V(glycan)$sub[r]
  motif_subs <- igraph::V(motif)$sub
  all(purrr::map2_lgl(glycan_subs, motif_subs, match_sub))
}


match_sub <- function(glycan_sub, motif_sub) {
  if (stringr::str_sub(motif_sub, 1, 1) == "?") {
    stringr::str_sub(glycan_sub, 2) == stringr::str_sub(motif_sub, 2)
  } else {
    glycan_sub == motif_sub
  }
}


linkage_check <- function(r, glycan, motif, ...) {
  edges <- get_corresponding_edges(r, glycan, motif)
  glycan_linkages <- edges$glycan$linkage
  motif_linkages <- edges$motif$linkage
  all(purrr::map2_lgl(glycan_linkages, motif_linkages, match_linkage))
}


get_corresponding_edges <- function(r, glycan, motif) {
  motif_edge_list <- igraph::as_edgelist(motif, names = FALSE)

  glycan_edge_ids <- vector("integer", length = nrow(motif_edge_list))
  for (i in seq_len(nrow(motif_edge_list))) {
    motif_edge <- motif_edge_list[i, ]  # c(node_id_1, node_id_2)
    glycan_edge <- r[motif_edge]
    glycan_edge_ids[[i]] <- igraph::get.edge.ids(glycan, glycan_edge)
  }

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
    pos2   = stringr::str_sub(linkage, 4, 4)
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
