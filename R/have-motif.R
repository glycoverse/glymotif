#' Check if the Glycans have the Given Motif
#'
#' @description
#' This function checks if the given `glycan`s have the given `motif`.
#' Technically speaking, it performs a subgraph isomorphism test to
#' determine if the `motif` is a subgraph of the `glycan`.
#' Monosaccharides, linkages , and substituents are all considered.
#'
#' @details
#' # Graph mode and monosaccharide type
#'
#' Both `glycan` and `motif` should be 'glyrepr_structure' objects
#' (see [glyrepr::glycan_structure()]).
#'
#' Also, they can have different monosaccharide types
#' ("concrete" or "generic", see [glyrepr::get_mono_type()]).
#' However, a "concrete" motif cannot be matched to a "generic" glycan.
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
#' Some examples:
#' - "b1-?" in motif will match "b1-4" in glycan.
#' - "b1-?" in motif will match "b1-?" in glycan.
#' - "b1-4" in motif will NOT match "b1-?" in glycan.
#' - "a1-?" in motif will NOT match "b1-4" in glycan.
#' - "a1-?" in motif will NOT match "a?-4" in glycan.
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
#' When using known motifs in the GlycoMotif GlyGen Collection,
#' the best practice is to not provide the `alignment` argument,
#' and let the function decide the alignment based on the motif name.
#' However, it is still possible to override the default alignments.
#' In this case, the user-provided alignments will be used,
#' but a warning will be issued.
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
#' @param glycans A 'glyrepr_structure' object, or an IUPAC-condensed structure string.
#' @param motif A 'glyrepr_structure' object, an IUPAC-condensed structure string,
#' or a known motif name (use [available_motifs()] to see all available motifs).
#' @param alignment A character string.
#' Possible values are "substructure", "core", "terminal" and "whole".
#' If not provided, the value will be decided based on the `motif` argument.
#' If `motif` is a motif name, the alignment in the database will be used.
#' Otherwise, "substructure" will be used.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical value indicating if the `glycan` has the `motif`.
#'
#' @seealso [count_motif()]
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' (glycan <- o_glycan_core_2(mono_type = "concrete"))
#'
#' # The glycan has the motif "Gal(b1-3)GalNAc"
#' have_motif(glycan, "Gal(b1-3)GalNAc")
#'
#' # But not "Gal(b1-4)GalNAc" (wrong linkage)
#' have_motif(glycan, "Gal(b1-4)GalNAc")
#'
#' # Set `ignore_linkages` to `TRUE` to ignore linkages
#' have_motif(glycan, "Gal(b1-4)GalNAc", ignore_linkages = TRUE)
#'
#' # Different monosaccharide types are allowed
#' have_motif(glycan, "Hex(b1-3)HexNAc")
#'
#' # Obscure linkages in the `motif` graph are allowed
#' have_motif(glycan, "Gal(b1-?)GalNAc")
#'
#' # However, obscure linkages in `glycan` will only match "?" in the `motif` graph
#' glycan_2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' have_motif(glycan_2, "Gal(b1-3)GalNAc")
#' have_motif(glycan_2, "Gal(b1-?)GalNAc")
#'
#' # The anomer of the motif will be matched to linkages in the glycan
#' have_motif(glycan_2, "GlcNAc(b1-")
#'
#' # Alignment types
#' # The default type is "substructure", which means the motif can be anywhere in the glycan.
#' # Other options include "core", "terminal" and "whole".
#' glycan_3 <- parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
#' motifs <- c(
#'   "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal",
#'   "Gal(a1-3)Gal(a1-4)Gal",
#'            "Gal(a1-4)Gal(a1-6)Gal",
#'            "Gal(a1-4)Gal"
#' )
#'
#' purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "whole"))
#' purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "core"))
#' purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "terminal"))
#' purrr::map_lgl(motifs, ~ have_motif(glycan_3, .x, alignment = "substructure"))
#'
#' # Substituents
#' glycan_4 <- "Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc"
#' glycan_5 <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc"
#'
#' have_motif(glycan_4, glycan_5)
#' have_motif(glycan_5, glycan_4)
#' have_motif(glycan_4, glycan_4)
#' have_motif(glycan_5, glycan_5)
#'
#' # Vectorization
#' glycans <- c(glycan, glycan_2, glycan_3)
#' motif <- "Gal(b1-3)GalNAc"
#' have_motif(glycans, motif)
#'
#' @export
have_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motif_args(glycans, motif, alignment, ignore_linkages)
  rlang::exec("have_motif_", !!!params)
}


have_motif_ <- function(glycans, motif, alignment, ignore_linkages = FALSE) {
  # This function vectorizes `has_motif_()`.
  motif_graph <- glyrepr::get_structure_graphs(motif)
  glyrepr::smap_lgl(glycans, has_motif_, motif_graph, alignment, ignore_linkages)
}


has_motif_ <- function(glycan_graph, motif_graph, alignment, ignore_linkages = FALSE) {
  # This function is the logic part of `have_motif()`.
  c_graphs <- colorize_graphs(glycan_graph, motif_graph)
  glycan_graph <- c_graphs$glycan
  motif_graph <- c_graphs$motif
  res <- perform_vf2(glycan_graph, motif_graph)
  any(purrr::map_lgl(
    res, is_vaild_result, glycan = glycan_graph, motif = motif_graph,
    alignment = alignment, ignore_linkages = ignore_linkages
  ))
}

#' Check if the Glycans have Multiple Motifs
#'
#' @description
#' This function checks if the given `glycan`s have the given `motif`s.
#' It is a vectorized version of [have_motif()] that accepts multiple motifs.
#' The function returns a logical matrix where rows represent glycans and
#' columns represent motifs.
#'
#' @param glycans A 'glyrepr_structure' object, or an IUPAC-condensed structure string vector.
#' @param motifs A character vector of motif names, IUPAC-condensed structure strings,
#' or a list of 'glyrepr_structure' objects.
#' @param alignments A character vector specifying alignment types for each motif.
#' Possible values are "substructure", "core", "terminal" and "whole".
#' If not provided, the values will be decided based on each motif.
#' If a motif is a known motif name, the alignment in the database will be used.
#' Otherwise, "substructure" will be used.
#' Can be a single value (applied to all motifs) or a vector of the same length as motifs.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A tibble where the first column 'glycan' contains glycan identifiers
#' (names if available, otherwise IUPAC structure strings), and subsequent
#' columns contain logical values indicating whether each glycan has each motif.
#'
#' @seealso [have_motif()]
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' # Create some glycans
#' glycan1 <- o_glycan_core_2(mono_type = "concrete")
#' glycan2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' glycans <- c(glycan1, glycan2)
#'
#' # Define multiple motifs
#' motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc", "GlcNAc(b1-6)GalNAc")
#'
#' # Check which glycans have which motifs
#' result <- have_motifs(glycans, motifs)
#' print(result)
#'
#' # With different alignment types
#' alignments <- c("substructure", "substructure", "core")
#' result2 <- have_motifs(glycans, motifs, alignments = alignments)
#' print(result2)
#'
#' @export
have_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  # Validate inputs
  valid_glycans_arg(glycans)
  if (!is.character(motifs) && !is.list(motifs)) {
    rlang::abort("`motifs` must be a character vector or a list of 'glyrepr_structure' objects.")
  }
  
  if (length(motifs) == 0) {
    rlang::abort("`motifs` cannot be empty.")
  }
  
  # Handle alignments parameter
  if (!is.null(alignments)) {
    if (length(alignments) == 1) {
      alignments <- rep(alignments, length(motifs))
    } else if (length(alignments) != length(motifs)) {
      rlang::abort("`alignments` must be NULL, a single value, or have the same length as `motifs`.")
    }
    # Validate each alignment
    purrr::walk(alignments, valid_alignment_arg)
  }
  
  valid_ignore_linkages_arg(ignore_linkages)
  
  # Save names of the original input before processing
  glycan_names <- names(glycans)
  
  # Ensure glycans are structures
  glycans <- ensure_glycans_are_structures(glycans)
  
  # Prepare motif names for column names
  if (is.character(motifs)) {
    motif_names <- motifs
  } else {
    # For structure objects, create names
    motif_names <- paste0("motif_", seq_along(motifs))
  }
  
  # Create the glycan column similar to describe_n_glycans
  if (!is.null(glycan_names)) {
    glycan_col <- glycan_names
  } else {
    glycan_col <- as.character(glycans)
  }
  
  # Create base tibble with glycan column
  result_tibble <- tibble::tibble(glycan = glycan_col)
  
  # Add columns for each motif
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    alignment <- if (is.null(alignments)) NULL else alignments[i]
    
    motif_results <- have_motif(glycans, motif, alignment = alignment, ignore_linkages = ignore_linkages)
    
    # Add column with motif name
    result_tibble[[motif_names[i]]] <- motif_results
  }
  
  result_tibble
}
