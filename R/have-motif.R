#' Check if the Glycans have the Given Motif(s)
#'
#' @description
#' These functions check if the given `glycan`s have the given `motif`(s).
#' - `have_motif()` checks a single motif against multiple glycans
#' - `have_motifs()` checks multiple motifs against multiple glycans
#' 
#' Technically speaking, they perform subgraph isomorphism tests to
#' determine if the `motif`(s) are subgraphs of the `glycan`s.
#' Monosaccharides, linkages, and substituents are all considered.
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
#' @param glycans A 'glyrepr_structure' object, or an IUPAC-condensed structure string vector.
#' @param motif A 'glyrepr_structure' object, an IUPAC-condensed structure string,
#' or a known motif name (use [available_motifs()] to see all available motifs).
#' @param motifs A character vector of motif names, IUPAC-condensed structure strings,
#' or a 'glyrepr_structure' object.
#' @param alignment A character string.
#' Possible values are "substructure", "core", "terminal" and "whole".
#' If not provided, the value will be decided based on the `motif` argument.
#' If `motif` is a motif name, the alignment in the database will be used.
#' Otherwise, "substructure" will be used.
#' @param alignments A character vector specifying alignment types for each motif.
#' Can be a single value (applied to all motifs) or a vector of the same length as motifs.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return 
#' - `have_motif()`: A logical vector indicating if each `glycan` has the `motif`.
#' - `have_motifs()`: A logical matrix where rows correspond to glycans and columns correspond to motifs.
#'   Row names contain glycan identifiers and column names contain motif identifiers.
#'
#' @seealso [count_motif()], [count_motifs()]
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
#' # Vectorization with single motif
#' glycans <- c(glycan, glycan_2, glycan_3)
#' motif <- "Gal(b1-3)GalNAc"
#' have_motif(glycans, motif)
#'
#' # Multiple motifs with have_motifs()
#' glycan1 <- o_glycan_core_2(mono_type = "concrete")
#' glycan2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' glycans <- c(glycan1, glycan2)
#'
#' motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc", "GlcNAc(b1-6)GalNAc")
#' result <- have_motifs(glycans, motifs)
#' print(result)
#'
#' @export
have_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motif_args(glycans, motif, alignment, ignore_linkages)
  rlang::exec("have_motif_", !!!params)
}


have_motif_ <- function(glycans, motif, alignment, ignore_linkages = FALSE) {
  # This function is a simpler version of `have_motif()`.
  # It performs the logic directly without argument validations and conversions.
  # It is a necessary abstraction for other functions.
  apply_single_motif_to_glycans(
    glycans = glycans,
    motif = motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    single_glycan_func = .have_motif_single,
    smap_func = glyrepr::smap_lgl
  )
}


.have_motif_single <- function(glycan_graph, motif_graph, alignment, ignore_linkages = FALSE) {
  # Optimized version with early termination
  # Check if any match is valid, returning immediately on first valid match
  c_graphs <- colorize_graphs(glycan_graph, motif_graph)
  glycan_graph <- c_graphs$glycan
  motif_graph <- c_graphs$motif
  res <- perform_vf2(glycan_graph, motif_graph)
  purrr::some(
    res, is_valid_result, glycan = glycan_graph, motif = motif_graph,
    alignment = alignment, ignore_linkages = ignore_linkages
  )
}

# ----- Generic function for single motif mapping -----
apply_single_motif_to_glycans <- function(glycans, motif, alignment, ignore_linkages, single_glycan_func, smap_func) {
  # Generic function to apply a single motif to multiple glycans
  # single_glycan_func should be either .have_motif_single or .count_motif_single
  # smap_func should be either glyrepr::smap_lgl or glyrepr::smap_int
  motif_graph <- glyrepr::get_structure_graphs(motif)
  smap_func(glycans, single_glycan_func, motif_graph, alignment, ignore_linkages)
}
