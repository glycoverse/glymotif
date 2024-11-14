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
#' @param glycan A 'glycan_graph' object, or an IUPAC-condensed structure string.
#' @param motif A 'glycan_graph' object, an IUPAC-condensed structure string,
#' or a known motif name (use [available_motifs()] to see all available motifs).
#' @param alignment A character string. Possible values are "substructure", "core", "terminal" and "whole".
#' Default is "substructure".
#' When `motif` is a known motif name and `alignment` is not provided,
#' the alignment type in the database will be used.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'
#' @return A logical value indicating if the `glycan` has the `motif`.
#'
#' @seealso [has_motifs()], [have_motif()], [have_motifs()], [counts_motif()]
#'
#' @examples
#' library(glyparse)
#' library(glyrepr)
#'
#' (glycan <- o_glycan_core_2(mode = "ne", mono_type = "concrete"))
#'
#' # The glycan has the motif "Gal(b1-3)GalNAc"
#' has_motif(glycan, "Gal(b1-3)GalNAc")
#'
#' # But not "Gal(b1-4)GalNAc" (wrong linkage)
#' has_motif(glycan, "Gal(b1-4)GalNAc")
#'
#' # Set `ignore_linkages` to `TRUE` to ignore linkages
#' has_motif(glycan, "Gal(b1-4)GalNAc", ignore_linkages = TRUE)
#'
#' # Different monosaccharide types are allowed
#' has_motif(glycan, "Hex(b1-3)HexNAc")
#'
#' # However, the monosaccharide type of `glycan` cannot be obscurer than that of `motif`
#' glycan_simple <- convert_glycan_mono_type(glycan, "simple")
#' try(has_motif(glycan_simple, "Gal(b1-3)GalNAc"))
#'
#' # Obscure linkages in the `motif` graph are allowed
#' has_motif(glycan, "Gal(b1-?)GalNAc")
#'
#' # However, obscure linkages in `glycan` will only match "?" in the `motif` graph
#' glycan_2 <- parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
#' has_motif(glycan_2, "Gal(b1-3)GalNAc")
#' has_motif(glycan_2, "Gal(b1-?)GalNAc")
#'
#' # The anomer of the motif will be matched to linkages in the glycan
#' has_motif(glycan_2, "GlcNAc(b1-")
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
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "whole"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "core"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "terminal"))
#' purrr::map_lgl(motifs, ~ has_motif(glycan_3, .x, alignment = "substructure"))
#'
#' # Substituents
#' glycan_4 <- "Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc"
#' glycan_5 <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc"
#'
#' has_motif(glycan_4, glycan_5)
#' has_motif(glycan_5, glycan_4)
#' has_motif(glycan_4, glycan_4)
#' has_motif(glycan_5, glycan_5)
#'
#' @export
has_motif <- function(glycan, motif, alignment = "substructure", ignore_linkages = FALSE) {
  alignment_provided <- !missing(alignment)
  params <- prepare_has_motif_args(glycan, motif, alignment, alignment_provided, ignore_linkages)
  rlang::exec("has_motif_", !!!params)
}


has_motif_ <- function(glycan, motif, alignment = "substructure", ignore_linkages = FALSE) {
  # This function is the logic part of `has_motif()`.
  c_graphs <- colorize_graphs(glycan, motif)
  glycan <- c_graphs$glycan
  motif <- c_graphs$motif
  res <- perform_vf2(glycan, motif)
  any(purrr::map_lgl(
    res, is_vaild_result, glycan = glycan, motif = motif,
    alignment = alignment, ignore_linkages = ignore_linkages
  ))
}
