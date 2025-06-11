#' Count How Many Times Glycans have the Given Motif
#'
#' This function is closely related to [have_motif()].
#' However, instead of returning a logical value, it returns the number of times
#' the `glycans` have the `motif`.
#'
#' @inheritParams have_motif
#'
#' @details
#' This function actually perform v2f algorithm to get all possible matches
#' between `glycans` and `motif`.
#' However, the result is not necessarily the number of matches.
#'
#' Think about the following example:
#' - glycan: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-4)GlcNAc`
#' - motif: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-`
#'
#' To draw the glycan out:
#' ```
#' Gal 1
#'    \ b1-? b1-4
#'     GlcNAc -- GlcNAc
#'    / b1-?
#' Gal 2
#' ```
#' To draw the motif out:
#' ```
#' Gal 1
#'    \ b1-?
#'     GlcNAc b1-
#'    / b1-?
#' Gal 2
#' ```
#'
#' To differentiate the galactoses, we number them as "Gal 1" and "Gal 2"
#' in both the glycan and the motif.
#' The v2f subisomorphic algorithm will return two matches:
#' - Gal 1 in the glycan matches Gal 1 in the motif, and Gal 2 matches Gal 2.
#' - Gal 1 in the glycan matches Gal 2 in the motif, and Gal 2 matches Gal 1.
#'
#' However, from a biological perspective, the two matches are the same.
#' This function will take care of this, and return the "unique" number of matches.
#'
#' For other details about the handling of monosaccharide, linkages, alignment,
#' substituents, and implementation, see [have_motif()].
#'
#' @seealso [have_motif()]
#'
#' @return An integer vector indicating how many times each `glycan` has the `motif`.
#'
#' @examples
#' count_motif("Gal(b1-3)Gal(b1-3)GalNAc", "Gal(b1-")
#' count_motif(
#'   "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc",
#'   "Man(b1-?)[Man(b1-?)]GalNAc"
#' )
#' count_motif("Gal(b1-3)Gal", "Man")
#' 
#' # Vectorized usage
#' count_motif(c("Gal(b1-3)Gal(b1-3)GalNAc", "Gal(b1-3)GalNAc"), "Gal(b1-")
#'
#' @export
count_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motif_args(glycans, motif, alignment, ignore_linkages)
  rlang::exec("count_motif_", !!!params)
}


count_motif_ <- function(glycans, motif, alignment, ignore_linkages = FALSE) {
  # This function vectorizes `count_single_motif_()`.
  motif_graph <- glyrepr::get_structure_graphs(motif)
  glyrepr::structure_map_int(glycans, count_single_motif_, motif_graph, alignment, ignore_linkages)
}


count_single_motif_ <- function(glycan_graph, motif_graph, alignment, ignore_linkages = FALSE) {
  # This function is the logic part of `count_motif()`.
  c_graphs <- colorize_graphs(glycan_graph, motif_graph)
  glycan_graph <- c_graphs$glycan
  motif_graph <- c_graphs$motif
  res <- perform_vf2(glycan_graph, motif_graph)
  valid_mask <- purrr::map_lgl(
    res, is_vaild_result, glycan = glycan_graph, motif = motif_graph,
    alignment = alignment, ignore_linkages = ignore_linkages
  )
  valid_res <- res[valid_mask]
  count_set_unique(valid_res)
}


count_set_unique <- function(lst) {
  # Given a list of integer vectors, return the number of unique vectors.
  # Uniqueness is determined by `setequal()`.
  unique_list <- list()
  for (vec in lst) {
    if (purrr::none(unique_list, ~ setequal(.x, vec))) {
      unique_list <- append(unique_list, list(vec))
    }
  }
  length(unique_list)
}
