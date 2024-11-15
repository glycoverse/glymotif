#' Count How Many Times a Glycan has the Given Motif
#'
#' This function is closely related to [has_motif()].
#' However, instead of returning a logical value, it returns the number of times
#' the `glycan` has the `motif`.
#'
#' @inheritParams has_motif
#'
#' @details
#' This function actually perform v2f algorithm to get all possible matches
#' between `glycan` and `motif`.
#' However, the result is not necessarily the number of matches.
#'
#' Think about the following example:
#' - glycan: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-4)GlcNAc`
#' - motif: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-)`
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
#' substituents, and implementation, see [has_motif()].
#'
#' @seealso [has_motif()]
#'
#' @return An integer indicating how many times the `glycan` has the `motif`.
#'
#' @examples
#' counts_motif("Gal(b1-3)Gal(b1-3)GalNAc", "Gal(b1-")
#' counts_motif(
#'   "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc",
#'   "Man(b1-?)[Man(b1-?)]GalNAc"
#' )
#' counts_motif("Gal(b1-3)Gal", "Man")
#'
#' @export
counts_motif <- function(glycan, motif, alignment = NULL, ignore_linkages = FALSE) {
  params <- prepare_has_motif_args(glycan, motif, alignment, ignore_linkages)
  rlang::exec("counts_motif_", !!!params)
}


counts_motif_ <- function(glycan, motif, alignment, ignore_linkages = FALSE) {
  # This function is the logic part of `counts_motif()`.
  c_graphs <- colorize_graphs(glycan, motif)
  glycan <- c_graphs$glycan
  motif <- c_graphs$motif
  res <- perform_vf2(glycan, motif)
  valid_mask <- purrr::map_lgl(
    res, is_vaild_result, glycan = glycan, motif = motif,
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
