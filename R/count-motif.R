#' Count How Many Times Glycans have the Given Motif(s)
#'
#' These functions are closely related to [have_motif()].
#' However, instead of returning logical values, they return the number of times
#' the `glycans` have the `motif`(s).
#' - `count_motif()` counts a single motif in multiple glycans
#' - `count_motifs()` counts multiple motifs in multiple glycans
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
#' @seealso [have_motif()], [have_motifs()]
#'
#' @return 
#' - `count_motif()`: An integer vector indicating how many times each `glycan` has the `motif`.
#' - `count_motifs()`: A tibble where the first column 'glycan' contains glycan identifiers,
#'   and subsequent columns contain integer values indicating how many times each glycan has each motif.
#'
#' @examples
#' count_motif("Gal(b1-3)Gal(b1-3)GalNAc", "Gal(b1-")
#' count_motif(
#'   "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc",
#'   "Man(b1-?)[Man(b1-?)]GalNAc"
#' )
#' count_motif("Gal(b1-3)Gal", "Man")
#' 
#' # Vectorized usage with single motif
#' count_motif(c("Gal(b1-3)Gal(b1-3)GalNAc", "Gal(b1-3)GalNAc"), "Gal(b1-")
#'
#' # Multiple motifs with count_motifs()
#' glycan1 <- parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
#' glycan2 <- parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc")
#' glycans <- c(glycan1, glycan2)
#' names(glycans) <- c("double_gal", "complex_man")
#'
#' motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-", "Man(b1-")
#' result <- count_motifs(glycans, motifs)
#' print(result)
#'
#' @export
count_motif <- function(glycans, motif, alignment = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motif_args(glycans, motif, alignment, ignore_linkages)
  rlang::exec("count_motif_", !!!params)
}


count_motif_ <- function(glycans, motif, alignment, ignore_linkages = FALSE) {
  # This function vectorizes `count_single_motif_()`.
  motif_graph <- glyrepr::get_structure_graphs(motif)
  glyrepr::smap_int(glycans, count_single_motif_, motif_graph, alignment, ignore_linkages)
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


#' @rdname count_motif
#' @export
count_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
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
    
    motif_counts <- count_motif(glycans, motif, alignment = alignment, ignore_linkages = ignore_linkages)
    
    # Add column with motif name
    result_tibble[[motif_names[i]]] <- motif_counts
  }
  
  result_tibble
}
