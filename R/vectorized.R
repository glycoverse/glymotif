#' @rdname have_motif
#' @export
have_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motifs_args(glycans, motifs, alignments, ignore_linkages)
  rlang::exec("have_motifs_", !!!params)
}

have_motifs_ <- function(glycans, motifs, alignments, ignore_linkages = FALSE) {
  apply_motifs_to_glycans(glycans, motifs, alignments, ignore_linkages, have_motif_)
}

#' @rdname count_motif
#' @export
count_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motifs_args(glycans, motifs, alignments, ignore_linkages)
  rlang::exec("count_motifs_", !!!params)
}

count_motifs_ <- function(glycans, motifs, alignments, ignore_linkages = FALSE) {
  apply_motifs_to_glycans(glycans, motifs, alignments, ignore_linkages, count_motif_)
}

# ----- Generic function for multiple motifs -----
apply_motifs_to_glycans <- function(glycans, motifs, alignments, ignore_linkages, single_motif_func) {
  # Generic function to apply multiple motifs to multiple glycans
  # single_motif_func should be either have_motif_ or count_motif_

  glycan_names <- names(glycans)
  motif_names <- names(motifs)

  # Create base tibble with glycan column
  result_tibble <- tibble::tibble(glycan = glycan_names)
  
  # Add columns for each motif
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
    alignment <- if (is.null(alignments)) NULL else alignments[i]
    
    motif_results <- single_motif_func(glycans, motif, alignment = alignment, ignore_linkages = ignore_linkages)
    
    # Add column with motif name
    result_tibble[[motif_names[i]]] <- motif_results
  }
  
  result_tibble
}
