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

  # Apply each motif to all glycans using purrr
  motif_results_list <- purrr::map2(
    motifs,
    alignments,
    ~ single_motif_func(
      glycans, 
      .x, 
      alignment = .y, 
      ignore_linkages = ignore_linkages
    )
  )
  
  # Set names for the results
  names(motif_results_list) <- motif_names
  
  # Convert results to matrix format
  # Each element in motif_results_list should be a vector of results for all glycans
  result_matrix <- do.call(cbind, motif_results_list)
  rownames(result_matrix) <- glycan_names
  colnames(result_matrix) <- motif_names
  
  return(result_matrix)
}
