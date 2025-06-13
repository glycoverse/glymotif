#' @rdname have_motif
#' @export
have_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motifs_args(glycans, motifs, alignments, ignore_linkages)
  glycan_names <- prepare_struc_names(glycans, params$glycans)
  motif_names <- prepare_struc_names(motifs, params$motifs)
  rlang::exec("have_motifs_", !!!params, glycan_names = glycan_names, motif_names = motif_names)
}

have_motifs_ <- function(glycans, motifs, alignments, glycan_names, motif_names, ignore_linkages = FALSE) {
  apply_motifs_to_glycans(glycans, motifs, alignments, ignore_linkages, have_motif_, glycan_names, motif_names)
}

#' @rdname count_motif
#' @export
count_motifs <- function(glycans, motifs, alignments = NULL, ignore_linkages = FALSE) {
  params <- prepare_have_motifs_args(glycans, motifs, alignments, ignore_linkages)
  glycan_names <- prepare_struc_names(glycans, params$glycans)
  motif_names <- prepare_struc_names(motifs, params$motifs)
  rlang::exec("count_motifs_", !!!params, glycan_names = glycan_names, motif_names = motif_names)
}

count_motifs_ <- function(glycans, motifs, alignments, glycan_names, motif_names, ignore_linkages = FALSE) {
  apply_motifs_to_glycans(glycans, motifs, alignments, ignore_linkages, count_motif_, glycan_names, motif_names)
}

# ----- Generic function for multiple motifs -----
prepare_struc_names <- function(x, strucs) {
  if (glyrepr::is_glycan_structure(x)) {
    return(as.character(x))
  } else {  # must be a character vector
    if (is.null(names(x))) {
      return(x)
    } else {
      return(names(x))
    }
  }
}

apply_motifs_to_glycans <- function(glycans, motifs, alignments, ignore_linkages, single_motif_func, glycan_names, motif_names) {
  # Generic function to apply multiple motifs to multiple glycans
  # single_motif_func should be either have_motif_ or count_motif_

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
