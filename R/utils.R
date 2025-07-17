# ----- Prepare arguments -----
prepare_have_motifs_args <- function(glycans, motifs, alignments, ignore_linkages) {
  # for `have_motifs()`, `count_motifs()`
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  motif_type <- get_motif_type(motifs)
  alignments <- decide_alignments(motifs, motif_type, alignments)

  glycans <- ensure_glycans_are_structures(glycans)
  motifs <- ensure_motifs_are_structures(motifs, motif_type)

  list(glycans = glycans, motifs = motifs, alignments = alignments, ignore_linkages = ignore_linkages)
}

prepare_have_motif_args <- function(glycans, motif, alignment, ignore_linkages) {
  # for `have_motif()`, `count_motif()`
  valid_alignments_arg(alignment, motif)
  valid_ignore_linkages_arg(ignore_linkages)

  motif_type <- get_motif_type(motif)
  alignment <- decide_alignments(motif, motif_type, alignment)

  glycans <- ensure_glycans_are_structures(glycans)
  motif <- ensure_motif_is_structure(motif, motif_type)

  list(glycans = glycans, motif = motif, alignment = alignment, ignore_linkages = ignore_linkages)
}

# ----- Argument validation -----
valid_alignment_arg <- function(x) {
  checkmate::assert(
    checkmate::test_null(x),
    checkmate::test_choice(x, c("substructure", "core", "terminal", "whole"))
  )
}

valid_alignments_arg <- function(alignments, motifs) {
  # Validate alignments parameter without modification
  if (!is.null(alignments)) {
    if (length(alignments) != 1 && length(alignments) != length(motifs)) {
      rlang::abort("`alignments` must be NULL, a single value, or have the same length as `motifs`.")
    }
    # Validate each alignment
    purrr::walk(alignments, valid_alignment_arg)
  }
}

valid_ignore_linkages_arg <- function(x) {
  checkmate::assert_flag(x)
}

glycan_type_err_msg <- paste(
  "`glycans` must be a 'glyrepr_structure' object",
  "or an IUPAC-condensed structure character."
)


motif_type_err_msg <- paste(
  "`motif` must be either a 'glyrepr_structure' object with length 1,",
  "an IUPAC-condensed structure character scalar,",
  "or a known motif name."
)


motifs_type_err_msg <- paste(
  "`motifs` must be either a `glyrepr_structure` object,",
  "a character vector of IUPAC-condensed structure strings,",
  "or a character vector of known motif names."
)

# ----- Decide motif type -----
# Decide if the `motifs` argument is known motifs,
# an IUPAC-condensed structure character vector, or a 'glyrepr_structure' object.
get_motif_type <- function(motifs) {
  # If it is neither a glycan graph or a known motif, it is assumed to
  # be an IUPAC-condensed structure string.
  # This assumption may not be correct, for it is possible that a wrong
  # motif name is passed in.
  if (glyrepr::is_glycan_structure(motifs)) {
    return("structure")
  } else if (is.character(motifs)) {
    known_motif_idx <- is_known_motif(motifs)
    if (all(known_motif_idx)) {
      return("known")
    } else if (any(known_motif_idx)) {
      unknown_motifs <- motifs[!known_motif_idx]
      cli::cli_abort("Unknown motif: {.val {unknown_motifs}}.")
    } else {
      return("iupac")
    }
  } else {
    return(NA_character_)
  }
}

# ----- Decide alignment type -----
decide_alignments <- function(motifs, motif_type, alignments) {
  # Decide which alignment type to use.
  # 1. If the motif is a known motif name: use the alignment type in the database unless provided.
  # 2. If the motif is an IUPAC string or a 'glyrepr_structure' object: use "substructure" unless provided.
  # In the first case, if the user-provided alignment is different from that in the database,
  # issue a warning.

  # Handle NA motif_type (invalid motif type)
  if (is.na(motif_type)) {
    if (is.null(alignments)) {
      alignments <- "substructure"
    }
  } else if (motif_type == "known") {
    db_alignments <- get_motif_alignment(motifs)
    if (!is.null(alignments)) {
      inconsistent <- alignments != db_alignments
      if (any(inconsistent)) {
        warn_msg <- paste(
          "The provided alignment type {.val {alignments[inconsistent]}} is different from",
          "the motif's alignment type {.val {db_alignments[inconsistent]}} in database",
          "for motif {.val {motifs[inconsistent]}}."
        )
        cli::cli_warn(warn_msg, class = "warning_custom_alignment")
      }
    } else {
      alignments <- db_alignments
    }
  } else {
    if (is.null(alignments)) {
      alignments <- "substructure"
    }
  }
  vctrs::vec_recycle(alignments, length(motifs))
}

# ----- Make sure glycans and motifs are structures -----
# Make sure `glycans` is a `glyrepr_structure` object.
ensure_glycans_are_structures <- function(glycans, call = rlang::caller_env()) {
  base_err_msg <- "`glycans` must be a 'glyrepr_structure' object or an IUPAC-condensed structure character."

  # Case 1: `glycans` is already a `glyrepr_structure` object
  if (glyrepr::is_glycan_structure(glycans)) {
    return(glycans)
  }

  # Case 2: `glycans` is a character vector
  # We need to parse it as IUPAC-condensed structure strings.
  if (is.character(glycans)) {
    tryCatch(
      return(glyparse::parse_iupac_condensed(glycans)),
      error = function(cnd) {
        cli::cli_abort(c(
          base_err_msg,
          "x" = "Some glycans could not be parsed as valid IUPAC-condensed structures."
        ), call = call, parent = cnd)
      }
    )
  }

  # Case 3: `glycans` has other types
  cli::cli_abort(c(
    base_err_msg,
    "x" = "The input is of class {.cls {class(glycans)}}."
  ), call = call)
}

# Make sure `motifs` is a `glyrepr_structure` object.
ensure_motifs_are_structures <- function(motifs, motif_type, call = rlang::caller_env()) {
  base_err_msg <- paste0(
    "`motifs` must be a 'glyrepr_structure' object,",
    "a character vector of IUPAC-condensed structure strings,",
    "or a character vector of known motif names."
  )

  # Case 1: `motifs` is of other types
  if (is.na(motif_type)) {
    cli::cli_abort(c(
      base_err_msg,
      "x" = "The input is of class {.cls {class(motifs)}}."
    ), call = call)
  }

  # Case 2: `motifs` is already a `glyrepr_structure` object
  if (motif_type == "structure") {
    return(motifs)
  }

  # Case 3: `motifs` is a vector of known motif names
  if (motif_type == "known") {
    tryCatch(
      return(get_motif_structure(motifs)),
      error = function(cnd) {
        cli::cli_abort(c(
          base_err_msg,
          "i" = "Use {.fn available_motifs} to see all valid motif names."
        ), call = call, parent = cnd)
      }
    )
  }

  # Case 4: `motifs` is a vector of IUPAC-condensed structure strings
  if (motif_type == "iupac") {
    tryCatch(
      return(glyparse::parse_iupac_condensed(motifs)),
      error = function(cnd) {
        cli::cli_abort(c(
          base_err_msg,
          "x" = "Some motifs are neither valid IUPAC-condensed structures nor known motif names.",
          "i" = "Use {.fn available_motifs} to see all valid motif names."
        ), call = call, parent = cnd)
      }
    )
  }
}

ensure_motif_is_structure <- function(motif, motif_type, call = rlang::caller_env()) {
  if (!length(motif) == 1) {
    cli::cli_abort(c(
      "The `motif` argument must be a scalar vector.",
      "x" = "The input is of length {.val {length(motif)}}."
    ), call = call)
  }
  ensure_motifs_are_structures(motif, motif_type, call)
}

# ----- Check and align mono types -----
ensure_glycans_mono_type <- function(glycans, motif) {
  # Make sure all glycans have the same monosaccharide type as the motif.
  # This function assumes that all glycans have the same monosaccharide type.
  glycan_type <- glyrepr::get_mono_type(glycans[[1]])
  motif_type <- glyrepr::get_mono_type(motif)

  if (glycan_type == "generic" && motif_type == "concrete") {
    cli::cli_abort("`generic` glycans cannot be compared with `concrete` motifs.")
  }
  glyrepr::convert_mono_type(glycans, to = motif_type)
}


same_mono_types <- function(structures) {
  # Check if all structures have the same monosaccharide type.
  mono_types <- glyrepr::get_mono_type(structures)
  dplyr::n_distinct(mono_types) <= 1
}

# ----- Generic function for single motif mapping -----
apply_single_motif_to_glycans <- function(glycans, motif, alignment, ignore_linkages, single_glycan_func, smap_func) {
  # Generic function to apply a single motif to multiple glycans
  # single_glycan_func should be either .have_motif_single or .count_motif_single
  # smap_func should be either glyrepr::smap_lgl or glyrepr::smap_int
  
  # Handle mono type conversion based on motif type
  motif_type <- glyrepr::get_mono_type(motif)
  if (motif_type == "generic") {
    # For generic motifs, convert glycans to generic to allow matching concrete glycans
    glycans_to_use <- glyrepr::convert_mono_type(glycans, to = "generic")
  } else {
    # For concrete motifs, use glycans as-is 
    # (generic glycans will naturally not match due to different mono names)
    glycans_to_use <- glycans
  }
  
  motif_graph <- glyrepr::get_structure_graphs(motif)
  smap_func(glycans_to_use, single_glycan_func, motif_graph, alignment, ignore_linkages)
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

  # Handle empty motifs case
  if (length(motifs) == 0) {
    cli::cli_abort("`motifs` cannot be empty.")
  }

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
