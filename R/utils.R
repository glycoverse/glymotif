# ----- Prepare arguments -----
# Unified function to prepare arguments for both single and multiple motifs
prepare_motif_args <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  single_motif = FALSE,
  strict_sub = TRUE,
  call = rlang::caller_env()
) {
  # Unified validation logic
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  motif_type <- get_motif_type(motifs, call = call)
  alignments <- decide_alignments(motifs, motif_type, alignments)

  glycans <- ensure_glycans_are_structures(glycans, call = call)
  motifs <- ensure_motifs_are_structures(motifs, motif_type, require_scalar = single_motif, call = call)

  # Return appropriate format based on single_motif flag
  if (single_motif) {
    return(list(
      glycans = glycans,
      motif = motifs,
      alignment = alignments,
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub
    ))
  } else {
    return(list(
      glycans = glycans,
      motifs = motifs,
      alignments = alignments,
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub
    ))
  }
}

# Helper function to check for duplicate motifs
# Works with both character vectors and glyrepr_structure objects
has_duplicate_motifs <- function(motifs) {
  if (length(motifs) <= 1) {
    return(FALSE)
  }

  # For glyrepr_structure, use unique() which is natively supported
  if (glyrepr::is_glycan_structure(motifs)) {
    length(unique(motifs)) < length(motifs)
  } else {
    # For character vectors, use base R unique
    length(unique(motifs)) < length(motifs)
  }
}

# Legacy wrapper functions for backward compatibility
prepare_have_motif_args <- function(glycans, motif, alignment, ignore_linkages, strict_sub) {
  prepare_motif_args(glycans, motif, alignment, ignore_linkages, single_motif = TRUE, strict_sub = strict_sub)
}

prepare_have_motifs_args <- function(glycans, motifs, alignments, ignore_linkages, strict_sub) {
  prepare_motif_args(glycans, motifs, alignments, ignore_linkages, single_motif = FALSE, strict_sub = strict_sub)
}

# ----- Argument validation -----
valid_alignment_arg <- function(x) {
  checkmate::assert_choice(x, c("substructure", "core", "terminal", "whole", "exact"), null.ok = TRUE)
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

# Error messages are now centralized in errors.R


motifs_type_err_msg <- paste(
  "`motifs` must be either a `glyrepr_structure` object,",
  "a character vector of IUPAC-condensed structure strings,",
  "or a character vector of known motif names."
)

# ----- Decide motif type -----
# Decide if the `motifs` argument is known motifs,
# an IUPAC-condensed structure character vector, or a 'glyrepr_structure' object.
get_motif_type <- function(motifs, call = rlang::caller_env()) {
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
      cli::cli_abort("Unknown motif: {.val {unknown_motifs}}.", call = call)
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
  # Case 1: `glycans` is already a `glyrepr_structure` object
  if (glyrepr::is_glycan_structure(glycans)) {
    return(glycans)
  }

  # Case 2: `glycans` is a character vector
  # We need to parse it as glycan structure strings (auto-detect format).
  if (is.character(glycans)) {
    tryCatch(
      return(glyparse::auto_parse(glycans)),
      error = function(cnd) {
        cli::cli_abort(c(
          "`glycans` must be a 'glyrepr_structure' object or a glycan structure character vector.",
          "x" = "Some glycans could not be parsed as valid glycan structures."
        ), call = call, parent = cnd)
      }
    )
  }

  # Case 3: `glycans` has other types
  cli::cli_abort(c(
    "`glycans` must be a 'glyrepr_structure' object or a glycan structure character vector.",
    "x" = "The input is of class {.cls {class(glycans)}}."
  ), call = call)
}

# Make sure `motifs` is a `glyrepr_structure` object.
# This function handles both single and multiple motifs
ensure_motifs_are_structures <- function(motifs, motif_type, require_scalar = FALSE, call = rlang::caller_env()) {
  # Check scalar requirement first if needed
  if (require_scalar && length(motifs) != 1) {
    cli::cli_abort(c(
      "The `motif` argument must be a scalar vector.",
      "x" = "The input is of length {.val {length(motifs)}}."
    ), call = call)
  }

  # Determine which error message to use based on scalar requirement
  base_err_msg <- if (require_scalar) {
    paste(
      "`motif` must be either a 'glyrepr_structure' object with length 1,",
      "a glycan structure character scalar,",
      "or a known motif name."
    )
  } else {
    paste0(
      "`motifs` must be a 'glyrepr_structure' object,",
      "a character vector of glycan structure strings,",
      "or a character vector of known motif names."
    )
  }

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

  # Case 4: `motifs` is a vector of glycan structure strings
  if (motif_type == "iupac") {
    tryCatch(
      return(glyparse::auto_parse(motifs)),
      error = function(cnd) {
        cli::cli_abort(c(
          base_err_msg,
          "x" = "Some motifs are neither valid glycan structures nor known motif names.",
          "i" = "Use {.fn available_motifs} to see all valid motif names."
        ), call = call, parent = cnd)
      }
    )
  }
}

# Simplified wrapper function for single motifs
ensure_motif_is_structure <- function(motif, motif_type, call = rlang::caller_env()) {
  ensure_motifs_are_structures(motif, motif_type, require_scalar = TRUE, call = call)
}

# ----- Generic function for single motif mapping -----
apply_single_motif_to_glycans <- function(glycans, motif, alignment, ignore_linkages, strict_sub, single_glycan_func, smap_func) {
  # Generic function to apply a single motif to multiple glycans
  # single_glycan_func should be either .have_motif_single or .count_motif_single
  # smap_func should be either glyrepr::smap_lgl or glyrepr::smap_int

  # Handle mono type conversion based on motif type
  # glyrepr 0.9.0.9000 guarantees all elements in a glyrepr_structure vector
  # have the same mono_type, so get_mono_type() returns a scalar
  motif_type <- glyrepr::get_mono_type(motif)
  if (motif_type == "generic") {
    # For generic motifs, convert glycans to generic to allow matching concrete glycans
    glycans_to_use <- fast_convert_to_generic(glycans)
  } else {
    # For concrete motifs, use glycans as-is
    # (generic glycans will naturally not match due to different mono names)
    glycans_to_use <- glycans
  }

  motif_graph <- glyrepr::get_structure_graphs(motif)
  smap_func(glycans_to_use, single_glycan_func, motif_graph, alignment, ignore_linkages, strict_sub)
}

#' Fast Convert Mono Types
#'
#' `apply_single_motif_to_glycans` used to use `glyrepr::convert_mono_type()` (now `convert_to_generic()`)
#' to ensure that mono types of glycans and motifs are the same.
#' That function uses `glyrepr::spmap_structure()` under the hood,
#' which calls time-consuming `glyrepr:::.structure_to_iupac_single()`
#' to ensure the new `glyrepr::glycan_structure()` has the correct IUPACs.
#' However, in motif mathcing functions of `glymotif`,
#' we don't need to return the converted `glyrepr::glycan_structure()` object.
#' All we need to do is to make sure the underlying glycan graphs are of correct mono types.
#' Therefore, we implement a fast version of `convert_to_generic()` here
#' for this special use case.
#' It only converts the underlying graphs without modifying anything else.
#' The implementation peeks under the hood of `glyrepr`,
#' which is generally a bad decision.
#' But currently `glrepr` cannot fulfill our needs here,
#' and the performance gain is worth the sacrifice.
#'
#' @param glycans A `glyrepr_structure` object.
#' @return A `glyrepr_structure` object with the same IUPACs but generic mono types.
#' @noRd
fast_convert_to_generic <- function(glycans) {
  iupacs <- as.character(glycans)
  unique_iupacs <- names(attr(glycans, "graphs"))
  convert_one_graph <- function(graph) {
    igraph::V(graph)$mono <- glyrepr::convert_to_generic(igraph::V(graph)$mono)
    graph
  }
  new_graphs <- purrr::map(attr(glycans, "graphs"), convert_one_graph)
  glyrepr:::new_glycan_structure(iupacs, new_graphs)
}

# ----- Generic function for multiple motifs -----
prepare_struc_names <- function(x, strucs) {
  if (glyrepr::is_glycan_structure(x)) {
    if (is.null(names(x))) {
      return(as.character(x))
    } else {
      return(names(x))
    }
  } else {  # must be a character vector
    if (is.null(names(x))) {
      return(x)
    } else {
      return(names(x))
    }
  }
}

# Helper to determine motif names for matrix columns
# Returns: explicit names if present, known motif names if applicable,
#          NULL otherwise (IUPAC strings or structure input without names)
prepare_motif_names <- function(motifs_input) {
  # If motifs_input has explicit names, use them
  if (!is.null(names(motifs_input))) {
    return(names(motifs_input))
  }

  # If motifs_input is a character vector of known motif names, use those
  if (is.character(motifs_input) && all(is_known_motif(motifs_input))) {
    return(motifs_input)
  }

  # Otherwise, no names (IUPAC strings or glyrepr_structure without names)
  NULL
}

apply_motifs_to_glycans <- function(glycans, motifs, alignments, ignore_linkages, single_motif_func, glycan_names, motif_names, strict_sub) {
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
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub
    )
  )

  # Set names for the results if provided
  if (!is.null(motif_names)) {
    names(motif_results_list) <- motif_names
  }

  # Convert results to matrix format
  # Each element in motif_results_list should be a vector of results for all glycans
  result_matrix <- do.call(cbind, motif_results_list)
  
  # Set rownames if provided
  if (!is.null(glycan_names)) {
    rownames(result_matrix) <- glycan_names
  }
  
  # Set colnames if provided
  if (!is.null(motif_names)) {
    colnames(result_matrix) <- motif_names
  }

  return(result_matrix)
}
