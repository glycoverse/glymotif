# ----- Prepare arguments -----
#' Prepare Motif Arguments
#'
#' Normalizes glycan, motif, alignment, linkage, substituent, and match-degree
#' arguments for public motif-matching entry points.
#'
#' @param glycans A glycan structure vector or parseable character vector.
#' @param motifs Motifs, database motif names, or a motif specification object.
#' @param alignments Optional alignment values.
#' @param ignore_linkages Whether linkage matching should be ignored.
#' @param match_degree Optional degree-matching mask.
#' @param single_motif Whether the caller expects exactly one motif.
#' @param strict_sub Whether substituent matching should be strict.
#' @param call The call environment used for errors.
#'
#' @return A normalized argument list for single- or multiple-motif internals.
#' @noRd
prepare_motif_args <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  match_degree = NULL,
  single_motif = FALSE,
  strict_sub = TRUE,
  call = rlang::caller_env()
) {
  prepared <- if (is_motif_spec(motifs)) {
    prepare_spec_motif_args(
      glycans,
      motifs,
      alignments,
      match_degree,
      strict_sub,
      ignore_linkages,
      call = call
    )
  } else {
    prepare_regular_motif_args(
      glycans,
      motifs,
      alignments,
      match_degree,
      ignore_linkages,
      call = call
    )
  }

  if (!prepared$allow_duplicate_motifs) {
    validate_duplicate_motifs(prepared$motifs, call = call)
  }

  motif_type <- get_motif_type(prepared$motifs, call = call)
  alignments <- resolve_prepared_alignments(
    prepared$motifs,
    motif_type,
    prepared$alignments,
    prepared$match_degree,
    prepared$from_spec
  )
  glycans <- if (prepared$from_spec) {
    prepared$glycans
  } else {
    ensure_glycans_are_structures(prepared$glycans, call = call)
  }

  motifs <- ensure_motifs_are_structures(
    prepared$motifs,
    motif_type,
    require_scalar = single_motif,
    call = call
  )
  match_degree <- validate_match_degree(
    prepared$match_degree,
    motifs,
    single_motif,
    call = call
  )
  warn_mismatched_structure_levels(glycans, motifs)

  new_prepared_motif_args(
    glycans = glycans,
    motifs = motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    single_motif = single_motif
  )
}

#' Warn About Mismatched Structure Levels
#'
#' Informs users when lower-information `glycans` are matched against
#' higher-information `motifs`, a common mistake that usually returns no matches.
#'
#' @param glycans A normalized glycan structure vector.
#' @param motifs A normalized motif structure vector.
#'
#' @return `NULL`, invisibly.
#' @noRd
warn_mismatched_structure_levels <- function(glycans, motifs) {
  glycan_level <- suppressWarnings(glyrepr::get_structure_level(glycans))
  motif_level <- suppressWarnings(glyrepr::get_structure_level(motifs))

  if (
    glycan_level == "topological" &&
      motif_level %in% c("intact", "partial")
  ) {
    warn_structure_level_mismatch(glycan_level, motif_level)
  }

  if (glycan_level == "basic" && motif_level != "basic") {
    warn_structure_level_mismatch(glycan_level, motif_level)
  }

  invisible(NULL)
}

#' Emit Structure-Level Mismatch Warning
#'
#' @param glycan_level The aggregate structure level of `glycans`.
#' @param motif_level The aggregate structure level of `motifs`.
#'
#' @return `NULL`, invisibly.
#' @noRd
warn_structure_level_mismatch <- function(glycan_level, motif_level) {
  cli::cli_warn(
    c(
      "Matching lower-level {.arg glycans} against higher-level {.arg motifs} usually returns no matches.",
      "i" = "{.arg glycans} have {.val {glycan_level}} structure level, while {.arg motifs} have {.val {motif_level}} structure level.",
      "i" = "Use motifs at the same structure level as the glycans, or reduce motif structure levels before matching.",
      "i" = "See {.code ?get_structure_level} for details."
    ),
    class = "warning_mismatched_structure_level"
  )
  invisible(NULL)
}

#' Test if an Object is a Motif Specification
#'
#' @param motifs An object passed to the `motifs` argument.
#'
#' @return `TRUE` if `motifs` is a motif specification.
#' @noRd
is_motif_spec <- function(motifs) {
  inherits(motifs, "dynamic_motifs_spec") ||
    inherits(motifs, "branch_motifs_spec") ||
    inherits(motifs, "db_motifs_spec")
}

#' Resolve Motif Specification Arguments
#'
#' @inheritParams prepare_motif_args
#'
#' @return A partially normalized motif argument list.
#' @noRd
prepare_spec_motif_args <- function(
  glycans,
  motifs,
  alignments,
  match_degree,
  strict_sub,
  ignore_linkages,
  call = rlang::caller_env()
) {
  glycans <- ensure_glycans_are_structures(glycans, call = call)
  resolved <- resolve_motif_spec(
    glycans,
    motifs,
    alignments,
    match_degree,
    strict_sub,
    ignore_linkages
  )
  motifs <- name_resolved_motifs(resolved$motifs)

  list(
    glycans = glycans,
    motifs = motifs,
    alignments = resolved$alignments,
    match_degree = resolved$match_degree,
    from_spec = TRUE,
    allow_duplicate_motifs = isTRUE(resolved$allow_duplicate_motifs)
  )
}

#' Normalize Regular Motif Arguments
#'
#' @inheritParams prepare_motif_args
#'
#' @return A partially normalized motif argument list.
#' @noRd
prepare_regular_motif_args <- function(
  glycans,
  motifs,
  alignments,
  match_degree,
  ignore_linkages,
  call = rlang::caller_env()
) {
  if (is.null(match_degree)) {
    valid_alignments_arg(alignments, motifs)
  }
  valid_ignore_linkages_arg(ignore_linkages)

  list(
    glycans = glycans,
    motifs = motifs,
    alignments = alignments,
    match_degree = match_degree,
    from_spec = FALSE,
    allow_duplicate_motifs = FALSE
  )
}

#' Resolve Prepared Alignments
#'
#' @param motifs Motifs after any specification resolution.
#' @param motif_type The motif input type.
#' @param alignments User-provided or resolved alignments.
#' @param match_degree Optional degree-matching mask.
#' @param from_spec Whether the motifs came from a motif specification.
#'
#' @return A character vector of alignments.
#' @noRd
resolve_prepared_alignments <- function(
  motifs,
  motif_type,
  alignments,
  match_degree,
  from_spec
) {
  if (from_spec) {
    return(alignments)
  }
  if (!is.null(match_degree)) {
    return(vctrs::vec_recycle("substructure", length(motifs)))
  }
  decide_alignments(motifs, motif_type, alignments)
}

#' Name Resolved Motifs
#'
#' Adds IUPAC names to resolved motifs when the resolver did not already assign
#' names for downstream matrix or list labeling.
#'
#' @param motifs A resolved motif vector.
#'
#' @return `motifs`, possibly with names.
#' @noRd
name_resolved_motifs <- function(motifs) {
  if (length(motifs) > 0 && is.null(names(motifs))) {
    names(motifs) <- as.character(motifs)
  }
  motifs
}

#' Validate Duplicate Motifs
#'
#' @param motifs Motifs to check.
#' @param call The call environment used for errors.
#'
#' @return `NULL`, invisibly.
#' @noRd
validate_duplicate_motifs <- function(motifs, call = rlang::caller_env()) {
  if (!has_duplicate_motifs(motifs)) {
    return(invisible(NULL))
  }

  dupes <- unique(motifs[duplicated(motifs)])
  cli::cli_abort(
    c(
      "`motifs` cannot have duplications.",
      "x" = "Duplicate motifs: {.val {dupes}}.",
      "i" = "Consider using {.fn unique}."
    ),
    call = call
  )
}

#' Create Prepared Motif Argument List
#'
#' @inheritParams prepare_motif_args
#'
#' @return A normalized argument list for downstream motif internals.
#' @noRd
new_prepared_motif_args <- function(
  glycans,
  motifs,
  alignments,
  ignore_linkages,
  strict_sub,
  match_degree,
  single_motif
) {
  common_args <- list(
    glycans = glycans,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree
  )

  if (single_motif) {
    return(c(
      common_args["glycans"],
      list(motif = motifs, alignment = alignments),
      common_args[c("ignore_linkages", "strict_sub", "match_degree")]
    ))
  }

  c(
    common_args["glycans"],
    list(motifs = motifs, alignments = alignments),
    common_args[c("ignore_linkages", "strict_sub", "match_degree")]
  )
}

# Works with both character vectors and glyrepr_structure objects
has_duplicate_motifs <- function(motifs) {
  if (length(motifs) <= 1) {
    return(FALSE)
  }
  length(unique(motifs)) < length(motifs)
}

#' Validate Match Degree
#'
#' @param match_degree A logical vector or a list of logical vectors.
#' @param motifs A `glyrepr_structure` object.
#' @param single_motif A logical scalar indicating if `motifs` is a single motif.
#' @param call The call environment.
#'
#' @return A normalized logical vector or list of logical vectors.
#' @noRd
validate_match_degree <- function(
  match_degree,
  motifs,
  single_motif,
  call = rlang::caller_env()
) {
  if (is.null(match_degree)) {
    return(NULL)
  }

  if (single_motif) {
    if (!is.logical(match_degree)) {
      cli::cli_abort(
        "`match_degree` must be a logical vector or NULL.",
        call = call
      )
    }
    if (length(match_degree) == 0) {
      cli::cli_abort("`match_degree` cannot be empty.", call = call)
    }
    motif_nodes <- igraph::vcount(glyrepr::get_structure_graphs(motifs))
    if (length(match_degree) == 1) {
      return(rep(match_degree, motif_nodes))
    }
    if (length(match_degree) != motif_nodes) {
      cli::cli_abort(
        "`match_degree` must have length 1 or match the number of motif nodes.",
        call = call
      )
    }
    return(match_degree)
  }

  if (!is.list(match_degree)) {
    cli::cli_abort(
      "`match_degree` must be a list of logical vectors or NULL.",
      call = call
    )
  }
  if (length(match_degree) != length(motifs)) {
    cli::cli_abort(
      "`match_degree` must have the same length as `motifs`.",
      call = call
    )
  }
  if (any(vapply(match_degree, is.null, logical(1)))) {
    cli::cli_abort(
      "Each element of the `match_degree` list cannot be NULL.",
      call = call
    )
  }

  motif_graphs <- glyrepr::get_structure_graphs(motifs)
  if (inherits(motif_graphs, "igraph")) {
    motif_graphs <- list(motif_graphs)
  }
  purrr::map2(match_degree, motif_graphs, function(mask, graph) {
    if (!is.logical(mask)) {
      cli::cli_abort(
        "Each `match_degree` element must be a logical vector.",
        call = call
      )
    }
    if (length(mask) == 0) {
      cli::cli_abort(
        "Each `match_degree` element cannot be empty.",
        call = call
      )
    }
    if (length(mask) == 1) {
      return(rep(mask, igraph::vcount(graph)))
    }
    if (length(mask) != igraph::vcount(graph)) {
      cli::cli_abort(
        "Each `match_degree` element must have length 1 or match the number of motif nodes.",
        call = call
      )
    }
    mask
  })
}

# ----- Argument validation -----
valid_alignment_arg <- function(x) {
  checkmate::assert_choice(
    x,
    c("substructure", "core", "terminal", "whole"),
    null.ok = TRUE
  )
}

valid_alignments_arg <- function(alignments, motifs) {
  # Validate alignments parameter without modification
  if (!is.null(alignments)) {
    if (length(alignments) != 1 && length(alignments) != length(motifs)) {
      rlang::abort(
        "`alignments` must be NULL, a single value, or have the same length as `motifs`."
      )
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
  "or a character vector of database motif names."
)

# ----- Decide motif type -----
# Decide if the `motifs` argument is database motif names,
# an IUPAC-condensed structure character vector, or a 'glyrepr_structure' object.
get_motif_type <- function(motifs, call = rlang::caller_env()) {
  # If it is neither a glycan graph or a database motif name, it is assumed to
  # be an IUPAC-condensed structure string.
  # This assumption may not be correct, for it is possible that a wrong
  # motif name is passed in.
  if (glyrepr::is_glycan_structure(motifs)) {
    return("structure")
  } else if (is.character(motifs)) {
    known_motif_idx <- is_db_motif_name(motifs)
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
  # 1. If the motif is a database motif name: use the alignment type in the database unless provided.
  # 2. If the motif is an IUPAC string or a 'glyrepr_structure' object: use "substructure" unless provided.
  # In the first case, if the user-provided alignment is different from that in the database,
  # issue a warning.

  # Handle NA motif_type (invalid motif type)
  if (is.na(motif_type)) {
    if (is.null(alignments)) {
      alignments <- "substructure"
    }
  } else if (motif_type == "known") {
    db_alignments <- db_motif_alignment_by_name(motifs)
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
        cli::cli_abort(
          c(
            "`glycans` must be a 'glyrepr_structure' object or a glycan structure character vector.",
            "x" = "Some glycans could not be parsed as valid glycan structures."
          ),
          call = call,
          parent = cnd
        )
      }
    )
  }

  # Case 3: `glycans` has other types
  cli::cli_abort(
    c(
      "`glycans` must be a 'glyrepr_structure' object or a glycan structure character vector.",
      "x" = "The input is of class {.cls {class(glycans)}}."
    ),
    call = call
  )
}

# Make sure `motifs` is a `glyrepr_structure` object.
# This function handles both single and multiple motifs
ensure_motifs_are_structures <- function(
  motifs,
  motif_type,
  require_scalar = FALSE,
  call = rlang::caller_env()
) {
  # Check scalar requirement first if needed
  if (require_scalar && length(motifs) != 1) {
    cli::cli_abort(
      c(
        "The `motif` argument must be a scalar vector.",
        "x" = "The input is of length {.val {length(motifs)}}."
      ),
      call = call
    )
  }

  # Determine which error message to use based on scalar requirement
  base_err_msg <- if (require_scalar) {
    paste(
      "`motif` must be either a 'glyrepr_structure' object with length 1,",
      "a glycan structure character scalar,",
      "or a database motif name."
    )
  } else {
    paste0(
      "`motifs` must be a 'glyrepr_structure' object,",
      "a character vector of glycan structure strings,",
      "or a character vector of database motif names."
    )
  }

  # Case 1: `motifs` is of other types
  if (is.na(motif_type)) {
    cli::cli_abort(
      c(
        base_err_msg,
        "x" = "The input is of class {.cls {class(motifs)}}."
      ),
      call = call
    )
  }

  # Case 2: `motifs` is already a `glyrepr_structure` object
  if (motif_type == "structure") {
    return(motifs)
  }

  # Case 3: `motifs` is a vector of database motif names
  if (motif_type == "known") {
    tryCatch(
      return(db_motif_structure_by_name(motifs)),
      error = function(cnd) {
        cli::cli_abort(
          c(
            base_err_msg,
            "i" = "Use {.fn db_motif_info} to see all valid motif names."
          ),
          call = call,
          parent = cnd
        )
      }
    )
  }

  # Case 4: `motifs` is a vector of glycan structure strings
  if (motif_type == "iupac") {
    tryCatch(
      return(glyparse::auto_parse(motifs)),
      error = function(cnd) {
        cli::cli_abort(
          c(
            base_err_msg,
            "x" = "Some motifs are neither valid glycan structures nor database motif names.",
            "i" = "Use {.fn db_motif_info} to see all valid motif names."
          ),
          call = call,
          parent = cnd
        )
      }
    )
  }
}

# Simplified wrapper function for single motifs
ensure_motif_is_structure <- function(
  motif,
  motif_type,
  call = rlang::caller_env()
) {
  ensure_motifs_are_structures(
    motif,
    motif_type,
    require_scalar = TRUE,
    call = call
  )
}

# ----- Generic function for single motif mapping -----
apply_single_motif_to_glycans <- function(
  glycans,
  motif,
  alignment,
  ignore_linkages,
  strict_sub,
  match_degree,
  single_glycan_func,
  smap_func
) {
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
  motif_has_linkages <- glyrepr::has_linkages(motif)[[1]]
  motif_composition_profile <- new_motif_composition_profile(motif_graph)
  smap_func(
    glycans_to_use,
    single_glycan_func,
    motif_graph,
    motif_has_linkages,
    motif_composition_profile,
    alignment,
    ignore_linkages,
    strict_sub,
    match_degree
  )
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
  } else {
    # must be a character vector
    if (is.null(names(x))) {
      return(x)
    } else {
      return(names(x))
    }
  }
}

# Helper to determine motif names for matrix columns
# Returns: explicit names if present, database motif names if applicable,
#          NULL otherwise (IUPAC strings or structure input without names)
prepare_motif_names <- function(motifs_input) {
  # Handle motif spec objects
  # These should not use their internal list names as motif names
  if (
    inherits(motifs_input, "dynamic_motifs_spec") ||
      inherits(motifs_input, "branch_motifs_spec") ||
      inherits(motifs_input, "db_motifs_spec")
  ) {
    return(NULL)
  }

  # If motifs_input has explicit names, use them
  if (!is.null(names(motifs_input))) {
    return(names(motifs_input))
  }

  # If motifs_input is a character vector of database motif names, use those
  if (is.character(motifs_input) && all(is_db_motif_name(motifs_input))) {
    return(motifs_input)
  }

  # Otherwise, no names (IUPAC strings or structure input without names)
  NULL
}

apply_motifs_to_glycans <- function(
  glycans,
  motifs,
  alignments,
  ignore_linkages,
  single_motif_func,
  glycan_names,
  motif_names,
  strict_sub,
  match_degree
) {
  # Generic function to apply multiple motifs to multiple glycans
  # single_motif_func should be either have_motif_ or count_motif_

  # Handle empty motifs case
  if (length(motifs) == 0) {
    cli::cli_abort("`motifs` cannot be empty.")
  }

  match_degree_list <- if (is.null(match_degree)) {
    rep(list(NULL), length(motifs))
  } else {
    match_degree
  }

  # Apply each motif to all glycans using purrr
  motif_results_list <- purrr::map2(
    motifs,
    seq_along(motifs),
    ~ single_motif_func(
      glycans,
      .x,
      alignment = alignments[[.y]],
      ignore_linkages = ignore_linkages,
      strict_sub = strict_sub,
      match_degree = match_degree_list[[.y]]
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
