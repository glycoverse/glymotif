# ----- Prepare arguments for `count_motif()` and `have_motif()` -----
prepare_have_motif_args <- function(glycans, motif, alignment, ignore_linkages) {
  # Check input arguments
  valid_glycans_arg(glycans)
  valid_motif_arg(motif)
  valid_alignment_arg(alignment)
  valid_ignore_linkages_arg(ignore_linkages)

  glycan_names <- get_structure_names(glycans)

  motif_type <- get_motif_type(motif)

  # Deal with `alignment` when `motif` is a known motif name
  if (motif_type == "known") {
    alignment <- decide_alignment(motif, alignment)
  } else if (is.null(alignment)) {
    alignment <- "substructure"
  }

  # Make sure `glycan` and `motif` are structures
  glycans <- ensure_glycans_are_structures(glycans)
  motif <- ensure_motif_is_structure(motif, motif_type)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycans_mono_type()`.
  glycans <- ensure_glycans_mono_type(glycans, motif)

  list(glycans = glycans, motif = motif, alignment = alignment, ignore_linkages = ignore_linkages)
}


# ----- Prepare arguments for `have_motifs()` and `count_motifs()` -----
prepare_have_motifs_args <- function(glycans, motifs, alignments, ignore_linkages) {
  valid_glycans_arg(glycans)
  valid_motifs_arg(motifs)
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  glycan_names <- get_structure_names(glycans)
  motif_names <- get_structure_names(motifs)

  motif_type <- get_motif_type(motifs)

  # Deal with `alignments` when `motif` is a known motif name
  if (motif_type == "known") {
    # For multiple known motifs, handle alignments for each
    if (is.null(alignments)) {
      alignments <- purrr::map_chr(motifs, get_motif_alignment)
    } else {
      # Handle alignment parameter
      if (length(alignments) == 1) {
        alignments <- rep(alignments, length(motifs))
      }
      # Check each alignment against database
      db_alignments <- purrr::map_chr(motifs, get_motif_alignment)
      for (i in seq_along(alignments)) {
        if (alignments[i] != db_alignments[i]) {
          cli::cli_warn(
            "The provided alignment type {.val {alignments[i]}} is different from the motif's alignment type {.val {db_alignments[i]}} in database for motif {.val {motifs[i]}}.",
            class = "warning_custom_alignment"
          )
        }
      }
    }
  } else if (is.null(alignments)) {
    alignments <- rep("substructure", length(motifs))
  } else {
    # Handle alignment parameter
    if (length(alignments) == 1) {
      alignments <- rep(alignments, length(motifs))
    }
  }

  # Make sure `glycan` and `motif` are structures
  glycans <- ensure_glycans_are_structures(glycans)
  motifs <- ensure_motifs_are_structures(motifs, motif_type)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycans_mono_type()`.
  glycans <- ensure_glycans_mono_type(glycans, motifs[[1]])

  # Restore names after conversion
  names(glycans) <- glycan_names
  names(motifs) <- motif_names

  list(glycans = glycans, motifs = motifs, alignments = alignments, ignore_linkages = ignore_linkages)
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


valid_glycans_arg <- function(x) {
  # Must be a structure or a character vector
  if (!glyrepr::is_glycan_structure(x) && !is.character(x)) {
    rlang::abort(glycan_type_err_msg)
  }
}


valid_motif_arg <- function(x) {
  if (
    (!glyrepr::is_glycan_structure(x) && length(x) != 1) ||
    (is.character(x) && length(x) != 1)
  ) {
    rlang::abort(motif_type_err_msg)
  }
}


valid_motifs_arg <- function(x) {
  if (length(x) == 0) {
    rlang::abort("`motifs` cannot be empty.")
  }
  if ((!glyrepr::is_glycan_structure(x)) && (!is.character(x))) {
    rlang::abort(motifs_type_err_msg)
  }
}


valid_ignore_linkages_arg <- function(x) {
  checkmate::assert_flag(x)
}


# ----- Get names -----
get_structure_names <- function(x) {
  if (is.null(names(x))) {
    as.character(x)
  } else {
    names(x)
  }
}


# ----- Decide motif type -----
get_motif_type <- function(motif) {
  # Decide if the `motif` argument is a known motif,
  # an IUPAC-condensed structure string, or a 'glyrepr_structure' object.
  # If it is neither a glycan graph or a known motif, it is assumed to
  # be an IUPAC-condensed structure string.
  # This assumption may not be correct, for it is possible that a wrong
  # motif name is passed in.
  if (glyrepr::is_glycan_structure(motif)) {
    return("structure")
  } else if (is.character(motif)) {
    if (all(is_known_motif(motif))) {
      return("known")
    } else if (any(is_known_motif(motif))) {
      rlang::abort("Some motifs are known, but some are not.")
    } else {
      return("iupac")
    }
  } else {
    rlang::abort(motif_type_err_msg)
  }
}


# ----- Decide alignment type -----
decide_alignment <- function(motif_name, alignment) {
  # Decide which alignment type to use.
  # If the alignment is provided, check if it is the same as the motif's
  # alignment type in the database and issue a warning if not.
  # If not provided, use the motif's alignment type in the database.
  db_alignment <- get_motif_alignment(motif_name)
  if (!is.null(alignment)) {
    if (alignment != db_alignment) {
      cli::cli_warn(
        "The provided alignment type {.val {alignment}} is different from the motif's alignment type {.val {db_alignment}} in database.",
        class = "warning_custom_alignment"
      )
    }
  } else {
    alignment <- db_alignment
  }
  alignment
}

# ----- Make sure glycans and motifs are structures -----
ensure_glycans_are_structures <- function(glycans) {
  # Make sure `glycans` is a `glyrepr_structure` object.
  if (is.character(glycans)) {
    tryCatch(
      glycans <- glyparse::parse_iupac_condensed(glycans),
      error = function(e) {
        rlang::abort("`glycan` could not be parsed as a valid IUPAC-condensed structure.")
      }
    )
  }
  glycans
}


ensure_motif_is_structure <- function(motif, motif_type) {
  # Make sure `motif` is a structure.
  if (motif_type == "known") {
    motif <- get_motif_structure(motif)
  } else if (motif_type == "iupac") {
    tryCatch(
      motif <- glyparse::parse_iupac_condensed(motif),
      error = function(e) {
        rlang::abort(motif_type_err_msg)
      }
    )
  } else {  # motif_type == "structure"
    motif <- motif
  }

  motif
}


ensure_motifs_are_structures <- function(motifs, motif_type) {
  if (motif_type == "known") {
    motifs <- get_motif_structure(motifs)
  } else if (motif_type == "iupac") {
    tryCatch(
      motifs <- glyparse::parse_iupac_condensed(motifs),
      error = function(e) {
        rlang::abort(motif_type_err_msg)
      }
    )
  } else {  # motif_type == "structure"
    motifs <- motifs
  }
  motifs
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
