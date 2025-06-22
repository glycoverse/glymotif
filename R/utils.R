# ----- Prepare arguments -----
prepare_have_motifs_args <- function(glycans, motifs, alignments, ignore_linkages) {
  # for `have_motifs()`, `count_motifs()`
  valid_glycans_arg(glycans)
  valid_motifs_arg(motifs)
  valid_alignments_arg(alignments, motifs)
  valid_ignore_linkages_arg(ignore_linkages)

  motif_type <- get_motif_type(motifs)
  alignments <- decide_alignments(motifs, motif_type, alignments)
  
  # Make sure `glycan` and `motif` are structures
  glycans <- ensure_glycans_are_structures(glycans)
  motifs <- ensure_motifs_are_structures(motifs, motif_type)

  # Note: We no longer enforce global mono type consistency here.
  # Each motif will handle its own type conversion in apply_single_motif_to_glycans.

  list(glycans = glycans, motifs = motifs, alignments = alignments, ignore_linkages = ignore_linkages)
}

prepare_have_motif_args <- function(glycans, motif, alignment, ignore_linkages) {
  # for `have_motif()`, `count_motif()`
  if (length(motif) != 1) {
    rlang::abort("`motif` must be a single structure.")
  }
  params <- prepare_have_motifs_args(glycans, motif, alignment, ignore_linkages)
  params$motif <- params$motifs
  params$alignment <- params$alignments
  params$motifs <- NULL
  params$alignments <- NULL
  params
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
decide_alignments <- function(motifs, motif_type, alignments) {
  # Decide which alignment type to use.
  # 1. If the motif is a known motif name: use the alignment type in the database unless provided.
  # 2. If the motif is an IUPAC string or a 'glyrepr_structure' object: use "substructure" unless provided.
  # In the first case, if the user-provided alignment is different from that in the database,
  # issue a warning.
  if (motif_type == "known") {
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
ensure_glycans_are_structures <- function(glycans) {
  # Make sure `glycans` is a `glyrepr_structure` object.
  # Fixed: Check if it's character AND not already a glyrepr_structure
  # This prevents unnecessary parsing of already-parsed structures
  if (is.character(glycans) && !glyrepr::is_glycan_structure(glycans)) {
    tryCatch(
      glycans <- glyparse::parse_iupac_condensed(glycans),
      error = function(e) {
        rlang::abort("`glycan` could not be parsed as a valid IUPAC-condensed structure.")
      }
    )
  }
  glycans
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
