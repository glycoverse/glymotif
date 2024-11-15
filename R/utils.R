# ----- Prepare arguments for `counts_motif()` and `has_motif()` -----
prepare_has_motif_args <- function(glycan, motif, alignment, ignore_linkages) {
  # Check input arguments
  valid_glycan_arg(glycan)
  valid_motif_arg(motif)
  valid_alignment_arg(alignment)
  valid_ignore_linkages_arg(ignore_linkages)

  motif_type <- get_motif_type(motif)

  # Deal with `alignment` when `motif` is a known motif name
  if (motif_type == "known") {
    alignment <- decide_alignment(motif, alignment)
  } else if (is.null(alignment)) {
    alignment <- "substructure"
  }

  # Make sure `glycan` and `motif` are graphs
  glycan <- ensure_glycan_is_graph(glycan)
  motif <- ensure_motif_is_graph(motif, motif_type)

  # Ensure that `glycan` and `motif` are all "NE" type
  glycan <- ensure_ne_graph(glycan)
  motif <- ensure_ne_graph(motif)

  # Ensure that `glycan` and `motif` have the same monosaccharide type
  # To ensure strict comparison, if the glycan type is lower than the motif type,
  # an error will be raised by `ensure_glycan_mono_type()`.
  glycan <- ensure_glycan_mono_type(glycan, motif)

  list(glycan = glycan, motif = motif, alignment = alignment, ignore_linkages = ignore_linkages)
}


# ----- Argument validation -----
valid_alignment_arg <- function(x) {
  checkmate::assert(
    checkmate::test_null(x),
    checkmate::test_choice(x, c("substructure", "core", "terminal", "whole"))
  )
}


test_graph_arg <- function(x) {
  (
    checkmate::test_class(x, "glycan_graph") ||
      checkmate::test_character(x, len = 1)
  )
}


glycan_type_err_msg <- paste(
  "`glycan` must be a 'glycan_graph' object",
  "or an IUPAC-condensed structure string."
)


motif_type_err_msg <- paste(
  "`motif` must be either a 'glycan_graph' object,",
  "an IUPAC-condensed structure string,",
  "or a known motif name."
)


valid_glycan_arg <- function(x) {
  if (!test_graph_arg(x)) {
    rlang::abort(glycan_type_err_msg)
  }
}


valid_motif_arg <- function(x) {
  if (!test_graph_arg(x)) {
    rlang::abort(motif_type_err_msg)
  }
}


valid_ignore_linkages_arg <- function(x) {
  checkmate::assert_flag(x)
}


valid_glycans_arg <- function(x) {
  error_msg <- paste(
    "`glycans` must be either a list of 'glycan_graph' objects",
    "or a character vector of IUPAC-condensed structure strings."
  )
  if (
    !checkmate::test_character(x) &&
    !checkmate::test_list(x, types = "glycan_graph")
  ) {
    rlang::abort(error_msg)
  }
}


valid_motifs_arg <- function(x) {
  error_msg <- paste(
    "`motifs` must be either a list of 'glycan_graph' objects,",
    "a character vector of IUPAC-condensed structure strings,",
    "or a character vector of known motif names."
  )
  if (
    !checkmate::test_null(x) &&
    !checkmate::test_character(x) &&
    !checkmate::test_list(x, types = "glycan_graph")
  ) {
    rlang::abort(error_msg)
  }
}


valid_alignments_arg <- function(x, motifs) {
  checkmate::assert(
    checkmate::test_null(x),
    checkmate::test_subset(x, c("substructure", "core", "terminal", "whole"))
  )
  if (!is.null(x) && !length(x) %in% c(1, length(motifs))) {
    cli::cli_abort(c(
      "`alignments` must be either a single character string or a character vector of the same length as `motifs`.",
      i = "`motif` length: {.val {length(motifs)}}, `alignments` length: {.val {length(x)}}"
    ))
  }
}


# ----- Decide motif type -----
get_motif_type <- function(motif) {
  # Decide if the `motif` argument is a known motif,
  # an IUPAC-condensed structure string, or a 'glycan_graph' object.
  # If it is neither a glycan graph or a known motif, it is assumed to
  # be an IUPAC-condensed structure string.
  # This assumption may not be correct, for it is possible that a wrong
  # motif name is passed in.
  if (glyrepr::is_glycan(motif)) {
    return("glycan_graph")
  } else if (is.character(motif)) {
    if (is_known_motif(motif)) {
      return("known")
    } else {
      return("iupac")
    }
  } else {
    rlang::abort(motif_type_err_msg)
  }
}


get_motifs_type <- function(motifs) {
  if (purrr::every(motifs, glyrepr::is_glycan)) {
    return("glycan_graph")
  } else if (all(is.character(motifs))) {
    if (purrr::every(motifs, is_known_motif)) {
      return("known")
    } else if (purrr::some(motifs, is_known_motif)) {
      unknown_motifs <- motifs[!purrr::map_lgl(motifs, is_known_motif)]
      cli::cli_abort("Motifs {.val {unknown_motifs}} are not known motif names.")
    } else {
      # If all motifs are characters but not known motif names,
      # assume they are IUPAC-condensed structure strings.
      return("iupac")
    }
  } else {
    rlang::abort("`motifs` must be either 'glycan_graph' objects, a character vector of IUPAC-condensed structure strings, or a character vector of known motif names.")
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


decide_alignments <- function(motif_names, alignments) {
  db_alignments <- get_motif_alignment(motif_names)
  if (!is.null(alignments)) {
    rlang::warn("Use user-provided alignments, not the ones in the database.")
    alignments
  } else {
    db_alignments
  }
}


# ----- Make sure glycans and motifs are graphs -----
ensure_glycan_is_graph <- function(glycan) {
  # Make sure `glycan` is a graph.
  if (is.character(glycan)) {
    tryCatch(
      glycan <- glyparse::parse_iupac_condensed(glycan),
      error = function(e) {
        rlang::abort("`glycan` could not be parsed as a valid IUPAC-condensed structure.")
      }
    )
  }
  glycan
}


ensure_motif_is_graph <- function(motif, motif_type) {
  # Make sure `motif` is a graph.
  if (motif_type == "known") {
    motif <- get_motif_graph(motif)
  } else if (motif_type == "iupac") {
    tryCatch(
      motif <- glyparse::parse_iupac_condensed(motif),
      error = function(e) {
        rlang::abort(motif_type_err_msg)
      }
    )
  } else {  # motif_type == "glycan_graph"
    motif <- motif
  }

  motif
}


ensure_motifs_are_graphs <- function(motifs, motif_type) {
  # Convert to graphs
  if (motif_type == "glycan_graph") {
    graph_list <- motifs
  } else if (motif_type == "iupac") {
    graph_list <- purrr::map(motifs, try_parse_iupac_condensed)
    if (any(is.na(graph_list))) {
      invalid_indices <- which(is.na(graph_list))
      msg <- paste(
        "Motifs at indices {.val {invalid_indices}} are neither known motif",
        "names or able to be parsed as IUPAC-condensed structure strings."
      )
      cli::cli_abort(msg)
    }
  } else if (motif_type == "known") {
    graph_list <- purrr::map(motifs, ~ get_motif_graph(.x))
  }

  # Add names
  if (is.null(names(motifs)) && is.character(motifs)) {
    names(graph_list) <- motifs
  }
  graph_list
}


ensure_glycans_are_graphs <- function(glycans) {
  if(is.character(glycans)) {
    graph_list <- purrr::map(glycans, try_parse_iupac_condensed)
    if (any(is.na(graph_list))) {
      invalid_indices <- which(is.na(graph_list))
      msg <- paste(
        "Glycans at indices {.val {invalid_indices}} are not able to be",
        "parsed as IUPAC-condensed structure strings."
      )
      cli::cli_abort(msg)
    }
    if (is.null(names(graph_list))) {
      names(graph_list) <- glycans
    }
  } else {
    graph_list <- glycans
  }
  graph_list
}


# ----- Ensure NE graphs -----
ensure_ne_graph <- function(graph) {
  glyrepr::convert_graph_mode(graph, to = "ne", strict = FALSE)
}


ensure_ne_graphs <- function(graphs) {
  purrr::map(graphs, ensure_ne_graph)
}


# ----- Check and align mono types -----
ensure_glycan_mono_type <- function(glycan, motif) {
  # Make the glycan the same monosaccharide type (concrete, generic or simple) as the motif.
  # If the glycan type is lower than the motif type, an error will be raised.
  glycan_type <- glyrepr::decide_glycan_mono_type(glycan)
  motif_type <- glyrepr::decide_glycan_mono_type(motif)

  conditions <- (
    (glycan_type == "concrete" && motif_type != "concrete") ||
      (glycan_type == "generic" && motif_type == "simple")
  )
  if (conditions) return(glyrepr::convert_glycan_mono_type(glycan, motif_type))

  if (glycan_type != motif_type) {
    cli::cli_abort(c(
      "The monosaccharide type of `glycan` cannot be obscurer than `motif`.",
      "x" = "{.val {glycan_type}} is obscurer than {.val {motif_type}}."
    ))
  }
  return(glycan)
}


ensure_glycans_mono_type <- function(glycans, motif) {
  # Make sure all glycans have the same monosaccharide type as the motif.
  # This function assumes that all glycans have the same monosaccharide type.
  glycan_type <- glyrepr::decide_glycan_mono_type(glycans[[1]])
  motif_type <- glyrepr::decide_glycan_mono_type(motif)

  conditions <- (
    (glycan_type == "concrete" && motif_type != "concrete") ||
      (glycan_type == "generic" && motif_type == "simple")
  )
  if (conditions) {
    return(purrr::map(glycans, glyrepr::convert_glycan_mono_type, to = motif_type))
  }

  if (glycan_type != motif_type) {
    cli::cli_abort(c(
      "The monosaccharide type of `glycans` cannot be obscurer than `motif`.",
      "x" = "{.val {glycan_type}} is obscurer than {.val {motif_type}}."
    ))
  }
  return(glycans)
}


same_mono_types <- function(graphs) {
  # Check if all graphs have the same monosaccharide type.
  mono_types <- purrr::map(graphs, glyrepr::decide_glycan_mono_type)
  dplyr::n_distinct(mono_types) <= 1
}
