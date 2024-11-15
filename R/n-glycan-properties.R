# ----- Interface -----
describe_n_glycans <- function(glycans, strict = FALSE) {
  # Validate the input
  valid_glycans_arg(glycans)
  checkmate::assert_flag(strict)
  glycans <- ensure_glycans_are_graphs(glycans)

  # Deal with strictness
  if (strict) {
    hf <- has_motif_
    cf <- counts_motif_
  } else {
    glycans <- purrr::map(
      glycans, glyrepr::convert_glycan_mono_type,
      to = "simple", strict = FALSE
    )
    hf <- lenient_has_motif_
    cf <- lenient_count_motif_
  }

  # Check if the glycans are N-glycans
  invalid_indices <- which(!purrr::map_lgl(glycans, .is_n_glycan, hf, cf))
  if (length(invalid_indices) > 0) {
    cli::cli_abort("Glycans at indices {.val {invalid_indices}} are not N-glycans.")
  }

  # Get the properties
  res <- tibble::tibble(
    glycan_type = purrr::map_chr(glycans, .n_glycan_type, hf, cf),
    bisecting = purrr::map_lgl(glycans, .has_bisecting, hf, cf),
    antennae = purrr::map_int(glycans, .n_antennae, hf, cf),
    core_fuc = purrr::map_int(glycans, .n_core_fuc, hf, cf),
    arm_fuc = purrr::map_int(glycans, .n_arm_fuc, hf, cf),
    terminal_gal = purrr::map_int(glycans, .n_terminal_gal, hf, cf)
  )

  # Add the glycan name column
  if (!is.null(names(glycans))) {
    res <- tibble::add_column(res, glycan = names(glycans), .before = 1)
  } else if (is.character(glycans)) {
    res <- tibble::add_column(res, glycan = glycans, .before = 1)
  }

  res
}


#' Check if a Glycan is an N-Glycan
#'
#' This function checks if a glycan has the N-glycan core motif:
#' ```
#' Man
#'    \
#'     Man - GlcNAc - GlcNAc
#'    /
#' Man
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return A logical value.
#' @export
is_n_glycan <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .is_n_glycan, check_n_glycan = FALSE)
}


#' Determine N-Glycan Type
#'
#' Four types of N-glycans are recognized: high mannose, hybrid, complex, and paucimannose.
#' For more information about N-glycan types,
#' see [Essentials of Glycobiology](https://www.ncbi.nlm.nih.gov/books/NBK579964/#_s9_2_).
#'
#' @param glycan A `glycan_graph` object, or a character string of IUPAC condensed format.
#' @param strict A logical value. If `TRUE`, the glycan must have "concrete"
#' monosaccharides (e.g. "GlcNAc", "Man", "Gal") and linkage information.
#' If `FALSE`, the function is more lenient,
#' checking monosaacharide identities on the "simple" level (e.g. "H", "N", "F")
#' and ignoring linkage information.
#' Default is `FALSE`. This is preferred because in most cases the
#' structural resolution could not be high, but we known for sure the glycans are indeed N-glycans.
#'
#' @return A character string of the N-glycan type,
#' either "highmannose", "hybrid", "complex", or "paucimannose".
#' If the glycan seems not to be any of the four types, an error is thrown
#' with a message "Not an N-glycan".
#' This doesn't necessarily mean the glycan is not an N-glycan.
#' Maybe you have used the strict mode with a glycan that is not well resolved.
#'
#' @export
n_glycan_type <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_glycan_type)
}


#' Does the Glycan have Bisecting GlcNAc?
#'
#' Bisecting GlcNAc is a GlcNAc residue attached to the core mannose of N-glycans.
#' ```
#'      Man
#'         \
#' GlcNAc - Man - GlcNAc - GlcNAc
#' ~~~~~~  /
#'      Man
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return A logical value.
#' @export
has_bisecting <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .has_bisecting)
}


#' Number of Antennae
#'
#' The number of antennae is the number of branching GlcNAc to the core mannoses
#' in a complex N-glycan. Bisecting GlcNAc is not counted as an antenna.
#' This functions returns NA_integer_ for non-complex N-glycans.
#'
#' @inheritParams n_glycan_type
#'
#' @return An integer of the number of antennae.
#' @export
n_antennae <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_antennae)
}


#' Number of Core Fucoses
#'
#' Core fucoses are those fucose residues attached to the core GlcNAc of an N-glycan.
#' ```
#' Man             Fuc  <- core fucose
#'    \             |
#'     Man - GlcNAc - GlcNAc
#'    /
#' Man
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return An integer of the number of core fucoses.
#' @export
n_core_fuc <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_core_fuc)
}


#' Number of Arm Focuses
#'
#' Arm focuses are thoses focuse residues attached to the branching GlcNAc
#' of an N-glycan.
#' ```
#'  Fuc  <- arm fucose
#'   |
#' GlcNAc - Man
#'             \
#'              Man - GlcNAc - GlcNAc
#'             /
#' GlcNAc - Man
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return An integer of the number of arm fucoses.
#' @export
n_arm_fuc <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_arm_fuc)
}


#' Number of Terminal Galactoses
#'
#' Terminal galactoses are those galactose residues on the non-reducing end
#' without sialic acid capping.
#' ```
#'          Gal - GlcNAc - Man
#'          ~~~               \
#'      terminal Gal           Man - GlcNAc - GlcNAc
#'                            /
#' Neu5Ac - Gal - GlcNAc - Man
#'          ~~~
#'    not terminal Gal
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return An integer of the number of terminal galactoses.
#' @export
n_terminal_gal <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_terminal_gal)
}


# ----- Implementation -----
# All functions in this section have three arguments:
# - `glycan`: A `glycan_graph` object. It is the caller's responsibility to ensure this.
# - `.has_motif`: A function that checks if a glycan has a motif.
#   It should have the signature `.has_motif(glycan, motif, alignment)`.
# - `.counts_motif`: A function that counts the number of motifs in a glycan.
#   It should have the signature `.counts_motif(glycan, motif, alignment)`.
#
# In these functions, use `.has_motif` and `.counts_motif` to check for motifs.

.is_n_glycan <- function(glycan, .has_motif, .counts_motif) {
  core_graph <- get_motif_graph("N-Glycan core basic")
  .has_motif(glycan, core_graph, alignment = "core")
}


.n_glycan_type <- function(glycan, .has_motif, .counts_motif) {
  H4N2_iupac <- "Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  H4N2_graph <- glyparse::parse_iupac_condensed(H4N2_iupac)
  core_graph <- get_motif_graph("N-Glycan core basic")
  hybrid_graph <- get_motif_graph("N-Glycan hybrid")
  highman_graph <- get_motif_graph("N-Glycan high mannose")
  if (.has_motif(glycan, core_graph, alignment = "whole") ||
      .has_motif(glycan, H4N2_graph, alignment = "whole")) {
    "paucimannose"
  } else if (.has_motif(glycan, hybrid_graph, alignment = "core")) {
    "hybrid"
  } else if (.has_motif(glycan, highman_graph, alignment = "core")) {
    "highmannose"
  } else {
    "complex"
  }
}


.has_bisecting <- function(glycan, .has_motif, .counts_motif) {
  bisect_graph <- get_motif_graph("N-glycan core, bisected")
  .has_motif(glycan, bisect_graph, alignment = "core")
}


.n_antennae <- function(glycan, .has_motif, .counts_motif) {
  ant2_graph <- get_motif_graph("N-Glycan biantennary")
  ant3_graph <- get_motif_graph("N-Glycan triantennary")
  ant4_graph <- get_motif_graph("N-Glycan tetraantennary")
  if (.n_glycan_type(glycan, .has_motif, .counts_motif) != "complex") {
    NA_integer_
  } else if (.has_motif(glycan, ant4_graph, alignment = "core")) {
    4L
  } else if (.has_motif(glycan, ant3_graph, alignment = "core")) {
    3L
  } else if (.has_motif(glycan, ant2_graph, alignment = "core")) {
    2L
  } else {
    1L
  }
}


.n_core_fuc <- function(glycan, .has_motif, .counts_motif) {
  core_fuc_graph <- get_motif_graph("N-Glycan core, core-fucosylated")
  .counts_motif(glycan, core_fuc_graph, alignment = "core")
}


.n_arm_fuc <- function(glycan, .has_motif, .counts_motif) {
  arm_fuc_graph <- get_motif_graph("N-Glycan core, arm-fucosylated")
  .counts_motif(glycan, arm_fuc_graph, alignment = "core")
}


.n_terminal_gal <- function(glycan, .has_motif, .counts_motif) {
  terminal_gal_graph <- glyparse::parse_iupac_condensed("Gal(?1-")
  .counts_motif(glycan, terminal_gal_graph, alignment = "terminal")
}


# ----- Utilities -----
n_glycan_property_wrapper <- function(glycan, strict, func, check_n_glycan = TRUE) {
  # This function encapsulates the common pattern of checking N-glycan properties.
  # To use it, write a function that takes a glycan graph,
  # a `.has_motif` argument taking a function that checks for a motif,
  # and a `.counts_motif` argument taking a function that counts motifs.
  # Inside the function, use `.has_motif` and `.count_motif` to check for motifs,
  # instead of calling `has_motif_()` and `counts_motif_()` directly.
  #
  # For example:
  # ```
  # .has_bisecting <- function(glycan, .has_motif, .counts_motif) {
  #   bisect_graph <- get_motif_graph("N-glycan core, bisected")
  #   .has_motif(glycan, bisect_graph, alignment = "core")
  # }
  # ```
  # Then pass this function to `n_glycan_property_wrapper()`:
  # ```
  # n_glycan_property_wrapper(glycan, strict, .has_bisecting)
  # ```
  # This function will take care of converting the glycan to a graph and
  # checking the motif with suitable strictness.
  valid_glycan_arg(glycan)
  checkmate::assert_flag(strict)
  glycan <- ensure_glycan_is_graph(glycan)
  if (strict) {
    has_motif_func <- has_motif_
    counts_motif_func <- counts_motif_
  } else {
    glycan <- glyrepr::convert_glycan_mono_type(glycan, to = "simple", strict = FALSE)
    has_motif_func <- lenient_has_motif_
    counts_motif_func <- lenient_count_motif_
  }
  if (check_n_glycan && !.is_n_glycan(glycan, has_motif_func, counts_motif_func)) {
    rlang::abort("Not an N-glycan.")
  }
  func(glycan, has_motif_func, counts_motif_func)
}


lenient_has_motif_ <- function(glycan, motif, alignment) {
  # This function differs from `has_motif_()` in these ways:
  # 1. It always converts the motif to "simple" monosaccharides.
  # 2. It ignores linkage information.
  motif <- glyrepr::convert_glycan_mono_type(motif, to = "simple", strict = FALSE)
  has_motif_(glycan, motif, alignment = alignment, ignore_linkages = TRUE)
}


lenient_count_motif_ <- function(glycan, motif, alignment) {
  motif <- glyrepr::convert_glycan_mono_type(motif, to = "simple", strict = FALSE)
  counts_motif_(glycan, motif, alignment = alignment, ignore_linkages = TRUE)
}
