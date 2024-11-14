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
    } else if (.has_motif(glycan, core_graph, alignment = "core")){
      "complex"
    } else {
      rlang::abort("Not an N-glycan")
    }
  }
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
  .has_bisecting <- function(glycan, .has_motif, .counts_motif) {
    bisect_graph <- get_motif_graph("N-glycan core, bisected")
    .has_motif(glycan, bisect_graph, alignment = "core")
  }
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
  .n_antennae <- function(glycan, .has_motif, .counts_motif) {
    ant2_graph <- get_motif_graph("N-Glycan biantennary")
    ant3_graph <- get_motif_graph("N-Glycan triantennary")
    ant4_graph <- get_motif_graph("N-Glycan tetraantennary")
    if (n_glycan_type(glycan, strict) != "complex") {
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
  n_glycan_property_wrapper(glycan, strict, .n_antennae)
}


#' Number of Core Fucoses
#'
#' The number of core fucoses is the number of fucose residues attached to the
#' core GlcNAc of an N-glycan.
#' ```
#' Man           Fuc  <- core fucose
#'    \           |
#'     GlcNAc - GlcNAc
#'    /
#' Man
#' ```
#'
#' @inheritParams n_glycan_type
#'
#' @return An integer of the number of core fucoses.
#' @export
n_core_fuc <- function(glycan, strict = FALSE) {
  .n_core_func <- function(glycan, .has_motif, .counts_motif) {
    core_fuc_graph <- get_motif_graph("N-Glycan core, core-fucosylated")
    .counts_motif(glycan, core_fuc_graph, alignment = "core")
  }
  n_glycan_property_wrapper(glycan, strict, .n_core_func)
}


n_glycan_property_wrapper <- function(glycan, strict, func) {
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
