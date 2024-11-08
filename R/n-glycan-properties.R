ng_motifs <- c(
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     └─Man (a1-3)
  pman1 = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3/6)
  #     └─Man (a1-3)
  pman2 = "Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-3)
  #     └─Man (a1-6)
  #       ├─Man (a1-6)
  #       └─Man (a1-3)
  hman = "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  hybrid = "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # same as pman1
  core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     ├─GlcNAc (b1-4)
  #     └─Man (a1-3)
  bisecting_core = "Man(a1-3)[GlcNAc(b1-4)][Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
)
ng_motifs <- purrr::map(ng_motifs, glyparse::parse_iupac_condensed)
ng_s_motifs <- purrr::map(ng_motifs, glyrepr::convert_glycan_mono_type, to = "simple")


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
  valid_glycan_arg(glycan)
  glycan <- ensure_glycan_is_graph(glycan)
  if (strict) {
    motifs <- ng_motifs
  } else {
    motifs <- ng_s_motifs
    glycan <- glyrepr::convert_glycan_mono_type(glycan, to = "simple", strict = FALSE)
  }
  has_motif_ <- purrr::partial(has_motif_, ignore_linkages = !strict)

  if (has_motif_(glycan, motifs$pman1, alignment = "whole") ||
      has_motif_(glycan, motifs$pman2, alignment = "whole")) {
    "paucimannose"
  } else if (has_motif_(glycan, motifs$hybrid, alignment = "core")) {
    "hybrid"
  } else if (has_motif_(glycan, motifs$hman, alignment = "core")) {
    "highmannose"
  } else if (has_motif_(glycan, motifs$core, alignment = "core")){
    "complex"
  } else {
    rlang::abort("Not an N-glycan")
  }
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
  valid_glycan_arg(glycan)
  glycan <- ensure_glycan_is_graph(glycan)
  if (strict) {
    motifs <- ng_motifs
  } else {
    motifs <- ng_s_motifs
    glycan <- glyrepr::convert_glycan_mono_type(glycan, to = "simple", strict = FALSE)
  }
  has_motif_(glycan, motifs$bisecting_core, alignment = "core", ignore_linkages = !strict)
}
