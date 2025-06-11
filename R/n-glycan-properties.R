# ----- Interface -----

#' Describe N-Glycans Properties
#'
#' Extract key properties of N-glycans, including:
#' - "glycan_type": N-glycan type: high mannose, hybrid, complex, or paucimannose.
#' - "bisecting": Bisecting GlcNAc presence.
#' - "antennae": Number of antennae.
#' - "core_fuc": Number of core fucoses.
#' - "arm_fuc": Number of arm fucoses.
#' - "terminal_gal": Number of terminal galactoses.
#'
#' This function is designed to work with N-glycans only.
#' If the glycans are not N-glycans, an error is thrown.
#'
#' @details
#' # Strictness
#'
#' By default (`strict = FALSE`), the function is very lenient for motif checking.
#' It only checks the monosaccharide types on the "simple" level (e.g. "H", "N", "F"),
#' and ignores linkage information.
#' This is preferred because in most cases the structural resolution could not be high,
#' e.g. in most glycoproteomics studies.
#' However, the glycans are guaranteed to be N-glycans by the glycosylation sites.
#' In this case, we could make some assumptions about the glycan structures,
#' and extract the key properties.
#' For example, an `H-N` terminal motif is considered a terminal galactose.
#' If you have high-resolution glycan structures, you can set `strict = TRUE`.
#'
#' # Enabling parallel processing
#'
#' This function can spend a lot of time on large datasets (e.g. > 500 glycans).
#' To speed up, you can enable parallel processing by setting `parallel = TRUE`.
#' However, changing the argument to `TRUE` only set the function "ready"
#' for parallel processing.
#' You still need to call [future::plan()] to change the parallel backend.
#' For example, to use the "multisession" backend:
#'
#' ```r
#' library(future)
#' old_plan <- future::plan("multisession")  # Save the old plan
#' describe_n_glycans(glycans, parallel = TRUE)
#' future::plan(old_plan)  # Restore the old plan
#' ```
#'
#' @param glycans A list of `glycan_graph` objects,
#' or a character vector of IUPAC-condensed structure strings.
#' @param strict A logical value. If `TRUE`, the glycan must have "concrete"
#' monosaccharides (e.g. "GlcNAc", "Man", "Gal") and linkage information.
#' If `FALSE`, the function is more lenient,
#' checking monosaacharide identities on the "simple" level (e.g. "H", "N", "F")
#' and ignoring linkage information.
#' Default is `FALSE`. This is preferred because in most cases the
#' structural resolution could not be high, but we known for sure the glycans are indeed N-glycans.
#' @param parallel A logical value. If `TRUE`, the function will use parallel processing.
#' Remember to call [future::plan()] before using this argument,
#' otherwise the function will still use sequential processing.
#'
#' @return A tibble with the following columns: "glycan_type", "bisecting",
#' "antennae", "core_fuc", "arm_fuc", "terminal_gal".
#' If the input glycans have names, the tibble will have a "glycan" column.
#' Otherwise, if IUPAC condensed strings are used, they will be used as the "glycan" column.
#'
#' @examples
#' library(glyparse)
#' library(purrr)
#'
#' glycans <- c(
#'   "(N(F)(N(H(H(N))(H(N(H))))))",
#'   "(N(N(H(H)(H(H)(H)))))",
#'   "(N(F)(N(H(H(N))(H(N(H(H)))))))",
#'   "(N(N(H(N)(H(N(H)(F)))(H(N(H)(F))(N(H)(F))))))",
#'   "(N(N(H(H(N(H(A))))(H(N(H(A)))))))"
#' )
#' glycans <- map(glycans, parse_pglyco_struc)
#' describe_n_glycans(glycans)
#'
#' @export
describe_n_glycans <- function(glycans, strict = FALSE, parallel = FALSE) {
  # Validate the input
  valid_glycans_arg(glycans)
  checkmate::assert_flag(strict)
  checkmate::assert(checkmate::check_null(parallel), checkmate::check_flag(parallel))
  glycans <- ensure_glycans_are_graphs(glycans)

  # Deal with strictness
  if (!strict) {
    glycans <- purrr::map(
      glycans, glyrepr::convert_glycan_mono_type,
      to = "simple", strict = FALSE
    )
  }
  hf <- purrr::partial(has_n_glycan_motif, strict = strict)
  cf <- purrr::partial(count_n_glycan_motif, strict = strict)

  # Check if the glycans are N-glycans
  invalid_indices <- which(!purrr::map_lgl(glycans, .is_n_glycan, hf, cf))
  if (length(invalid_indices) > 0) {
    cli::cli_abort("Glycans at indices {.val {invalid_indices}} are not N-glycans.")
  }

  # Get the properties
  glycan_type <- purrr::map_chr(glycans, .n_glycan_type, hf, cf)
  res <- tibble::tibble(
    glycan_type = glycan_type,
    bisecting = purrr::map_lgl(glycans, .has_bisecting, hf, cf),
    n_antennae = purrr::map2_int(
      glycans, glycan_type == "complex",
      ~ .n_antennae(.x, hf, cf, is_complex = .y)
    ),
    n_core_fuc = purrr::map_int(glycans, .n_core_fuc, hf, cf),
    n_arm_fuc = purrr::map_int(glycans, .n_arm_fuc, hf, cf),
    n_gal = purrr::map_int(glycans, .n_gal, hf, cf),
    n_terminal_gal = purrr::map_int(glycans, .n_terminal_gal, hf, cf)
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
#' @param glycan A `glycan_graph` object, or a character string of IUPAC condensed format.
#' @param strict A logical value. If `TRUE`, the glycan must have "concrete"
#' monosaccharides (e.g. "GlcNAc", "Man", "Gal") and linkage information.
#' If `FALSE`, the function is more lenient,
#' checking monosaacharide identities on the "simple" level (e.g. "H", "N", "F")
#' and ignoring linkage information.
#' Default is `FALSE`. This is preferred because in most cases the
#' structural resolution could not be high, but we known for sure the glycans are indeed N-glycans.
#'
#' @return A logical value.
#' @export
is_n_glycan <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .is_n_glycan, check_n_glycan = FALSE)
}


#' Determine N-Glycan Key Properties
#'
#' These functions check key properties of an N-glycan:
#' - `n_glycan_type()`: Determine the N-glycan type.
#' - `has_bisecting()`: Check if the glycan has a bisecting GlcNAc.
#' - `n_antennae()`: Count the number of antennae.
#' - `n_core_fuc()`: Count the number of core fucoses.
#' - `n_arm_fuc()`: Count the number of arm fucoses.
#' - `n_gal()`: Count the number of galactoses.
#' - `n_terminal_gal()`: Count the number of terminal galactoses.
#'
#' @details
#' # `n_glycan_type()`: N-Glycan Types
#'
#' Four types of N-glycans are recognized: high mannose, hybrid, complex, and paucimannose.
#' For more information about N-glycan types,
#' see [Essentials of Glycobiology](https://www.ncbi.nlm.nih.gov/books/NBK579964/#_s9_2_).
#'
#' # `has_bisecting()`: Bisecting GlcNAc
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
#' # `n_antennae()`: Number of Antennae
#'
#' The number of antennae is the number of branching GlcNAc to the core mannoses
#' in a complex N-glycan. Bisecting GlcNAc is not counted as an antenna.
#' This functions returns NA_integer_ for non-complex N-glycans.
#'
#' # `n_core_fuc()`: Number of Core Fucoses
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
#' # `n_arm_fuc()`: Number of Arm Fucoses
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
#' # `n_gal()`: Number of Galactoses
#'
#' This function seems useless and silly.
#' It is, if you have a well-structured glycan with concrete monosaccharides.
#' However, if you only have "Hex" or "H" at hand,
#' it is tricky to know how many of them are "Gal" and how many are "Man".
#' This function makes a simply assumption that all the rightmost "H" in a
#' "H-H-N-H" unit is a galactose.
#' The two "H" on the left are mannoses of the N-glycan core.
#' The "N" is a GlcNAc attached to one core mannose.
#'
#' This method is only used when `strict = FALSE`.
#' If `strict = TRUE`, the function will authentically check the number of "Gal".
#'
#' # `n_terminal_gal()`: Number of Terminal Galactoses
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
#' @param glycan A `glycan_graph` object, or a character string of IUPAC condensed format.
#' @param strict A logical value. If `TRUE`, the glycan must have "concrete"
#' monosaccharides (e.g. "GlcNAc", "Man", "Gal") and linkage information.
#' If `FALSE`, the function is more lenient,
#' checking monosaacharide identities on the "simple" level (e.g. "H", "N", "F")
#' and ignoring linkage information.
#' Default is `FALSE`. This is preferred because in most cases the
#' structural resolution could not be high, but we known for sure the glycans are indeed N-glycans.
#'
#' @returns
#' - `n_glycan_type()`: A character scalar indicating the N-glycan type,
#'   either "highmannose", "hybrid", "complex", or "paucimannose".
#' - `has_bisecting()`: A logical value indicating if the glycan has a bisecting GlcNAc.
#' - `n_antennae()`: An integer scalar indicating the number of antennae.
#' - `n_core_fuc()`: An integer scalar indicating the number of core fucoses.
#' - `n_arm_fuc()`: An integer scalar indicating the number of arm fucoses.
#' - `n_gal()`: An integer scalar indicating the number of galactoses.
#' - `n_terminal_gal()`: An integer scalar indicating the number of terminal galactoses.
#'
#' @export
n_glycan_type <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_glycan_type)
}


#' @rdname n_glycan_type
#' @export
has_bisecting <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .has_bisecting)
}


#' @rdname n_glycan_type
#' @export
n_antennae <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_antennae)
}


#' @rdname n_glycan_type
#' @export
n_core_fuc <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_core_fuc)
}


#' @rdname n_glycan_type
#' @export
n_arm_fuc <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_arm_fuc)
}


#' @rdname n_glycan_type
#' @export
n_gal <- function(glycan, strict = FALSE) {
  n_glycan_property_wrapper(glycan, strict, .n_gal)
}


#' @rdname n_glycan_type
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
  .has_motif(glycan, "core", alignment = "core")
}


.n_glycan_type <- function(glycan, .has_motif, .counts_motif) {
  if (.has_motif(glycan, "core", alignment = "whole") ||
      .has_motif(glycan, "pauciman", alignment = "whole")) {
    "paucimannose"
  } else if (.has_motif(glycan, "hybrid", alignment = "core")) {
    "hybrid"
  } else if (.has_motif(glycan, "highman", alignment = "core")) {
    "highmannose"
  } else {
    "complex"
  }
}


.has_bisecting <- function(glycan, .has_motif, .counts_motif) {
  .has_motif(glycan, "bisect", alignment = "core")
}


.n_antennae <- function(glycan, .has_motif, .counts_motif, is_complex = NULL) {
  if (is.null(is_complex)) {
    is_complex <- .n_glycan_type(glycan, .has_motif, .counts_motif) == "complex"
  }
  if (!is_complex) {
    return(NA_integer_)
  }
  if (.has_motif(glycan, "ant4", alignment = "core")) {
    4L
  } else if (.has_motif(glycan, "ant3", alignment = "core")) {
    3L
  } else if (.has_motif(glycan, "ant2", alignment = "core")) {
    2L
  } else {
    1L
  }
}


.n_core_fuc <- function(glycan, .has_motif, .counts_motif) {
  .counts_motif(glycan, "core_fuc", alignment = "core")
}


.n_arm_fuc <- function(glycan, .has_motif, .counts_motif) {
  .counts_motif(glycan, "arm_fuc", alignment = "core")
}


.n_gal <- function(glycan, .has_motif, .counts_motif) {
  .counts_motif(glycan, "gal", alignment = "substructure")
}


.n_terminal_gal <- function(glycan, .has_motif, .counts_motif) {
  .counts_motif(glycan, "gal", alignment = "terminal")
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
  if (!strict) {
    glycan <- glyrepr::convert_glycan_mono_type(glycan, to = "simple", strict = FALSE)
  }
  has_motif_func <- purrr::partial(has_n_glycan_motif, strict = strict)
  counts_motif_func <- purrr::partial(count_n_glycan_motif, strict = strict)
  if (check_n_glycan && !.is_n_glycan(glycan, has_motif_func, counts_motif_func)) {
    rlang::abort("Not an N-glycan.")
  }
  func(glycan, has_motif_func, counts_motif_func)
}


has_n_glycan_motif <- function(glycan, motif_name, alignment, strict) {
  motif <- get_n_glycan_motif(motif_name, simple = !strict)
  has_motif_(glycan, motif, alignment, ignore_linkages = !strict)
}


count_n_glycan_motif <- function(glycan, motif_name, alignment, strict) {
  motif <- get_n_glycan_motif(motif_name, simple = !strict)
  counts_motif_(glycan, motif, alignment, ignore_linkages = !strict)
}


get_n_glycan_motif <- function(name, simple = FALSE) {
  motifs <- list(
    core     = get_motif_structure("N-Glycan core basic"),
    pauciman = glyparse::parse_iupac_condensed("Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"),
    hybrid   = get_motif_structure("N-Glycan hybrid"),
    highman  = get_motif_structure("N-Glycan high mannose"),
    bisect   = get_motif_structure("N-glycan core, bisected"),
    ant2     = get_motif_structure("N-Glycan biantennary"),
    ant3     = get_motif_structure("N-Glycan triantennary"),
    ant4     = get_motif_structure("N-Glycan tetraantennary"),
    core_fuc = get_motif_structure("N-Glycan core, core-fucosylated"),
    arm_fuc  = get_motif_structure("N-Glycan core, arm-fucosylated")
  )
  # add gal
  if (simple) {
    # Here we use H-H-N-H, not just H.
    # This is for telling between Man and Gal.
    motifs[["gal"]] <- glyparse::parse_pglyco_struc("(H(H(N(H))))")
  } else {
    motifs[["gal"]] <- glyparse::parse_iupac_condensed("Gal")
  }
  motif <- motifs[[name]]
  if (simple) {
    motif <- glyrepr::convert_glycan_mono_type(motif, to = "simple", strict = FALSE)
  }
  motif
}
get_n_glycan_motif <- memoise(get_n_glycan_motif)
