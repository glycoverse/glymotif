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
#' It only checks the monosaccharide types on the "generic" level (e.g. "Hex", "HexNAc"),
#' and ignores linkage information.
#' This is preferred because in most cases the structural resolution could not be high,
#' e.g. in most glycoproteomics studies.
#' However, the glycans are guaranteed to be N-glycans by the glycosylation sites.
#' In this case, we could make some assumptions about the glycan structures,
#' and extract the key properties.
#' For example, an `H-N` terminal motif is considered a terminal galactose.
#' If you have high-resolution glycan structures, you can set `strict = TRUE`.
#'
#' @param glycans A `glyrepr_structure` object,
#' or a character vector of IUPAC-condensed structure strings.
#' @param strict A logical value. If `TRUE`, the glycan must have "concrete"
#' monosaccharides (e.g. "GlcNAc", "Man", "Gal") and linkage information.
#' If `FALSE`, the function is more lenient,
#' checking monosaacharide identities on the "generic" level (e.g. "Hex", "HexNAc")
#' and ignoring linkage information.
#' Default is `FALSE`. This is preferred because in most cases the
#' structural resolution could not be high, but we known for sure the glycans are indeed N-glycans.
#'
#' @return A tibble with the following columns: "glycan_type", "bisecting",
#' "antennae", "core_fuc", "arm_fuc", "terminal_gal".
#' If the input glycans have names, the tibble will have a "glycan" column.
#' Otherwise, IUPAC-condensed strings will be used as the "glycan" column.
#'
#' @examples
#' library(purrr)
#'
#' glycans <- c(
#'   "Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",
#'   "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
#' )
#' describe_n_glycans(glycans)
#'
#' @export
describe_n_glycans <- function(glycans, strict = FALSE) {
  checkmate::assert_flag(strict)

  glycan_structures <- ensure_glycans_are_structures(glycans)
  glycan_names <- prepare_struc_names(glycans, glycan_structures)

  # Deal with strictness
  if (!strict) {
    glycan_structures <- glyrepr::convert_mono_type(glycan_structures, to = "generic")
  }

  # Separate motifs that only need presence/absence from those that need counting
  have_motif_names <- c("core", "pauciman", "hybrid", "highman", "bisect", "ant2", "ant3", "ant4")
  have_motifs <- purrr::map(have_motif_names, ~ get_n_glycan_motif(.x, generic = !strict))

  have_alignments <- c(
    core = "core", pauciman = "whole", hybrid = "core", highman = "core",
    bisect = "core", ant2 = "core", ant3 = "core", ant4 = "core"
  )

  count_motif_names <- c("core_fuc", "arm_fuc", "gal", "terminal_gal")
  count_motifs <- purrr::map(count_motif_names, ~ get_n_glycan_motif(.x, generic = !strict))

  count_alignments <- c(core_fuc = "core", arm_fuc = "core", gal = "substructure", terminal_gal = "terminal")

  # Use vectorized functions
  motif_has_matrix <- have_motifs_(
    glycans = glycan_structures, 
    motifs = have_motifs, 
    alignments = have_alignments, 
    glycan_names = glycan_names, 
    motif_names = have_motif_names, 
    ignore_linkages = !strict
  )
  motif_count_matrix <- count_motifs_(
    glycans = glycan_structures, 
    motifs = count_motifs, 
    alignments = count_alignments, 
    glycan_names = glycan_names, 
    motif_names = count_motif_names, 
    ignore_linkages = !strict
  )
  rownames(motif_has_matrix) <- glycan_names
  rownames(motif_count_matrix) <- glycan_names
  colnames(motif_has_matrix) <- have_motif_names
  colnames(motif_count_matrix) <- count_motif_names

  # Check if the glycans are N-glycans
  is_n_glycan <- motif_has_matrix[, "core"]
  invalid_indices <- which(!is_n_glycan)
  if (length(invalid_indices) > 0) {
    cli::cli_abort("Glycans at indices {.val {invalid_indices}} are not N-glycans.")
  }

  # Determine glycan types based on motif presence
  # Extract logical vectors for each motif check
  pauciman_results <- motif_has_matrix[, "pauciman"]
  hybrid_results <- motif_has_matrix[, "hybrid"]
  highman_results <- motif_has_matrix[, "highman"]
  
  # Use vectorized operations instead of purrr for better performance
  glycan_type <- ifelse(pauciman_results, "paucimannose", 
                 ifelse(hybrid_results, "hybrid",
                 ifelse(highman_results, "highmannose", "complex")))
  
  # Determine number of antennae for complex glycans
  # Extract logical vectors for antenna checks
  ant4_results <- motif_has_matrix[, "ant4"]
  ant3_results <- motif_has_matrix[, "ant3"]
  ant2_results <- motif_has_matrix[, "ant2"]
  
  # Use vectorized operations for antenna counting
  n_antennae <- ifelse(glycan_type != "complex", NA_integer_,
                ifelse(ant4_results, 4L,
                ifelse(ant3_results, 3L,
                ifelse(ant2_results, 2L, 1L))))

  # Build the result tibble
  res <- tibble::tibble(
    glycan_type = glycan_type,
    bisecting = motif_has_matrix[, "bisect"],
    n_antennae = n_antennae,
    n_core_fuc = motif_count_matrix[, "core_fuc"],
    n_arm_fuc = motif_count_matrix[, "arm_fuc"],
    n_gal = motif_count_matrix[, "gal"],
    n_terminal_gal = as.integer(motif_count_matrix[, "terminal_gal"])
  )

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
#' checking monosaacharide identities on the "generic" level (e.g. "Hex", "HexNAc")
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
#' checking monosaacharide identities on the "generic" level (e.g. "Hex", "HexNAc")
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
# - `.count_motif`: A function that counts the number of motifs in a glycan.
#   It should have the signature `.count_motif(glycan, motif, alignment)`.
#
# In these functions, use `.has_motif` and `.count_motif` to check for motifs.
.is_n_glycan <- function(glycan, .has_motif, .count_motif) {
  .has_motif(glycan, "core", alignment = "core")
}


.n_glycan_type <- function(glycan, .has_motif, .count_motif) {
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


.has_bisecting <- function(glycan, .has_motif, .count_motif) {
  .has_motif(glycan, "bisect", alignment = "core")
}


.n_antennae <- function(glycan, .has_motif, .count_motif, is_complex = NULL) {
  if (is.null(is_complex)) {
    is_complex <- .n_glycan_type(glycan, .has_motif, .count_motif) == "complex"
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


.n_core_fuc <- function(glycan, .has_motif, .count_motif) {
  .count_motif(glycan, "core_fuc", alignment = "core")
}


.n_arm_fuc <- function(glycan, .has_motif, .count_motif) {
  .count_motif(glycan, "arm_fuc", alignment = "core")
}


.n_gal <- function(glycan, .has_motif, .count_motif) {
  .count_motif(glycan, "gal", alignment = "substructure")
}


.n_terminal_gal <- function(glycan, .has_motif, .count_motif) {
  .count_motif(glycan, "gal", alignment = "terminal")
}


# ----- Utilities -----
n_glycan_property_wrapper <- function(glycan, strict, func, check_n_glycan = TRUE) {
  # This function encapsulates the common pattern of checking N-glycan properties.
  # To use it, write a function that takes a glycan graph,
  # a `.has_motif` argument taking a function that checks for a motif,
  # and a `.count_motif` argument taking a function that counts motifs.
  # Inside the function, use `.has_motif` and `.count_motif` to check for motifs,
  # instead of calling `has_motif_()` and `count_motif_()` directly.
  #
  # For example:
  # ```
  # .has_bisecting <- function(glycan, .has_motif, .count_motif) {
  #   bisect_graph <- get_motif_graph("N-glycan core, bisected")
  #   .have_motif(glycan, bisect_graph, alignment = "core")
  # }
  # ```
  # Then pass this function to `n_glycan_property_wrapper()`:
  # ```
  # n_glycan_property_wrapper(glycan, strict, .has_bisecting)
  # ```
  # This function will take care of converting the glycan to a graph and
  # checking the motif with suitable strictness.
  checkmate::assert_flag(strict)
  glycan <- ensure_glycans_are_structures(glycan)
  if (!strict) {
    glycan <- glyrepr::convert_mono_type(glycan, to = "generic")
  }
  have_motif_func <- purrr::partial(have_n_glycan_motif, strict = strict)
  count_motif_func <- purrr::partial(count_n_glycan_motif, strict = strict)
  if (check_n_glycan && !.is_n_glycan(glycan, have_motif_func, count_motif_func)) {
    rlang::abort("Not an N-glycan.")
  }
  func(glycan, have_motif_func, count_motif_func)
}


have_n_glycan_motif <- function(glycan, motif_name, alignment, strict) {
  motif <- get_n_glycan_motif(motif_name, generic = !strict)
  # If glycan is an igraph object (from smap), convert it to glyrepr_structure
  if (inherits(glycan, "igraph")) {
    glycan <- glyrepr::as_glycan_structure(glycan)
  }
  have_motif_(glycan, motif, alignment, ignore_linkages = !strict)
}


count_n_glycan_motif <- function(glycan, motif_name, alignment, strict) {
  motif <- get_n_glycan_motif(motif_name, generic = !strict)
  # If glycan is an igraph object (from smap), convert it to glyrepr_structure
  if (inherits(glycan, "igraph")) {
    glycan <- glyrepr::as_glycan_structure(glycan)
  }
  count_motif_(glycan, motif, alignment, ignore_linkages = !strict)
}


get_n_glycan_motif <- function(name, generic = FALSE) {
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
  if (generic) {
    # Here we use H-H-N-H, not just H.
    # This is for telling between Man and Gal.
    motifs[["gal"]] <- glyparse::parse_pglyco_struc("(H(H(N(H))))")
    motifs[["terminal_gal"]] <- motifs[["gal"]]
  } else {
    motifs[["gal"]] <- glyparse::parse_iupac_condensed("Gal")
    motifs[["terminal_gal"]] <- motifs[["gal"]]
  }
  motif <- motifs[[name]]
  if (generic) {
    motif <- glyrepr::convert_mono_type(motif, to = "generic")
  }
  motif
}
get_n_glycan_motif <- memoise::memoise(get_n_glycan_motif)
