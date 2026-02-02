#' Dynamic Motifs Specification
#'
#' Create a specification for dynamic motif extraction.
#' This should be passed to the `motifs` argument of [have_motifs()],
#' [count_motifs()], [match_motifs()], [add_motifs_lgl()], or [add_motifs_int()].
#'
#' When used, motifs will be extracted dynamically from the input glycans
#' using [extract_motif()], and motif matching will use `alignment = "substructure"`.
#'
#' @param max_size The maximum number of monosaccharides in the extracted motifs.
#'   Default is 3. Passed to [extract_motif()].
#'
#' @returns A `dynamic_motifs_spec` object.
#'
#' @seealso [branch_motifs()], [have_motifs()], [count_motifs()]
#'
#' @examples
#' # Use in have_motifs()
#' # have_motifs(glycans, dynamic_motifs(max_size = 4))
#'
#' @export
dynamic_motifs <- function(max_size = 3) {
  checkmate::assert_int(max_size, lower = 1)
  structure(list(max_size = max_size), class = "dynamic_motifs_spec")
}

#' @export
print.dynamic_motifs_spec <- function(x, ...) {
  cli::cli_text("<{.cls dynamic_motifs_spec}>")
  cli::cli_text("This object should be passed to the {.arg motifs} argument of {.fn have_motifs}, {.fn count_motifs}, {.fn match_motifs}, {.fn add_motifs_lgl}, or {.fn add_motifs_int}.")
  cli::cli_text("Configuration: {.field max_size} = {.val {x$max_size}}")
  invisible(x)
}

#' Branch Motifs Specification
#'
#' Create a specification for branch motif extraction.
#' This should be passed to the `motifs` argument of [have_motifs()],
#' [count_motifs()], [match_motifs()], [add_motifs_lgl()], or [add_motifs_int()].
#'
#' When used, motifs will be extracted dynamically from the input glycans
#' using [extract_branch_motif()] with `including_core = TRUE`.
#' Motif matching will use a custom `match_degree` where the last 4 nodes
#' of each motif are not required to match degree exactly (to allow for
#' attachment to the glycan core).
#'
#' @returns A `branch_motifs_spec` object.
#'
#' @seealso [dynamic_motifs()], [have_motifs()], [count_motifs()]
#'
#' @examples
#' # Use in have_motifs()
#' # have_motifs(glycans, branch_motifs())
#'
#' @export
branch_motifs <- function() {
  structure(list(), class = "branch_motifs_spec")
}

#' @export
print.branch_motifs_spec <- function(x, ...) {
  cli::cli_text("<{.cls branch_motifs_spec}>")
  cli::cli_text("This object should be passed to the {.arg motifs} argument of {.fn have_motifs}, {.fn count_motifs}, {.fn match_motifs}, {.fn add_motifs_lgl}, or {.fn add_motifs_int}.")
  cli::cli_text("Configuration: extracts branch motifs with core included")
  invisible(x)
}

#' Resolve Motif Specification
#'
#' Internal function to resolve a motif specification into actual motifs
#' and matching parameters.
#'
#' @param glycans A `glyrepr_structure` object.
#' @param spec A `dynamic_motifs_spec` or `branch_motifs_spec` object.
#' @param alignments User-provided alignments (should be NULL).
#' @param match_degree User-provided match_degree (should be NULL).
#'
#' @returns A list with `motifs`, `alignments`, and `match_degree`.
#' @noRd
resolve_motif_spec <- function(glycans, spec, alignments, match_degree) {
  if (!is.null(alignments)) {
    cli::cli_abort(c(
      "Cannot specify {.arg alignments} when using {.fn dynamic_motifs} or {.fn branch_motifs}.",
      "i" = "Alignment is controlled automatically by the algorithm."
    ))
  }
  if (!is.null(match_degree)) {
    cli::cli_abort(c(
      "Cannot specify {.arg match_degree} when using {.fn dynamic_motifs} or {.fn branch_motifs}.",
      "i" = "Match degree is controlled automatically by the algorithm."
    ))
  }

  UseMethod("resolve_motif_spec", spec)
}

#' @export
resolve_motif_spec.dynamic_motifs_spec <- function(glycans, spec, alignments, match_degree) {
  motifs <- extract_motif(glycans, max_size = spec$max_size)

  list(
    motifs = motifs,
    alignments = rep("substructure", length(motifs)),
    match_degree = NULL
  )
}

#' Trim Branch Motif IUPAC String
#
# Removes the last 4 monosaccharides from a branch motif IUPAC string.
# The core (4 residues: Man-Man-GlcNAc-GlcNAc or Hex-Hex-HexNAc-HexNAc)
# is appended by extract_branch_motif(), but for display purposes we only want the branch part.
#
# @param iupac A character vector of IUPAC-condensed strings.
# @return Trimmed IUPAC strings.
# @noRd
.trim_branch_iupac <- function(iupac) {
  # Remove the core suffix pattern for concrete mono types using fixed string matching
  # Pattern: ?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-
  concrete_suffix <- "?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  result <- stringr::str_remove(iupac, stringr::fixed(concrete_suffix))

  # Remove the core suffix pattern for generic mono types
  generic_suffix <- "?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"
  result <- stringr::str_remove(result, stringr::fixed(generic_suffix))

  result
}

resolve_motif_spec.branch_motifs_spec <- function(glycans, spec, alignments, match_degree) {
  motifs <- extract_branch_motif(glycans, including_core = TRUE)

  # Construct match_degree: for each motif, last 4 nodes are FALSE, others TRUE
  motif_graphs <- glyrepr::get_structure_graphs(motifs)
  if (inherits(motif_graphs, "igraph")) {
    motif_graphs <- list(motif_graphs)
  }

  match_degree <- purrr::map(motif_graphs, function(g) {
    n_nodes <- igraph::vcount(g)
    if (n_nodes <= 4) {
      rep(FALSE, n_nodes)
    } else {
      c(rep(TRUE, n_nodes - 4), rep(FALSE, 4))
    }
  })

  # Set trimmed IUPAC strings as names for display purposes
  if (length(motifs) > 0) {
    iupac_full <- as.character(motifs)
    iupac_trimmed <- .trim_branch_iupac(iupac_full)
    names(motifs) <- iupac_trimmed
  }

  list(
    motifs = motifs,
    alignments = rep("substructure", length(motifs)),  # Ignored when match_degree is provided
    match_degree = match_degree
  )
}
