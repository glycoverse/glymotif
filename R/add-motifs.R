#' Add Motif Annotations
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `add_motifs_int()` and `add_motifs_lgl()` were deprecated in glymotif
#' 0.17.0. For data frames, use [dplyr::mutate()] with [tibble::as_tibble()]
#' and [count_motifs()] or [have_motifs()]. For [glyexp::GlycomicSE()] or
#' [glyexp::GlycoproteomicSE()] objects, use [glyexp::mutate_row()] with the
#' same tibble expression. The deprecated methods for legacy glyexp containers
#' remain available for compatibility.
#'
#' @section About Names:
#'
#' The naming rule for the new columns follows these priorities:
#' 1. If `motifs` is a named vector (character or `glyrepr::glycan_structure()`),
#'    the names are used directly as column names.
#' 2. If `motifs` is unnamed and contains GGM database motif names (e.g., "N-Glycan core"),
#'    the motif names are used as column names.
#' 3. If `motifs` is unnamed and contains `glyrepr::glycan_structure()` objects
#'    or IUPAC-condensed structure strings, the IUPAC-condensed strings are used
#'    as column names.
#' 4. If `motifs` is a motif spec from [dynamic_motifs()] or [branch_motifs()],
#'    the IUPAC-condensed strings of the extracted motifs are used as column names.
#'
#' Note: This behavior differs from [have_motifs()] and [count_motifs()], which
#' return matrices with NULL column names for unnamed IUPAC string or structure
#' motifs. The functions here always provide column names since they are designed
#' for adding motif annotations to data frames.
#'
#' @param x A legacy glyexp data container or a tibble with a structure column.
#' @param ... Additional arguments passed to the method.
#' @inheritParams have_motifs
#'
#' @return The legacy glyexp data container or data frame supplied in `x`, with
#'   motif annotations added to its variable information.
#'
#' @examples
#' library(glyexp)
#' library(dplyr)
#' library(tibble)
#'
#' glyco_se <- real_experiment2
#' motifs <- c(
#'   lacnac = "Gal(??-?)GlcNAc(??-",
#'   sia_lacnac = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-"
#' )
#'
#' glyco_se |>
#'   mutate_row(as_tibble(count_motifs(glycan_structure, motifs)))
#'
#' df <- tibble(
#'   glycan_structure = c(
#'     "Gal(b1-4)GlcNAc(b1-",
#'     "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
#'   )
#' )
#' df |>
#'   mutate(as_tibble(count_motifs(glycan_structure, motifs)))
#'
#' @seealso [glymotif::have_motifs()], [glymotif::count_motifs()],
#'   [glyexp::GlycomicSE()], [glyexp::GlycoproteomicSE()]
#'
#' @export
add_motifs_int <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  lifecycle::deprecate_warn(
    "0.17.0",
    "add_motifs_int()",
    details = c(
      "For data frames, use `dplyr::mutate(as_tibble(count_motifs(glycan_structure, motifs)))`.",
      "For glyexp containers, use `glyexp::mutate_row(as_tibble(count_motifs(glycan_structure, motifs)))`."
    )
  )
  UseMethod("add_motifs_int")
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  lifecycle::deprecate_warn(
    "0.17.0",
    "add_motifs_lgl()",
    details = c(
      "For data frames, use `dplyr::mutate(as_tibble(have_motifs(glycan_structure, motifs)))`.",
      "For glyexp containers, use `glyexp::mutate_row(as_tibble(have_motifs(glycan_structure, motifs)))`."
    )
  )
  UseMethod("add_motifs_lgl")
}

#' @rdname add_motifs_int
#' @export
add_motifs_int.glyexp_experiment <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  .add_motifs_anno_exp(
    x,
    count_motifs,
    motifs,
    alignments,
    ignore_linkages,
    strict_sub,
    match_degree
  )
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.glyexp_experiment <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  .add_motifs_anno_exp(
    x,
    have_motifs,
    motifs,
    alignments,
    ignore_linkages,
    strict_sub,
    match_degree
  )
}

.add_motifs_anno_exp <- function(
  exp,
  motif_anno_fn,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort(
      "The experiment must have a {.field glycan_structure} column."
    )
  }

  motif_anno <- motif_anno_fn(
    exp$var_info$glycan_structure,
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree
  )
  # have_motifs/count_motifs set colnames for motif specs, but may return NULL for regular motifs
  # We always need colnames for add_motifs, so generate them if missing
  if (is.null(colnames(motif_anno))) {
    colnames(motif_anno) <- .get_motif_colnames(motifs)
  }
  motif_anno <- tibble::as_tibble(motif_anno)
  exp$var_info <- dplyr::bind_cols(exp$var_info, motif_anno)
  exp
}

#' @rdname add_motifs_int
#' @export
add_motifs_int.data.frame <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  .add_motifs_anno_df(
    x,
    count_motifs,
    motifs,
    alignments,
    ignore_linkages,
    strict_sub,
    match_degree
  )
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.data.frame <- function(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
) {
  .add_motifs_anno_df(
    x,
    have_motifs,
    motifs,
    alignments,
    ignore_linkages,
    strict_sub,
    match_degree
  )
}

# Helper function to get column names for motif annotations
# Follows the rules:
# 1. If motifs is a named vector, use the names
# 2. If motifs is unnamed and a vector of GGM database motif names, use the motif names
# 3. If motifs is unnamed and a glyrepr_structure or character vector of structure strings, use IUPAC strings
.get_motif_colnames <- function(motifs) {
  motif_names <- prepare_motif_names(motifs)
  if (!is.null(motif_names)) {
    return(motif_names)
  }

  # Fallback to IUPAC strings for unnamed structures or IUPAC strings
  as.character(motifs)
}

.add_motifs_anno_df <- function(
  df,
  motif_anno_fn,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  if (!"glycan_structure" %in% colnames(df)) {
    cli::cli_abort("A {.field glycan_structure} column is required.")
  }
  motif_anno <- motif_anno_fn(
    df$glycan_structure,
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree
  )
  # have_motifs/count_motifs set colnames for motif specs, but may return NULL for regular motifs
  # We always need colnames for add_motifs, so generate them if missing
  if (is.null(colnames(motif_anno))) {
    colnames(motif_anno) <- .get_motif_colnames(motifs)
  }
  motif_anno <- tibble::as_tibble(motif_anno)
  dplyr::bind_cols(df, motif_anno)
}
