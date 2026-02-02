#' Add Motif Annotations
#'
#' This function adds motif annotations to the variable information
#' of a [glyexp::experiment()] or a tibble with a structure column.
#' `add_motifs_int()` adds integer annotations (how many motifs are present).
#' `add_motifs_lgl()` adds boolean annotations (whether the motif is present).
#'
#' @section About Names:
#'
#' The naming rule for the new columns follows these priorities:
#' 1. If `motifs` is a named vector (character or `glyrepr::glycan_structure()`),
#'    the names are used directly as column names.
#' 2. If `motifs` is unnamed and contains known motif names (e.g., "N-Glycan core"),
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
#' @section Why do we need these functions:
#'
#' Adding one motif annotation to a [glyexp::experiment()] is easy:
#'
#' ```r
#' exp |>
#'   mutate_var(has_hex = have_motif(glycan_structure, "Hex"))
#' ```
#'
#' However, adding multiple motifs is not as straightforward.
#' You can still use `mutate_var()` to add multiple motifs like this:
#'
#' ```r
#' exp |>
#'   mutate_var(
#'     n_hex = count_motif(glycan_structure, "Hex"),
#'     n_dhex = count_motif(glycan_structure, "dHex"),
#'     n_hexnac = count_motif(glycan_structure, "HexNAc"),
#'   )
#' ```
#'
#' This method has two problems:
#'
#' 1. it has a lot of boilerplate code (a lot of typing)
#' 2. it is not very efficient, as each call to `count_motif`
#'    performs validation and conversion on `glycan_structure`,
#'    which is a time-consuming process.
#'
#' Therefore, we think it would be better to have a function that
#' adds multiple motif annotations in a single call, in a more intuitive way.
#' That's why we provide these two functions.
#'
#' Under the hood, they use a more straightforward approach for [glyexp::experiment()] objects:
#'
#' 1. get the motif annotation matrix using `count_motifs()` or `have_motifs()`
#' 2. convert the matrix to a tibble
#' 3. use `dplyr::bind_cols()` to add the tibble to the variable information
#'
#' @param x A [glyexp::experiment()] object, or a tibble with a structure column.
#' @param ... Additional arguments passed to the method.
#' @inheritParams have_motifs
#'
#' @return An [glyexp::experiment()] object with motif annotations added to the variable information.
#'
#' @examples
#' library(glyexp)
#'
#' exp <- real_experiment2
#'
#' exp |>
#'   add_motifs_lgl(c(
#'     lacnac = "Gal(??-?)GlcNAc(??-",
#'     sia_lacnac = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-"
#'   )) |>
#'   get_var_info()
#'
#' exp |>
#'   add_motifs_int(c(
#'     lacnac = "Gal(??-?)GlcNAc(??-",
#'     sia_lacnac = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-"
#'   )) |>
#'   get_var_info()
#'
#' @seealso [glymotif::have_motifs()], [glymotif::count_motifs()], [glyexp::experiment()]
#'
#' @export
add_motifs_int <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  UseMethod("add_motifs_int")
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  UseMethod("add_motifs_lgl")
}

#' @rdname add_motifs_int
#' @export
add_motifs_int.glyexp_experiment <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  .add_motifs_anno_exp(x, count_motifs, motifs, alignments, ignore_linkages, strict_sub, match_degree)
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.glyexp_experiment <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  .add_motifs_anno_exp(x, have_motifs, motifs, alignments, ignore_linkages, strict_sub, match_degree)
}

.add_motifs_anno_exp <- function(exp, motif_anno_fn, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL) {
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort("The experiment must have a {.field glycan_structure} column.")
  }

  motif_anno <- motif_anno_fn(exp$var_info$glycan_structure, motifs, alignments, ignore_linkages, strict_sub, match_degree)
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
add_motifs_int.data.frame <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  .add_motifs_anno_df(x, count_motifs, motifs, alignments, ignore_linkages, strict_sub, match_degree)
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.data.frame <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL, ...) {
  .add_motifs_anno_df(x, have_motifs, motifs, alignments, ignore_linkages, strict_sub, match_degree)
}

# Helper function to get column names for motif annotations
# Follows the rules:
# 1. If motifs is a named vector, use the names
# 2. If motifs is unnamed and a vector of known motif names, use the motif names
# 3. If motifs is unnamed and a glyrepr_structure or character vector of structure strings, use IUPAC strings
.get_motif_colnames <- function(motifs) {
  motif_names <- prepare_motif_names(motifs)
  if (!is.null(motif_names)) {
    return(motif_names)
  }

  # Fallback to IUPAC strings for unnamed structures or IUPAC strings
  as.character(motifs)
}

.add_motifs_anno_df <- function(df, motif_anno_fn, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, match_degree = NULL) {
  if (!"glycan_structure" %in% colnames(df)) {
    cli::cli_abort("A {.field glycan_structure} column is required.")
  }
  motif_anno <- motif_anno_fn(df$glycan_structure, motifs, alignments, ignore_linkages, strict_sub, match_degree)
  # have_motifs/count_motifs set colnames for motif specs, but may return NULL for regular motifs
  # We always need colnames for add_motifs, so generate them if missing
  if (is.null(colnames(motif_anno))) {
    colnames(motif_anno) <- .get_motif_colnames(motifs)
  }
  motif_anno <- tibble::as_tibble(motif_anno)
  dplyr::bind_cols(df, motif_anno)
}
