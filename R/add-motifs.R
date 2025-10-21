#' Add Motif Annotations
#'
#' This function adds motif annotations to the variable information
#' of a [glyexp::experiment()] or a tibble with a structure column.
#' `add_motifs_int()` adds integer annotations (how many motifs are present).
#' `add_motifs_lgl()` adds boolean annotations (whether the motif is present).
#'
#' @section About Names:
#'
#' The naming rule for the new columns is similar to that of [have_motifs()].
#' Briefly, you can use named character vector to name the motifs,
#' and that will be used as the new column names.
#' The only catchup is that you cannot pass a named `glyrepr::glycan_structure()` to `motifs`.
#' This is a fundamental limitation of the `vctrs_rcrd` class,
#' which `glyrepr::glycan_structure()` is built on.
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
#' @seealso [glymotif::have_motifs()], [glymotif::count_motifs()], [glyexp::experiment()]
#'
#' @export
add_motifs_int <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  UseMethod("add_motifs_int")
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  UseMethod("add_motifs_lgl")
}

#' @rdname add_motifs_int
#' @export
add_motifs_int.glyexp_experiment <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  .add_motifs_anno_exp(x, count_motifs, motifs, alignments, ignore_linkages)
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.glyexp_experiment <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  .add_motifs_anno_exp(x, have_motifs, motifs, alignments, ignore_linkages)
}

.add_motifs_anno_exp <- function(exp, motif_anno_fn, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE) {
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort("The experiment must have a {.field glycan_structure} column.")
  }

  motif_anno <- motif_anno_fn(exp$var_info$glycan_structure, motifs, alignments, ignore_linkages, strict_sub)
  names(motif_anno) <- as.character(motifs)
  motif_anno <- tibble::as_tibble(motif_anno)
  exp$var_info <- dplyr::bind_cols(exp$var_info, motif_anno)
  exp
}

#' @rdname add_motifs_int
#' @export
add_motifs_int.data.frame <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  .add_motifs_anno_df(x, count_motifs, motifs, alignments, ignore_linkages)
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl.data.frame <- function(x, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE, ...) {
  .add_motifs_anno_df(x, have_motifs, motifs, alignments, ignore_linkages)
}

.add_motifs_anno_df <- function(df, motif_anno_fn, motifs, alignments = NULL, ignore_linkages = FALSE, strict_sub = TRUE) {
  if (!"glycan_structure" %in% colnames(df)) {
    cli::cli_abort("A {.field glycan_structure} column is required.")
  }
  motif_anno <- motif_anno_fn(df$glycan_structure, motifs, alignments, ignore_linkages, strict_sub)
  motif_anno <- tibble::as_tibble(motif_anno)
  dplyr::bind_cols(df, motif_anno)
}