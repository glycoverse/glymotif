#' Add Motif Annotations to an Experiment
#'
#' This function adds motif annotations to the variable information
#' of a [glyexp::experiment()].
#' `add_motifs_int()` adds integer annotations (how many motifs are present).
#' `add_motifs_lgl()` adds boolean annotations (whether the motif is present).
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
#' Advanced R users might want to use `count_motifs()`
#' (the plural cousin of `count_motif()`) with `!!!`:
#'
#' ```r
#' exp |>
#'   mutate_var(!!!count_motifs(glycan_structure, c("Hex", "dHex", "HexNAc")))
#' ```
#'
#' Sadly, this doesn't work.
#' Firstly, `count_motifs` returns a matrix, not a list.
#' Secondly, even if you use `as.data.frame()` to convert it to a list,
#' `!!!` triggers early evaluation of `glycan_structure` in the calling environment,
#' before passing it to `count_motifs()`.
#' This will raise an "object not found" error,
#' and there is no easy way to fix this, at least for now.
#'
#' Therefore, we think it would be better to have a function that
#' adds multiple motif annotations in a single call, in a more intuitive way.
#' That's why we provide these two functions.
#'
#' Under the hood, they use a more straightforward approach:
#'
#' 1. get the motif annotation matrix using `count_motifs()` or `have_motifs()`
#' 2. convert the matrix to a tibble
#' 3. use `dplyr::bind_cols()` to add the tibble to the variable information
#'
#' @param exp An [glyexp::experiment()] object.
#' @inheritParams have_motif
#'
#' @return An [glyexp::experiment()] object with motif annotations added to the variable information.
#' @seealso [glymotif::have_motifs()], [glymotif::count_motifs()], [glyexp::experiment()]
#'
#' @export
add_motifs_int <- function(exp, motifs, alignments = NULL, ignore_linkages = FALSE) {
  .add_motifs_anno(exp, count_motifs, motifs, alignments, ignore_linkages)
}

#' @rdname add_motifs_int
#' @export
add_motifs_lgl <- function(exp, motifs, alignments = NULL, ignore_linkages = FALSE) {
  .add_motifs_anno(exp, have_motifs, motifs, alignments, ignore_linkages)
}

.add_motifs_anno <- function(exp, motif_anno_fn, motifs, alignments = NULL, ignore_linkages = FALSE) {
  checkmate::assert_class(exp, "glyexp_experiment")
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort("The experiment must have a {.field glycan_structure} column.")
  }
  checkmate::assert_class(exp$var_info$glycan_structure, "glyrepr_structure")

  motif_anno <- motif_anno_fn(exp$var_info$glycan_structure, motifs, alignments, ignore_linkages)
  names(motif_anno) <- as.character(motifs)
  motif_anno <- tibble::as_tibble(motif_anno)
  exp$var_info <- dplyr::bind_cols(exp$var_info, motif_anno)
  exp
}
