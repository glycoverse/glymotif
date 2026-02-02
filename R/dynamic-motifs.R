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
