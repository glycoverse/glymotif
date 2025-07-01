#' Add Descriptions to Glycan Compositions
#'
#' This function adds the following columns to the
#' variable information tibble:
#' - 'n_hex': number of Hex
#' - 'n_hexnac': number of HexNAc
#' - 'n_fuc': number of Fuc
#' - 'n_neuac': number of NeuAc
#' - 'n_neugc': number of NeuGc
#' - 'n_sia': number of Sia
#'
#' `n_sia` is the sum of `n_neuac` and `n_neugc`.
#'
#' @details
#' The function relies on the `glycan_composition` column in
#' the `var_info` tibble, and it should be a `glyrepr_composition` object.
#' If you use the `glyread` package to directly read results from search engines
#' (such as pGlyco3), this column will be automatically added to `var_info`.
#' If you construct the `glyexp_experiment` object yourself, you need to ensure
#' that this column exists and its type is `glyrepr_composition`.
#' See [glyrepr::glycan_composition()] or [glyrepr::as_glycan_composition()]
#' for how to create such objects manually.
#'
#' @param exp An [glyexp::experiment()] object.
#'
#' @returns The experiment object with the new columns added.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
add_comp_descriptions <- function(exp) {
  # Check arguments
  checkmate::assert_class(exp, "glyexp_experiment")
  if (!.has_composition_column(exp)) {
    cli::cli_abort("Column {.field glycan_composition} not found in {.field var_info}.")
  }
  checkmate::assert_class(exp$var_info$glycan_composition, "glyrepr_composition")

  # Check if function has been called before
  if (!is.null(exp$meta_data$comp_descriptions_added)) {
    cli::cli_alert_info("Composition descriptions already added. Skipping.")
    return(exp)
  }

  # Add columns
  new_var_info <- dplyr::mutate(
    exp$var_info,
    n_hex = glyrepr::count_mono(.data$glycan_composition, "Hex"),
    n_hexnac = glyrepr::count_mono(.data$glycan_composition, "HexNAc"),
    n_fuc = glyrepr::count_mono(.data$glycan_composition, "dHex"),
    n_neuac = glyrepr::count_mono(.data$glycan_composition, "NeuAc"),
    n_neugc = glyrepr::count_mono(.data$glycan_composition, "NeuGc"),
    n_sia = .data$n_neuac + .data$n_neugc
  )

  # Update experiment and mark as added
  exp$var_info <- new_var_info
  exp$meta_data$comp_descriptions_added <- TRUE
  exp
}


#' Add Descriptions to Glycan Structures
#'
#' This function adds columns about glycan structural properties
#' to the variable information tibble.
#' Depending on the glycan type (N-glycan, O-glycan),
#' different columns are added.
#' Current, only N-glycan descriptions are implemented.
#' See [glymotif::describe_n_glycans()] for columns added.
#'
#' @details
#' The function relies on the `glycan_structure` column in
#' the `var_info` tibble, and it should be a `glyrepr_structure` object.
#' If you use the `glyread` package to directly read results from search engines
#' (such as pGlyco3), this column will be automatically added to `var_info`.
#' If you construct the `glyexp_experiment` object yourself, you need to ensure
#' that this column exists and its type is `glyrepr_structure`.
#' See [glyrepr::glycan_structure()] or [glyrepr::as_glycan_structure()]
#' for how to create such objects manually.
#'
#' @param exp An [glyexp::experiment()] object.
#'
#' @returns The experiment object with the new columns added.
#' @seealso [glymotif::describe_n_glycans()]
#' @export
#'
#' @importFrom rlang .data
add_struct_descriptions <- function(exp) {
  # Check arguments
  checkmate::assert_class(exp, "glyexp_experiment")
  if (!.has_structure_column(exp)) {
    cli::cli_abort("Column {.field glycan_structure} not found in {.field var_info}.")
  }
  if (is.null(exp$meta_data$glycan_type)) {
    cli::cli_abort("Column {.field glycan_type} not found in {.field meta_data}.")
  }
  if (exp$meta_data$glycan_type != "N") {
    cli::cli_abort("Only N-glycans are currently supported.")
  }
  checkmate::assert_class(exp$var_info$glycan_structure, "glyrepr_structure")

  # Check if function has been called before
  if (!is.null(exp$meta_data$struct_descriptions_added)) {
    cli::cli_alert_info("Structure descriptions already added. Skipping.")
    return(exp)
  }

  # Add descriptions (only N-glycans are supported)
  glycan_descriptions <- glymotif::describe_n_glycans(exp$var_info$glycan_structure)
  new_var_info <- dplyr::bind_cols(exp$var_info, glycan_descriptions)
  exp$var_info <- new_var_info
  exp$meta_data$struct_descriptions_added <- TRUE
  exp
}

#' Add Glycan Descriptions
#'
#' This function adds glycan description columns to the
#' variable information tibble of an [glyexp::experiment()] object.
#' If structure information is available,
#' both composition and structure descriptions are added.
#' Otherwise, only composition descriptions are added.
#'
#' @details
#' This function is a wrapper around [add_comp_descriptions()]
#' and [add_struct_descriptions()].
#'
#' @param exp An [glyexp::experiment()] object.
#'
#' @returns The experiment object with the new columns added.
#' @seealso [add_comp_descriptions()], [add_struct_descriptions()]
#' @export
add_glycan_descriptions <- function(exp) {
  # Check if descriptions have been added
  has_comp_desc <- !is.null(exp$meta_data$comp_descriptions_added) && exp$meta_data$comp_descriptions_added
  has_struct_desc <- !is.null(exp$meta_data$struct_descriptions_added) && exp$meta_data$struct_descriptions_added

  if (.has_structure_column(exp) && !has_struct_desc) {
    exp <- add_struct_descriptions(exp)
    cli::cli_alert_success("Structure descriptions added.")
  } else if (.has_structure_column(exp)) {
    cli::cli_alert_info("Structure descriptions already added. Skipping.")
  }

  if (!has_comp_desc) {
    exp <- add_comp_descriptions(exp)
    cli::cli_alert_success("Composition descriptions added.")
  } else {
    cli::cli_alert_info("Composition descriptions already added. Skipping.")
  }

  exp
}

.has_structure_column <- function(exp) {
  "glycan_structure" %in% colnames(exp$var_info)
}

.has_composition_column <- function(exp) {
  "glycan_composition" %in% colnames(exp$var_info)
}