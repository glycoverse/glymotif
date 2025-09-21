#' Quantify motifs in an experiment
#'
#' @description
#' This function quantifies motifs from glycomic or glycoproteomic profiles.
#' For glycomics data, it calculates the motif quantifications directly.
#' For glycoproteomics data, each glycosite is treated as a separate glycome,
#' and motif quantifications are calculated in a site-specific manner.
#'
#' The function takes a `glyexp::experiment()` and returns a new `glyexp::experiment()` 
#' with motif quantifications. Instead of containing the quantification of each glycan 
#' on each glycosite in each sample, the new experiment contains the quantification 
#' of each motif on each glycosite in each sample (for glycoproteomics data) or 
#' the motif quantification in each sample (for glycomics data).
#'
#' Due to the unified data structure of `glyexp::experiment()`,
#' the returned motif experiment can be passed to downstream `glycoverse`
#' functions like `glystats::gly_ttest()` for further statistical analysis.
#' Also, you can use `as_tibble()` to convert it to a "tidy" tibble for custom analysis.
#'
#' @param exp A [glyexp::experiment()] object. Before using this function,
#'   you should preprocess the data using the `glyclean` package.
#'   For glycoproteomics data, the data should be aggregated to the
#'   "gfs" (glycoforms with structures) level using `glyclean::aggregate()`.
#'   Also, please make sure that the `glycan_structure` column is present in the `var_info` table,
#'   as not all glycoproteomics identification softwares provide this information.
#'   The column can be a [glyrepr::glycan_structure()] vector,
#'   or a character vector of glycan structure strings supported by [glyparse::auto_parse()].
#'
#'   For glycoproteomics data, the `var_info` table must contain:
#'   - `protein`: protein ID
#'   - `protein_site`: the glycosite position on the protein
#'   The unique combination of `protein` and `protein_site` determines a glycosite.
#' @inheritParams have_motif
#'
#' @returns A new [glyexp::experiment()] object for motif quantifications.
#'   The new experiment contains the following columns in the `var_info` table:
#'   - `variable`: variable ID
#'   - `motif`: motif name
#'
#'   For glycoproteomics data, with additional columns:
#'   - `protein`: protein ID
#'   - `protein_site`: the glycosite position on the protein
#'
#'   Other columns in the `var_info` table (e.g. `gene`) are retained if they have "many-to-one"
#'   relationship with glycosites (unique combinations of `protein`, `protein_site`).
#'   That is, each glycosite cannot have multiple values for these columns.
#'   `gene` is a common example, as a glycosite can only be associate with one gene.
#'   Descriptions about glycans are not such a column, as a glycosite can have multiple glycans,
#'   thus having multiple descriptions.
#'   Columns not having this relationship with glycosites will be dropped.
#'   Don't worry if you cannot understand this logic,
#'   as long as you know that this function will try its best to preserve useful information.
#'
#'   `sample_info` and `meta_data` are not modified,
#'   except that the `exp_type` field of `meta_data` is set to "traitomics" for glycomics data,
#'   and "traitproteomics" for glycoproteomics data.
#'
#' @seealso [count_motifs()]
#' @export
quantify_motifs <- function(exp, motifs, alignments = NULL, ignore_linkages = FALSE) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort("The experiment must have a {.field glycan_structure} column.")
  }

  # Extract glycan structures
  glycan_structures <- ensure_glycans_are_structures(exp$var_info$glycan_structure)

  # Calculate motif counts for each glycan
  motif_counts_matrix <- count_motifs(
    glycans = glycan_structures,
    motifs = motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages
  )

  # Get motif names
  motif_names <- colnames(motif_counts_matrix)

  if (exp$meta_data$exp_type == "glycomics") {
    # Case: glycomics data
    result_data <- .quantify_glycomic_motifs(exp, motif_counts_matrix, motif_names)
  } else {
    # Case: glycoproteomics data
    result_data <- .quantify_glycoprot_motifs(exp, motif_counts_matrix, motif_names)
  }

  new_var_info <- result_data$var_info
  new_expr_mat <- result_data$expr_mat

  # Create new experiment with updated meta_data
  new_meta_data <- exp$meta_data
  new_meta_data$exp_type <- if (exp$meta_data$exp_type == "glycomics") "traitomics" else "traitproteomics"

  # Use new_experiment to bypass validation since "traitomics" is not a standard exp_type
  result <- glyexp:::new_experiment(
    expr_mat = new_expr_mat,
    sample_info = exp$sample_info,
    var_info = new_var_info,
    meta_data = new_meta_data
  )

  result
}


# Helper function for quantifying motifs in glycomics data
.quantify_glycomic_motifs <- function(exp, motif_counts_matrix, motif_names) {
  n_motifs <- length(motif_names)

  # Create new var_info with only motif column
  new_var_info <- tibble::tibble(
    variable = paste0("V", seq_len(n_motifs)),
    motif = motif_names
  )

  # Create new expression matrix by summing across all glycans for each motif
  n_samples <- ncol(exp$expr_mat)
  new_expr_mat <- matrix(0, nrow = n_motifs, ncol = n_samples)
  colnames(new_expr_mat) <- colnames(exp$expr_mat)
  rownames(new_expr_mat) <- new_var_info$variable

  # Sum expression values weighted by motif counts for each motif
  for (j in seq_len(n_motifs)) {
    motif_weighted_values <- colSums(exp$expr_mat * motif_counts_matrix[, j])
    new_expr_mat[j, ] <- motif_weighted_values
  }

  list(var_info = new_var_info, expr_mat = new_expr_mat)
}


# Helper function for quantifying motifs in glycoproteomics data
.quantify_glycoprot_motifs <- function(exp, motif_counts_matrix, motif_names) {
  # Check required columns for glycoproteomics data
  if (!all(c("protein", "protein_site") %in% colnames(exp$var_info))) {
    cli::cli_abort("For glycoproteomics data, {.field protein} and {.field protein_site} columns are required.")
  }

  n_motifs <- length(motif_names)

  # Create glycosite identifiers
  glycosites <- stringr::str_c(exp$var_info[["protein"]], exp$var_info[["protein_site"]], sep = "@")
  splits <- split(seq_along(glycosites), glycosites)

  # Function to quantify motifs for one glycosite
  quantify_one_site <- function(site_idx) {
    site_expr_mat <- exp$expr_mat[site_idx, , drop = FALSE]
    site_motif_counts <- motif_counts_matrix[site_idx, , drop = FALSE]

    # Create result matrix for this site (one row per motif)
    site_res_mat <- matrix(0, nrow = n_motifs, ncol = ncol(site_expr_mat))
    colnames(site_res_mat) <- colnames(site_expr_mat)

    for (j in seq_len(n_motifs)) {
      # Sum expression values weighted by motif counts for this motif
      aggregated_value <- colSums(site_expr_mat * site_motif_counts[, j])
      site_res_mat[j, ] <- aggregated_value
    }

    # Create var_info for this site
    site_var_info <- tibble::tibble(
      variable = NA_character_,  # Placeholder, will be replaced later
      protein = unique(exp$var_info[["protein"]][site_idx]),
      protein_site = unique(exp$var_info[["protein_site"]][site_idx]),
      motif = motif_names
    )

    list(res_mat = site_res_mat, res_var_info = site_var_info)
  }

  # Process all glycosites
  res_list <- purrr::map(splits, quantify_one_site)

  # Combine results
  new_expr_mat <- do.call(rbind, purrr::map(res_list, ~ .x$res_mat))
  new_var_info <- do.call(rbind, purrr::map(res_list, ~ .x$res_var_info))

  # Set variable names
  new_var_info$variable <- paste0("V", seq_len(nrow(new_var_info)))
  rownames(new_expr_mat) <- new_var_info$variable

  # Get glycosite descriptive columns (many-to-one relationship with glycosites)
  glycosite_descriptions <- .get_glycosite_var_info(exp$var_info)
  new_var_info <- dplyr::left_join(new_var_info, glycosite_descriptions, by = c("protein", "protein_site"))

  list(var_info = new_var_info, expr_mat = new_expr_mat)
}

#' Get Glycosite Descriptive Columns for quantify_motifs
#'
#' This function is adapted from glydet's .get_glycosite_var_info function.
#' Some columns are descriptive columns of glycosites.
#' That is, for each combination of `protein` and `protein_site`,
#' there is only one value for these columns.
#' A common example is `gene`.
#' This function returns the distinct values of these columns,
#' along with `protein` and `protein_site`.
#' It is used to be left joined with the resulting var_info table in `quantify_motifs()`.
#'
#' @param var_info A tibble with the variable information.
#'
#' @returns A tibble with the distinct values of the descriptive columns,
#'   along with `protein` and `protein_site`.
#'
#' @noRd
.get_glycosite_var_info <- function(var_info) {
  cols <- var_info |>
    dplyr::group_by(.data$protein, .data$protein_site) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::n_distinct)) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::where(~ all(.x == 1))) |>
    colnames()

  var_info |>
    dplyr::select(dplyr::all_of(c("protein", "protein_site", cols))) |>
    dplyr::distinct()
}
