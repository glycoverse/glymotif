#' Quantify motifs in an experiment
#'
#' @description
#' This function takes a `glyexp::experiment()`,
#' and returns a new `glyexp::experiment()` with motif quantifications.
#'
#' The new experiment is different from a normal `glyexp::experiment()` in three ways:
#' 1. It doesn't have a "glycan_structure" or "glycan_composition" column,
#' but a "motif" column instead.
#' 2. The metadata field "exp_type" is "motifomics", not "glycoproteomics" or glycomics.
#'
#' You can understand a "motif experiment" in this way:
#' For glycoproteomics data,
#' instead of containing the quantification of each glycan on each glycosite,
#' it now contains the quantification of each motif on each glycosite in each sample.
#' For glycomics data,
#' it contains the motif quantification in each sample.
#'
#' Due to the unified data structure of `glyexp::experiment()`,
#' the returned motif experiment can be passed to downstream `glycoverse`
#' functions like `glystats::gly_ttest()` for further statistical analysis.
#' Also, you can use `as_tibble()` to convert it to a "tidy" tibble for custom analysis.
#'
#' @details
#' Here is a comprehensive comparison of the input experiment and the returned experiment:
#'
#' **Sample information:**
#' Sample information doesn't change.
#'
#' **Variable information:**
#' Firstly, we define the "base variables" as the unique combination of all columns
#' in the variable information tibble except for "variable", "glycan_composition", and "glycan_structure".
#' The new var_info would be a Cartesian product of all "base variables" and all motifs.
#' For glycomics data, no "base variable" is defined,
#' so the resulted var_info only contains the "motif" column.
#'
#' **Expression matrix:**
#' Rows are variables, columns are samples, same as before.
#' Filling the cells are the quantifications for motifs on different glycosites in each sample,
#' instead of quantifications of glycoforms or glycopeptides.
#' Or just motif quantification in each sample, for glycomics data.
#'
#' @param exp A `glyexp::experiment()` object containing glycoproteomics data.
#' @inheritParams have_motif
#'
#' @returns A `glyexp::experiment()` object containing motif quantifications.
#' @seealso [count_motifs()]
#' @export
quantify_motifs <- function(exp, motifs, alignments = NULL, ignore_linkages = FALSE) {
  # Validate inputs
  checkmate::assert_class(exp, "glyexp_experiment")
  if (!"glycan_structure" %in% colnames(exp$var_info)) {
    cli::cli_abort("The experiment must have a {.field glycan_structure} column.")
  }
  checkmate::assert_class(exp$var_info$glycan_structure, "glyrepr_structure")

  # Extract glycan structures
  glycan_structures <- exp$var_info$glycan_structure

  # Calculate motif counts for each glycan
  motif_counts_matrix <- count_motifs(
    glycans = glycan_structures,
    motifs = motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages
  )

  # Create base variables (all columns except variable, glycan_composition, glycan_structure)
  base_var_cols <- setdiff(
    colnames(exp$var_info),
    c("variable", "glycan_composition", "glycan_structure")
  )

  # Get motif names
  motif_names <- colnames(motif_counts_matrix)
  n_motifs <- length(motif_names)

  if (length(base_var_cols) == 0) {
    # Case: glycomics data with no base variables
    result_data <- .quantify_glycomic_motifs(exp, motif_counts_matrix, motif_names)
  } else {
    # Case: glycoproteomics data with base variables
    result_data <- .quantify_glycoprot_motifs(exp, motif_counts_matrix, motif_names, base_var_cols)
  }

  new_var_info <- result_data$var_info
  new_expr_mat <- result_data$expr_mat

  # Create new experiment with updated meta_data
  new_meta_data <- exp$meta_data
  new_meta_data$exp_type <- "motifomics"

  # Use new_experiment to bypass validation since "motifomics" is not a standard exp_type
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
.quantify_glycoprot_motifs <- function(exp, motif_counts_matrix, motif_names, base_var_cols) {
  n_motifs <- length(motif_names)

  # Group by unique combinations of base variables
  base_variables <- exp$var_info[, base_var_cols, drop = FALSE]

  # Find unique base variable combinations
  unique_base_vars <- unique(base_variables)
  n_base_vars <- nrow(unique_base_vars)

  # Create new var_info as Cartesian product of unique base variables and motifs
  new_var_info <- tibble::tibble(
    variable = paste0("V", seq_len(n_base_vars * n_motifs))
  )

  # Add base variable columns
  for (col_name in base_var_cols) {
    new_var_info[[col_name]] <- rep(unique_base_vars[[col_name]], each = n_motifs)
  }

  # Add motif column
  new_var_info$motif <- rep(motif_names, times = n_base_vars)

  # Create new expression matrix
  n_samples <- ncol(exp$expr_mat)
  new_expr_mat <- matrix(0, nrow = n_base_vars * n_motifs, ncol = n_samples)
  colnames(new_expr_mat) <- colnames(exp$expr_mat)
  rownames(new_expr_mat) <- new_var_info$variable

  # Fill the new expression matrix by aggregating within each base variable group
  for (i in seq_len(n_base_vars)) {
    # Find rows that match this unique base variable combination
    # Fixed: Use proper tibble comparison instead of apply
    target_row <- unique_base_vars[i, ]
    matching_conditions <- rep(TRUE, nrow(base_variables))

    for (col_name in base_var_cols) {
      matching_conditions <- matching_conditions &
        (base_variables[[col_name]] == target_row[[col_name]])
    }

    matching_rows <- which(matching_conditions)

    for (j in seq_len(n_motifs)) {
      row_idx <- (i - 1) * n_motifs + j

      # Sum expression values weighted by motif counts for matching rows
      aggregated_value <- colSums(exp$expr_mat[matching_rows, , drop = FALSE] *
                                    motif_counts_matrix[matching_rows, j])
      new_expr_mat[row_idx, ] <- aggregated_value
    }
  }

  list(var_info = new_var_info, expr_mat = new_expr_mat)
}
