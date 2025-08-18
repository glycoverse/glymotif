make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  if (mono_type == "generic") {
    glycan <- glyrepr::convert_to_generic(glycan)
  }
  if (!linkage) {
    glycan <- glyrepr::remove_linkages(glycan)
  }
  glycan
}

create_expr_mat <- function(samples, variables) {
  n_row <- length(variables)
  n_col <- length(samples)
  mat <- matrix(1:(n_row*n_col), nrow = n_row)
  rownames(mat) <- variables
  colnames(mat) <- samples
  mat
}

create_sample_info <- function(samples) {
  tibble::tibble(
    sample = samples,
    group = rep("A", length(samples))
  )
}

create_var_info <- function(variables) {
  tibble::tibble(
    variable = variables,
    type = rep("B", length(variables))
  )
}

create_test_exp <- function(samples, variables, exp_type = "glycomics", glycan_type = "N") {
  expr_mat <- create_expr_mat(samples, variables)
  sample_info <- create_sample_info(samples)
  var_info <- create_var_info(variables)
  glyexp::experiment(expr_mat, sample_info, var_info, exp_type, glycan_type)
}