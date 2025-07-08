test_that("quantify_motifs works in a most common glycoproteomics situation", {
  # ----- Create an experiment for test -----
  # Create expression matrix
  expr_mat <- matrix(
    c(1, 2, 1,
      2, 1, 2,
      1, 2, 1),
    nrow = 3, byrow = TRUE
  )
  colnames(expr_mat) <- c("S1", "S2", "S3")
  rownames(expr_mat) <- c("V1", "V2", "V3")

  # Create sample info
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))

  # Create variable info
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    protein = c("PRO1", "PRO2", "PRO2"),  # V2 and V3 are the same protein and gene
    gene = c("GENE1", "GENE2", "GENE2"),
    protein_site = c(10L, 20L, 20L),
    glycan_composition = glyrepr::as_glycan_composition(c(
      "GalNAc(1)", "Gal(1)GalNAc(1)", "Gal(1)GalNAc(1)Fuc(1)"
    )),
    glycan_structure = glyrepr::as_glycan_structure(c(
      "GalNAc(b1-",
      "Gal(b1-3)GalNAc(b1-",
      "Gal(b1-3)[Fuc(a1-3)]GalNAc(b1-"
    ))
  )

  # Create the experiment
  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N"
  )

  # ----- Test the function -----
  #     Motif 1    Motif 2
  # V1     0          0
  # V2     1          0
  # V3     1          1
  motifs <- c(motif1 = "Gal(b1-3)GalNAc(b1-", motif2 = "Fuc(a1-3)GalNAc(b1-")
  result <- quantify_motifs(exp, motifs, alignments = "substructure")

  # ----- Check the results -----
  # Check var_info
  expected_var_info <- tibble::tibble(
    variable = paste0("V", 1:4),
    protein = rep(c("PRO1", "PRO2"), each = 2),
    gene = rep(c("GENE1", "GENE2"), each = 2),
    protein_site = rep(c(10L, 20L), each = 2),
    motif = rep(c("motif1", "motif2"), times = 2)
  )
  expect_equal(result$var_info, expected_var_info)

  # Check sample_info
  expect_equal(result$sample_info, exp$sample_info)

  # Check expr_mat
  expected_expr_mat <- matrix(
    # S1 S2 S3
    c(0, 0, 0,   # PRO1
      0, 0, 0,   # PRO1
      3, 3, 3,   # PRO2
      1, 2, 1),  # PRO2
    nrow = 4, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- paste0("V", 1:4)
  expect_equal(result$expr_mat, expected_expr_mat)

  # Check experiment type
  expect_equal(result$meta_data$exp_type, "motifomics")
})

test_that("quantify_motifs works in a common glycomics situation", {
  # ----- Create an experiment for test -----
  # Create expression matrix
  expr_mat <- matrix(
    c(1, 2, 1,
      2, 1, 2,
      1, 2, 1),
    nrow = 3, byrow = TRUE
  )
  colnames(expr_mat) <- c("S1", "S2", "S3")
  rownames(expr_mat) <- c("V1", "V2", "V3")

  # Create sample info
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))

  # Create variable info
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    glycan_composition = glyrepr::as_glycan_composition(c(
      "GalNAc(1)", "Gal(1)GalNAc(1)", "Gal(1)GalNAc(1)Fuc(1)"
    )),
    glycan_structure = glyrepr::as_glycan_structure(c(
      "GalNAc(b1-",
      "Gal(b1-3)GalNAc(b1-",
      "Gal(b1-3)[Fuc(a1-3)]GalNAc(b1-"
    ))
  )

  # Create the experiment
  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N"
  )

  # ----- Test the function -----
  #     Motif 1    Motif 2
  # V1     0          0
  # V2     1          0
  # V3     1          1
  motifs <- c(motif1 = "Gal(b1-3)GalNAc(b1-", motif2 = "Fuc(a1-3)GalNAc(b1-")
  result <- quantify_motifs(exp, motifs, alignments = "substructure")

  # ----- Check the results -----
  # Check var_info
  expected_var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    motif = c("motif1", "motif2")
  )
  expect_equal(result$var_info, expected_var_info)

  # Check sample_info
  expect_equal(result$sample_info, exp$sample_info)

  # Check expr_mat
  expected_expr_mat <- matrix(
    c(3, 3, 3,
      1, 2, 1),
    nrow = 2, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- paste0("V", 1:2)
  expect_equal(result$expr_mat, expected_expr_mat)

  # Check experiment type
  expect_equal(result$meta_data$exp_type, "motifomics")
})
