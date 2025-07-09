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

test_that("quantify_motifs correctly handles complex glycoproteomics data with tibble comparisons", {
  # This test specifically targets the bug where tibble comparisons failed
  # in the aggregation logic, causing most motif values to be zero

  # ----- Create a more complex experiment -----
  expr_mat <- matrix(
    c(100, 200, 150,
      300, 400, 350,
      500, 600, 550,
      700, 800, 750,
      900, 1000, 950),
    nrow = 5, byrow = TRUE
  )
  colnames(expr_mat) <- c("S1", "S2", "S3")
  rownames(expr_mat) <- c("V1", "V2", "V3", "V4", "V5")

  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))

  # Use realistic protein IDs and gene names to test complex string comparisons
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3", "V4", "V5"),
    protein = c("P08185", "P04196", "P04196", "P10909", "P13671"),
    gene = c("SERPINA6", "HRG", "HRG", "CLU", "C6"),
    protein_site = c(176L, 344L, 345L, 291L, 855L),
    glycan_structure = glyrepr::as_glycan_structure(c(
      "NeuAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",
      "Hex(??-?)HexNAc(??-?)HexNAc(??-",
      "HexNAc(??-?)HexNAc(??-",
      "NeuAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",
      "dHex(??-?)HexNAc(??-?)HexNAc(??-"
    ))
  )

  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N"
  )

  # ----- Test the function -----
  motifs <- c(
    "NeuAc(??-?)Hex(??-?)HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-"
  )

  result <- quantify_motifs(exp, motifs, alignments = "terminal", ignore_linkages = TRUE)

  # ----- Check the results -----
  # Each protein site should have 3 motif variables (one per motif)
  expect_equal(nrow(result$var_info), 5 * 3)  # 5 unique sites × 3 motifs
  expect_equal(ncol(result$expr_mat), 3)  # 3 samples

  # Check that we have the expected number of non-zero values
  # Based on motif counts: V1 and V4 have motif1, V2 has motif2, V3 has motif3
  expected_non_zero_count <- 4 * 3  # 4 non-zero motif matches × 3 samples
  actual_non_zero_count <- sum(result$expr_mat != 0)
  expect_equal(actual_non_zero_count, expected_non_zero_count)

  # Check specific aggregated values
  # V1 (P08185_176) should have motif1 = 100,200,150 and others = 0
  v1_motif1_idx <- which(result$var_info$protein == "P08185" &
                         result$var_info$protein_site == 176 &
                         result$var_info$motif == motifs[1])
  expect_equal(as.numeric(result$expr_mat[v1_motif1_idx, ]), c(100, 200, 150))

  # V2 (P04196_344) should have motif2 = 300,400,350 and others = 0
  v2_motif2_idx <- which(result$var_info$protein == "P04196" &
                         result$var_info$protein_site == 344 &
                         result$var_info$motif == motifs[2])
  expect_equal(as.numeric(result$expr_mat[v2_motif2_idx, ]), c(300, 400, 350))

  # V3 (P04196_345) should have motif3 = 500,600,550 and others = 0
  v3_motif3_idx <- which(result$var_info$protein == "P04196" &
                         result$var_info$protein_site == 345 &
                         result$var_info$motif == motifs[3])
  expect_equal(as.numeric(result$expr_mat[v3_motif3_idx, ]), c(500, 600, 550))

  # V4 (P10909_291) should have motif1 = 700,800,750 and others = 0
  v4_motif1_idx <- which(result$var_info$protein == "P10909" &
                         result$var_info$protein_site == 291 &
                         result$var_info$motif == motifs[1])
  expect_equal(as.numeric(result$expr_mat[v4_motif1_idx, ]), c(700, 800, 750))

  # Check that the aggregation logic works correctly by ensuring
  # the total expression is preserved for each motif
  structures <- exp$var_info$glycan_structure
  motif_counts_matrix <- count_motifs(structures, motifs, alignments = "terminal", ignore_linkages = TRUE)

  for (j in 1:length(motifs)) {
    # Calculate expected total for this motif across all samples
    expected_total <- colSums(exp$expr_mat * motif_counts_matrix[, j])

    # Get actual total from result
    motif_rows <- which(result$var_info$motif == motifs[j])
    actual_total <- colSums(result$expr_mat[motif_rows, , drop = FALSE])

    expect_equal(actual_total, expected_total,
                 info = paste("Motif", j, "total expression mismatch"))
  }
})

test_that("quantify_motifs correctly aggregates multiple glycopeptides per protein site", {
  # This test ensures that when multiple glycopeptides belong to the same protein site,
  # they are correctly aggregated according to their motif counts

  # ----- Create experiment with multiple glycopeptides per site -----
  expr_mat <- matrix(
    c(100, 200, 150,   # V1: P08185_176, motif1 count = 2
      50, 100, 75,     # V2: P08185_176, motif1 count = 1
      300, 400, 350,   # V3: P04196_344, motif2 count = 1
      200, 300, 250),  # V4: P04196_344, motif2 count = 2
    nrow = 4, byrow = TRUE
  )
  colnames(expr_mat) <- c("S1", "S2", "S3")
  rownames(expr_mat) <- c("V1", "V2", "V3", "V4")

  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))

  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3", "V4"),
    protein = c("P08185", "P08185", "P04196", "P04196"),  # Two sites, each with 2 glycopeptides
    gene = c("SERPINA6", "SERPINA6", "HRG", "HRG"),
    protein_site = c(176L, 176L, 344L, 344L),  # Same sites repeated
    glycan_structure = glyrepr::as_glycan_structure(c(
      "NeuAc(??-?)Hex(??-?)HexNAc(??-?)NeuAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",  # 2 motif1
      "NeuAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",                                  # 1 motif1
      "Hex(??-?)HexNAc(??-?)HexNAc(??-",                                             # 1 motif2
      "Hex(??-?)HexNAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"                        # 2 motif2
    ))
  )

  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N"
  )

  # ----- Test the function -----
  motifs <- c("NeuAc(??-?)Hex(??-?)HexNAc(??-", "Hex(??-?)HexNAc(??-")
  result <- quantify_motifs(exp, motifs, alignments = "terminal", ignore_linkages = TRUE)

  # ----- Check the results -----
  # Should have 2 unique sites × 2 motifs = 4 variables
  expect_equal(nrow(result$var_info), 4)

  # Check aggregation for P08185_176 motif1:
  # V1 contributes: 100*1, 200*1, 150*1 = 100, 200, 150
  # V2 contributes: 50*1, 100*1, 75*1 = 50, 100, 75
  # Total: 150, 300, 225
  p08185_motif1_idx <- which(result$var_info$protein == "P08185" &
                             result$var_info$protein_site == 176 &
                             result$var_info$motif == motifs[1])
  expect_equal(as.numeric(result$expr_mat[p08185_motif1_idx, ]), c(150, 300, 225))

  # Check aggregation for P04196_344 motif2:
  # V3 contributes: 300*1, 400*1, 350*1 = 300, 400, 350
  # V4 contributes: 200*1, 300*1, 250*1 = 200, 300, 250
  # Total: 500, 700, 600
  p04196_motif2_idx <- which(result$var_info$protein == "P04196" &
                             result$var_info$protein_site == 344 &
                             result$var_info$motif == motifs[2])
  expect_equal(as.numeric(result$expr_mat[p04196_motif2_idx, ]), c(500, 700, 600))

  # Check that motifs with zero counts result in zero values
  p08185_motif2_idx <- which(result$var_info$protein == "P08185" &
                             result$var_info$protein_site == 176 &
                             result$var_info$motif == motifs[2])
  expect_equal(as.numeric(result$expr_mat[p08185_motif2_idx, ]), c(0, 0, 0))

  p04196_motif1_idx <- which(result$var_info$protein == "P04196" &
                             result$var_info$protein_site == 344 &
                             result$var_info$motif == motifs[1])
  expect_equal(as.numeric(result$expr_mat[p04196_motif1_idx, ]), c(0, 0, 0))
})
