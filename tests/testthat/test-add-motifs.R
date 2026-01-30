test_that("add_motifs_int works", {
  # Create an experiment
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Add motifs
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  suppressWarnings(exp <- add_motifs_int(exp, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(`O-Glycan core 1` = c(1L, 0L), `O-Glycan core 2` = c(0L, 1L))
  expect_identical(exp$var_info[, c("O-Glycan core 1", "O-Glycan core 2")], expected)
})

test_that("add_motifs_lgl works", {
  # Create an experiment
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Add motifs
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  suppressWarnings(exp <- add_motifs_lgl(exp, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(`O-Glycan core 1` = c(TRUE, FALSE), `O-Glycan core 2` = c(FALSE, TRUE))
  expect_identical(exp$var_info[, c("O-Glycan core 1", "O-Glycan core 2")], expected)
})

test_that("add_motifs_int works with named motif names", {
  # Create an experiment
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Add motifs
  motifs <- c(motif1 = "O-Glycan core 1", motif2 = "O-Glycan core 2")
  suppressWarnings(exp <- add_motifs_int(exp, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(motif1 = c(1L, 0L), motif2 = c(0L, 1L))
  expect_identical(exp$var_info[, c("motif1", "motif2")], expected)
})

test_that("add_motifs_int works with named IUPAC strings", {
  # Create an experiment
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Add motifs
  motifs <- c(motif1 = "Gal(b1-3)GalNAc(a1-", motif2 = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  suppressWarnings(exp <- add_motifs_int(exp, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(motif1 = c(1L, 0L), motif2 = c(0L, 1L))
  expect_identical(exp$var_info[, c("motif1", "motif2")], expected)
})

test_that("add_motifs_lgl works for data frames", {
  df <- tibble::tibble(glycan_structure = c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2()))
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  suppressWarnings(df <- add_motifs_lgl(df, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(`O-Glycan core 1` = c(TRUE, FALSE), `O-Glycan core 2` = c(FALSE, TRUE))
  expect_identical(df[, c("O-Glycan core 1", "O-Glycan core 2")], expected)
})

test_that("add_motifs_int works for data frames", {
  df <- tibble::tibble(glycan_structure = c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2()))
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  suppressWarnings(df <- add_motifs_int(df, motifs, alignments = c("whole", "whole")))

  # Check the results
  expected <- tibble::tibble(`O-Glycan core 1` = c(1L, 0L), `O-Glycan core 2` = c(0L, 1L))
  expect_identical(df[, c("O-Glycan core 1", "O-Glycan core 2")], expected)
})

# Tests for new naming behavior

test_that("add_motifs_int uses IUPAC strings for unnamed glyrepr_structure motifs", {
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Use unnamed glyrepr_structure motifs
  motifs <- glyparse::parse_iupac_condensed(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  ))
  suppressWarnings(exp <- add_motifs_int(exp, motifs, alignments = c("whole", "whole")))

  # Check that IUPAC strings are used as column names
  expected <- tibble::tibble(
    `Gal(b1-3)GalNAc(a1-` = c(1L, 0L),
    `Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-` = c(0L, 1L)
  )
  expect_identical(exp$var_info[, c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")], expected)
})

test_that("add_motifs_lgl uses IUPAC strings for unnamed IUPAC string motifs", {
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Use unnamed IUPAC string motifs (not known motif names)
  motifs <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  suppressWarnings(exp <- add_motifs_lgl(exp, motifs, alignments = c("whole", "whole")))

  # Check that IUPAC strings are used as column names
  expected <- tibble::tibble(
    `Gal(b1-3)GalNAc(a1-` = c(TRUE, FALSE),
    `Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-` = c(FALSE, TRUE)
  )
  expect_identical(exp$var_info[, c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")], expected)
})

test_that("add_motifs_int for data frame uses IUPAC strings for unnamed glyrepr_structure motifs", {
  df <- tibble::tibble(glycan_structure = c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2()))

  # Use unnamed glyrepr_structure motifs
  motifs <- glyparse::parse_iupac_condensed(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  ))
  suppressWarnings(df <- add_motifs_int(df, motifs, alignments = c("whole", "whole")))

  # Check that IUPAC strings are used as column names
  expected <- tibble::tibble(
    `Gal(b1-3)GalNAc(a1-` = c(1L, 0L),
    `Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-` = c(0L, 1L)
  )
  expect_identical(df[, c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")], expected)
})

test_that("add_motifs_lgl for data frame uses explicit names when provided", {
  df <- tibble::tibble(glycan_structure = c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2()))

  # Use named IUPAC string motifs
  motifs <- c(motif_a = "Gal(b1-3)GalNAc(a1-", motif_b = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  suppressWarnings(df <- add_motifs_lgl(df, motifs, alignments = c("whole", "whole")))

  # Check that explicit names are used as column names
  expected <- tibble::tibble(motif_a = c(TRUE, FALSE), motif_b = c(FALSE, TRUE))
  expect_identical(df[, c("motif_a", "motif_b")], expected)
})

test_that("add_motifs_int uses motif names for unnamed known motif names", {
  exp <- create_test_exp(c("S1", "S2"), c("V1", "V2"))
  structures <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  exp$var_info$glycan_structure <- structures

  # Use unnamed known motif names
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  suppressWarnings(exp <- add_motifs_int(exp, motifs, alignments = c("whole", "whole")))

  # Check that motif names are used as column names
  expected <- tibble::tibble(
    `O-Glycan core 1` = c(1L, 0L),
    `O-Glycan core 2` = c(0L, 1L)
  )
  expect_identical(exp$var_info[, c("O-Glycan core 1", "O-Glycan core 2")], expected)
})