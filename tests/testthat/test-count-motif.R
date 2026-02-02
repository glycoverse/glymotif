test_that("count motifs in glycan", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-"
  expect_equal(count_motif(glycan, motif), 2L)
})


test_that("count_motif supports multiple glycan structure formats", {
  # Test IUPAC-condensed format
  glycan_iupac <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif_iupac <- "Gal(b1-"
  expect_equal(count_motif(glycan_iupac, motif_iupac), 2L)

  # Test IUPAC-short format
  glycan_short <- "Galb3Galb3GalNAca-"
  motif_short <- "Galb-"
  expect_equal(count_motif(glycan_short, motif_short), 2L)
})


test_that("count motifs with branching", {
  glycan <- "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-"
  motif <- "Man(b1-?)[Man(b1-?)]GalNAc(b1-"
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count symmetrical motif", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-3)Gal(b1-"
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count 0 motif", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-4)GalNAc(b1-"
  expect_equal(count_motif(glycan, motif), 0L)
})


# ========== count_motifs ==========
test_that("count_motifs works with multiple motifs", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  names(glycans) <- c("double_gal", "single_gal")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  # IUPAC strings without names should have no colnames
  expect_null(colnames(result))
  expect_equal(rownames(result), c("double_gal", "single_gal"))
  
  # Check specific counts
  expect_equal(as.integer(result[1, 1]), 1L)  # Gal(b1-3)GalNAc(b1-
  expect_equal(as.integer(result[1, 2]), 2L)  # Gal(b1-
  expect_equal(as.integer(result[2, 1]), 1L)  # Gal(b1-3)GalNAc(b1-
  expect_equal(as.integer(result[2, 2]), 1L)  # Gal(b1-
})


test_that("count_motifs works without glycan names", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  # IUPAC strings without names should have no colnames
  expect_null(colnames(result))
  
  # Check that row names contain IUPAC strings
  expect_true(grepl("Gal\\(b1-3\\)", rownames(result)[1]))
  expect_true(grepl("Gal\\(b1-3\\)", rownames(result)[2]))
})


test_that("count_motifs works with complex branching motifs", {
  glycans <- c("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-", "Gal(b1-3)Gal(b1-")
  names(glycans) <- c("complex", "simple")
  
  motifs <- c("Man(b1-?)[Man(b1-?)]GalNAc(b1-", "Man(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  
  # Check complex structure counts
  expect_equal(as.integer(result[1, 1]), 1L)  # Man(b1-?)[Man(b1-?)]GalNAc(b1-
  expect_equal(as.integer(result[1, 2]), 2L)  # Man(b1-
  expect_equal(as.integer(result[1, 3]), 0L)  # Gal(b1-
  
  # Check simple structure counts
  expect_equal(as.integer(result[2, 1]), 0L)  # Man(b1-?)[Man(b1-?)]GalNAc(b1-
  expect_equal(as.integer(result[2, 2]), 0L)  # Man(b1-
  expect_equal(as.integer(result[2, 3]), 2L)  # Gal(b1-
})


test_that("count_motifs works with different alignments", {
  glycans <- c("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(a1-", "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(a1-")
  names(glycans) <- c("glycan1", "glycan2")
  
  motifs <- c("Gal(a1-3)Gal(a1-4)Gal(a1-", "Gal(a1-4)Gal(a1-6)Gal(a1-")
  alignments <- c("core", "terminal")
  
  result <- count_motifs(glycans, motifs, alignments = alignments)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("count_motifs handles zero counts", {
  glycans <- c("Gal(b1-3)GalNAc(b1-", "GlcNAc(b1-4)GlcNAc(b1-")
  
  motifs <- c("Man(b1-", "Fuc(a1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_true(all(result[, 1] == 0L))  # Man(b1-
  expect_true(all(result[, 2] == 0L))  # Fuc(a1-
})


test_that("count_motifs handles empty motifs", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- character(0)
  
  expect_error(count_motifs(glycans, motifs), "`motifs` cannot be empty")
})

test_that("count_motifs raises error for duplicate motifs", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- c("Gal(b1-", "Gal(b1-")

  expect_error(
    count_motifs(glycans, motifs),
    "cannot have duplications"
  )
})


test_that("count_motifs handles invalid motifs argument", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- 123

  expect_error(count_motifs(glycans, motifs), "`motifs` must be a 'glyrepr_structure' object")
})


test_that("count_motifs handles mismatched alignments length", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-", "GalNAc(b1-")
  alignments <- c("substructure", "core")  # length 2, motifs length 3
  
  expect_error(count_motifs(glycans, motifs, alignments = alignments), 
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})


test_that("count_motifs works with single alignment for all motifs", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs, alignments = "substructure")
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("count_motifs works with ignore_linkages", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-4)GalNAc(b1-")
  
  result_normal <- count_motifs(glycans, motifs, ignore_linkages = FALSE)
  result_ignore <- count_motifs(glycans, motifs, ignore_linkages = TRUE)
  
  # Use numeric indices since IUPAC strings without names have no colnames
  expect_equal(as.integer(result_normal[1, 2]), 0L)  # Gal(b1-4)GalNAc(b1-
  expect_equal(as.integer(result_ignore[1, 2]), 1L)  # Gal(b1-4)GalNAc(b1-
})


# ========== Name preservation tests ==========
test_that("count_motif preserves names from glycan_structure input", {
  glycan1 <- glyrepr::o_glycan_core_1()
  glycan2 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("core1", "core2")

  result <- count_motif(glycans, "Gal(b1-")
  expect_equal(names(result), c("core1", "core2"))
})

test_that("count_motif preserves names from character vector input", {
  glycans <- c("Gal(b1-3)GalNAc(?1-", "Gal(b1-3)Gal(b1-3)GalNAc(?1-")
  names(glycans) <- c("simple", "complex")

  result <- count_motif(glycans, "Gal(b1-")
  expect_equal(names(result), c("simple", "complex"))
})

test_that("count_motif returns no names when input has no names", {
  glycans <- c("Gal(b1-3)GalNAc(?1-", "Gal(b1-3)Gal(b1-3)GalNAc(?1-")

  result <- count_motif(glycans, "Gal(b1-")
  expect_null(names(result))
})

# ========== match_degree ==========
test_that("count_motif respects match_degree", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  expect_equal(count_motif(glycan, motif), 1L)
  expect_equal(count_motif(glycan, motif, match_degree = c(FALSE, TRUE)), 0L)
})

test_that("count_motif differentiates branching GlcNAc and bisecting GlcNAc with match_degree", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
  motif <- "GlcNAc(b1-2)Man(a1-3/6)Man(b1-"

  expect_equal(count_motif(glycan, motif, match_degree = c(TRUE, FALSE, FALSE)), 2L)
})

test_that("count_motifs validates match_degree list", {
  glycans <- glyrepr::o_glycan_core_2()
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"))

  expect_error(count_motifs(glycans, motifs, match_degree = list(c(TRUE, FALSE))), "match_degree.*same length")
})

# ========== Known Motif Names as Column Names ==========
test_that("count_motifs uses known motif names as colnames when motifs unnamed", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  names(glycans) <- c("core1", "core2")

  # Use known motif names without explicit names
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")

  result <- count_motifs(glycans, motifs)
  expect_equal(rownames(result), c("core1", "core2"))
  expect_equal(colnames(result), c("O-Glycan core 1", "O-Glycan core 2"))
})

test_that("count_motifs preserves explicit motif names over known names", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  names(glycans) <- c("core1", "core2")

  # Use known motif names WITH explicit names
  motifs <- c(m1 = "O-Glycan core 1", m2 = "O-Glycan core 2")

  result <- count_motifs(glycans, motifs)
  expect_equal(colnames(result), c("m1", "m2"))
})

test_that("count_motifs has no colnames when motifs are IUPAC strings without names", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  names(glycans) <- c("core1", "core2")

  motifs <- c("Gal(b1-3)GalNAc(?1-", "Gal(b1-4)GalNAc(?1-")

  result <- count_motifs(glycans, motifs)
  expect_null(colnames(result))
})

# ========== Integration tests with motif specs ==========
test_that("count_motifs works with dynamic_motifs()", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))

  result <- count_motifs(glycans, dynamic_motifs(max_size = 2))

  expect_type(result, "integer")
  expect_equal(dim(result), c(2, length(extract_motif(glycans, max_size = 2))))
})

test_that("count_motifs works with branch_motifs()", {
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )

  result <- count_motifs(glycans, branch_motifs())

  expect_type(result, "integer")
  expect_equal(nrow(result), 1)
  expect_true(ncol(result) > 0)
})

# ========== IUPAC Column Names for dynamic_motifs and branch_motifs ==========
test_that("count_motifs returns IUPAC strings as colnames for dynamic_motifs", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))
  result <- count_motifs(glycans, dynamic_motifs(max_size = 2))
  
  # Column names should be the IUPAC strings of extracted motifs
  expect_type(colnames(result), "character")
  expect_equal(length(colnames(result)), ncol(result))
  # All column names should be valid IUPAC strings (non-empty)
  expect_true(all(nchar(colnames(result)) > 0))
})

test_that("count_motifs returns trimmed IUPAC strings as colnames for branch_motifs", {
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  result <- count_motifs(glycans, branch_motifs())
  
  # Column names should be present
  expect_type(colnames(result), "character")
  expect_equal(length(colnames(result)), ncol(result))
  expect_true(all(nchar(colnames(result)) > 0))
  
  # Column names should NOT contain the core suffix
  expect_false(any(grepl(")Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc", colnames(result), fixed = TRUE)))
  expect_false(any(grepl(")Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc", colnames(result), fixed = TRUE)))
  
  # Column names should end with the branch root linkage pattern
  expect_true(all(grepl("\\([a-z]1-$", colnames(result))))
})
