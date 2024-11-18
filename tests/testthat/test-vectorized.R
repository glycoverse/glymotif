# ----- has_motifs() -----
test_that("has_motifs: IUPAC-condensed", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(Gal = TRUE, GlcNAc = TRUE, GalNAc = FALSE))
})


test_that("has_motifs: motif names", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  # `glycan` is actually o-glycan core 2, with an additional
  # GlcNAc at the reducing end.
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  # Note that these two motifs both have "core" alignment,
  # so when `alignments` is not provided, it will use the default alignment.
  expect_equal(result, c("O-Glycan core 1" = FALSE, "O-Glycan core 2" = FALSE))
})


test_that("has_motifs: motif names with custom alignments", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-3)GlcNAc"
  # Output should be c("O-Glycan core 1" = TRUE, "O-Glycan core 2" = TRUE)
  expect_snapshot(has_motifs(glycan, motifs, alignments = "substructure"))
})


test_that("has_motifs: custom alignments with wrong length", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycan <- "Gal(b1-3)GlcNAc"
  alignments <- c("substructure", "core")
  expect_snapshot(has_motifs(glycan, motifs, alignments = alignments), error = TRUE)
})


test_that("has_motifs: glycan graphs", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(TRUE, TRUE, FALSE))
})


test_that("has_motifs: some bad IUPAC", {
  motifs <- c("Gal", "bad1", "bad2")
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("has_motifs: some bad motif names", {
  motifs <- c("GM1", "bad1", "bad2")
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("has_motifs: wrong types", {
  motifs <- 1:3
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("has_motifs: empty motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  expected <- logical(0)
  names(expected) <- character(0)
  expect_equal(has_motifs(glycan, character(0)), expected)
})


test_that("has_motifs: motif names are kept for IUPAC-condensed motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  expect_equal(names(has_motifs(glycan, motifs)), names(motifs))
})


test_that("has_motifs: motif names are kept for glycan graph motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  expect_equal(names(has_motifs(glycan, motifs)), names(motifs))
})


test_that("has_motifs: motif names are kept for motif names", {
  glycan <- "Gal(b1-3)GlcNAc"
  motifs <- c(M1 = "O-Glycan core 1", M2 = "O-Glycan core 2")
  expect_equal(names(has_motifs(glycan, motifs)), names(motifs))
})


test_that("has_motifs: missing motifs uses default motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan)
  expect_snapshot(result[1:10])
})


test_that("has_motifs: motifs have different mono types", {
  glycan <- "Gal(b1-3)GlcNAc"
  motifs <- c("Gal", "Hex")
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("has_motifs: glycan has different mono types with motifs", {
  glycan <- "Gal"
  motifs <- c(M1 = "Hex", M2 = "HexNAc")
  expect_equal(has_motifs(glycan, motifs), c(M1 = TRUE, M2 = FALSE))
})


test_that("has_motifs: glycan has obscurer mono type than motifs", {
  glycan <- "Hex"
  motifs <- c("Gal", "GlcNAc")
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("has_motifs: glycan and motifs have different graph modes", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc", mode = "ne")
  motifs <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Gal(b1-4)GlcNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed, mode = "dn")
  expect_equal(has_motifs(glycan, motifs), c(G1 = TRUE, G2 = FALSE))
})


# ----- have_motif() -----
test_that("have_motif: IUPAC-condensed", {
  motif <- "Gal"
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  result <- have_motif(glycans, motif)
  expected <- c("Gal(b1-3)GlcNAc" = TRUE, "Man(b1-4)GlcNAc" = FALSE, "GlcNAc" = FALSE)
  expect_equal(result, expected)
})


test_that("have_motif: glycan graphs", {
  motif <- glyparse::parse_iupac_condensed("Gal")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed)
  expect_equal(have_motif(glycans, motif), c(TRUE, FALSE, FALSE))
})


test_that("have_motif: motif names", {
  glycans <- c(
    G1 = "Gal(b1-3)Gal(b1-3)GalNAc(a1-",
    G2 = "Gal(b1-3)GalNAc(a1-?)GlcNAc"
  )
  motif <- "O-Glycan core 1"
  result <- have_motif(glycans, motif)
  expected <- c(G1 = TRUE, G2 = FALSE)
  expect_equal(result, expected)
})


test_that("have_motif: motif names and custom alignment", {
  glycans <- c(
    G1 = "Gal(b1-3)Gal(b1-3)GalNAc(a1-",
    G2 = "Gal(b1-3)GalNAc(a1-?)GlcNAc"
  )
  motif <- "O-Glycan core 1"
  expect_snapshot(have_motif(glycans, motif, alignment = "whole"))
})


test_that("have_motif: some bad glycans", {
  motif <- "Gal"
  glycans <- c("Gal", "bad1", "bad2")
  expect_snapshot(have_motif(glycans, motif), error = TRUE)
})


test_that("have_motif: empty glycans", {
  motif <- "Gal"
  expected <- logical(0)
  names(expected) <- character(0)
  expect_equal(have_motif(character(0), motif), expected)
})


test_that("have_motif: glycan names are kept for IUPAC-condensed", {
  motif <- "Gal"
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  result <- have_motif(glycans, motif)
  expect_equal(names(result), names(glycans))
})


test_that("have_motif: glycan names are kept for glycan graphs", {
  motif <- glyparse::parse_iupac_condensed("Gal")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed)
  result <- have_motif(glycans, motif)
  expect_equal(names(result), names(glycans))
})


test_that("have_motif: glycans have different mono types", {
  motif <- "Gal"
  glycans <- c("Gal", "Hex")
  expect_snapshot(have_motif(glycans, motif), error = TRUE)
})


test_that("have_motif: glycans have different mono types with motif", {
  motif <- "Hex"
  glycans <- c(G1 = "Gal", G2 = "GlcNAc")
  expect_equal(have_motif(glycans, motif), c(G1 = TRUE, G2 = FALSE))
})


test_that("have_motif: glycans have obscurer mono types than motif", {
  motif <- "Gal"
  glycans <- c("Hex", "HexNAc")
  expect_snapshot(have_motif(glycans, motif), error = TRUE)
})


test_that("have_motif: glycans and motif have different graph modes", {
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal", mode = "ne")
  glycans <- c(G1 = "Gal(b1-3)Gal", G2 = "Gal(b1-4)Gal")
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed, mode = "dn")
  expect_equal(have_motif(glycans, motif), c(G1 = TRUE, G2 = FALSE))
})


# ----- have_motifs() -----
test_that("have_motifs: IUPAC-condensed (no names)", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    TRUE, TRUE, FALSE,
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- motifs
  rownames(expected) <- glycans
  expect_equal(result, expected)
})


test_that("have_motifs: IUPAC-condensed (with names)", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    TRUE, TRUE, FALSE,
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- c("M1", "M2", "M3")
  rownames(expected) <- c("G1", "G2", "G3")
  expect_equal(result, expected)
})


test_that("have_motifs: glycan graphs (no names)", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed)

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    TRUE, TRUE, FALSE,
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE)
  expect_equal(result, expected)
})


test_that("have_motifs: glycan graphs (with names)", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed)

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    TRUE, TRUE, FALSE,
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- c("M1", "M2", "M3")
  rownames(expected) <- c("G1", "G2", "G3")
  expect_equal(result, expected)
})


test_that("have_motifs: motif names (no names)", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  glycans <- c("Gal(b1-3)GalNAc(a1-3)GlcNAc(b1-", "Gal(b1-3)GalNAc(a1-")

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    FALSE, FALSE,
    TRUE, FALSE
  ), nrow = 2, byrow = TRUE)
  colnames(expected) <- motifs
  rownames(expected) <- glycans
  expect_equal(result, expected)
})


test_that("have_motifs: motif names (with names)", {
  motifs <- c(M1 = "O-Glycan core 1", M2 = "O-Glycan core 2")
  glycans <- c(G1 = "Gal(b1-3)GalNAc(a1-3)GlcNAc(b1-", G2 = "Gal(b1-3)GalNAc(a1-")

  result <- have_motifs(glycans, motifs)

  expected <- matrix(c(
    FALSE, FALSE,
    TRUE, FALSE
  ), nrow = 2, byrow = TRUE)
  colnames(expected) <- c("M1", "M2")
  rownames(expected) <- c("G1", "G2")
  expect_equal(result, expected)
})


test_that("have_motifs: custom alignments", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "GlcNAc(b1-4)GalNAc", G3 = "GlcNAc")

  result <- have_motifs(glycans, motifs, alignments = c("whole", "core", "terminal"))

  # If `alignments` are not specified, the result will be:
  # TRUE  TRUE FALSE
  # FALSE TRUE TRUE
  # FALSE TRUE FALSE

  expected <- matrix(c(
    FALSE, TRUE, FALSE,
    FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- c("M1", "M2", "M3")
  rownames(expected) <- c("G1", "G2", "G3")
  expect_equal(result, expected)
})


test_that("have_motifs: alignments provided for known motif names", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  glycans <- c("Gal(b1-3)GalNAc(a1-3)GlcNAc(b1-", "GalNAc(a1-3)Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")

  # should be all FALSE
  expect_snapshot(have_motifs(glycans, motifs, alignments = "whole"))
})


test_that("have_motifs: alignments have wrong length", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  alignments <- c("substructure", "core")
  expect_snapshot(have_motifs(glycans, motifs, alignments = alignments), error = TRUE)
})


test_that("have_motifs: bad motif names", {
  motifs <- c("GM1", "bad1", "bad2")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: bad IUPAC in motifs", {
  motifs <- c("Gal", "bad1", "bad2")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: bad IUPAC in glycans", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycans <- c("Gal(b1-3)GlcNAc", "bad1", "bad2")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: simplify", {
  motifs <- c("Gal", "GlcNAc", "Xyl")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  result <- have_motifs(glycans, motifs, simplify = TRUE)
  expected <- matrix(c(
    TRUE, TRUE,
    FALSE, TRUE,
    FALSE, TRUE
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- c("Gal", "GlcNAc")
  rownames(expected) <- glycans
  expect_equal(result, expected)
})


test_that("have_motifs: glycans have different mono types", {
  motifs <- c("Gal", "Gal")
  glycans <- c("Gal", "Hex")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: motifs have different mono types", {
  motifs <- c("Gal", "Hex")
  glycans <- c("Gal", "Gal")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: glycans have different mono types with motifs", {
  motifs <- c(M1 = "HexNAc", M2 = "Hex")
  glycans <- c(G1 = "Gal", G2 = "GlcNAc")
  expected <- matrix(c(
    FALSE, TRUE,
    TRUE, FALSE
  ), nrow = 2, byrow = TRUE)
  colnames(expected) <- c("M1", "M2")
  rownames(expected) <- c("G1", "G2")
  expect_equal(have_motifs(glycans, motifs), expected)
})


test_that("have_motifs: glycans have obscurer mono types than motifs", {
  motifs <- c("Gal", "GlcNAc")
  glycans <- c("Hex", "HexNAc")
  expect_snapshot(have_motifs(glycans, motifs), error = TRUE)
})


test_that("have_motifs: glycans and motifs have different graph modes", {
  motifs <- c(M1 = "Gal(b1-3)GlcNAc", M2 = "Gal(b1-4)GlcNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed, mode = "ne")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Gal(b1-4)GlcNAc")
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed, mode = "dn")
  expected <- matrix(c(
    TRUE, FALSE,
    FALSE, TRUE
  ), nrow = 2, byrow = TRUE)
  colnames(expected) <- c("M1", "M2")
  rownames(expected) <- c("G1", "G2")
  expect_equal(have_motifs(glycans, motifs), expected)
})


test_that("have_motifs: missing motifs", {
  glycans <- c(G1 = "Gal(b1-3)GalNAc", G2 = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc")
  expect_snapshot(have_motifs(glycans, simplify = TRUE))
})


test_that("have_motifs: simplify to dim 1", {
  # This is to test that, if only one motif is present,
  # and `simplify = TRUE`, the return value is still a matrix.
  # Note that subseting a matrix without `drop = FALSE` will reduce the result
  # to the lowest dimension. (e.g. matrix -> vector)
  # This is a bug in `have_motifs()` before.
  glycans <- c(G1 = "Gal", G2 = "GlcNAc")
  motifs <- c(M1 = "Gal", M2 = "Glc")
  expected <- matrix(c(TRUE, FALSE), ncol = 1)
  colnames(expected) <- "M1"
  rownames(expected) <- c("G1", "G2")
  expect_equal(have_motifs(glycans, motifs, simplify = TRUE), expected)
})


# ----- counts_motifs() -----
# As "counts_motif" functions use the same argument processing code
# as "have_motif" functions, we only need to test the main functionality.

test_that("counts_motifs: basic test", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycan <- c("Gal(b1-3)GlcNAc")
  result <- counts_motifs(glycan, motifs)
  expected <- c(M1 = 1, M2 = 1, M3 = 0)
  expect_equal(result, expected)
})


test_that("count_motif: basic test", {
  motif <- "Gal"
  glycans <- c(G1 = "Gal", G2 = "GlcNAc", G3 = "Gal(b1-3)Gal")
  result <- count_motif(glycans, motif)
  expected <- c(G1 = 1, G2 = 0, G3 = 2)
  expect_equal(result, expected)
})


test_that("count_motifs: basic test", {
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Gal(b1-3)Gal", G3 = "GlcNAc")
  result <- count_motifs(glycans, motifs)
  expected <- matrix(c(
    1, 1, 0,
    2, 0, 0,
    0, 1, 0
  ), nrow = 3, byrow = TRUE)
  colnames(expected) <- c("M1", "M2", "M3")
  rownames(expected) <- c("G1", "G2", "G3")
  expect_equal(result, expected)
})
