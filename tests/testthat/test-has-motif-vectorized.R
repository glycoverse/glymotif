# ----- has_motifs() -----
test_that("`has_motifs()` works for IUPAC-condensed", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(Gal = TRUE, GlcNAc = TRUE, GalNAc = FALSE))
})


test_that("`has_motifs()` works for motif names", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  # `glycan` is actually o-glycan core 2, with an additional
  # GlcNAc at the reducing end.
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  # Note that these two motifs both have "core" alignment,
  # so when `alignments` is not provided, it will use the default alignment.
  expect_equal(result, c("O-Glycan core 1" = FALSE, "O-Glycan core 2" = FALSE))
})


test_that("`has_motifs()` works for motif names with custom alignments", {
  motifs <- c("O-Glycan core 1", "O-Glycan core 2")
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-3)GlcNAc"
  # Output should be c("O-Glycan core 1" = TRUE, "O-Glycan core 2" = TRUE)
  expect_snapshot(has_motifs(glycan, motifs, alignments = "substructure"))
})


test_that("`has_motifs()` works for glycan graphs", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(TRUE, TRUE, FALSE))
})


test_that("some bad IUPAC", {
  motifs <- c("Gal", "bad1", "bad2")
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("some bad motif names", {
  motifs <- c("GM1", "bad1", "bad2")
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("wrong types", {
  motifs <- 1:3
  glycan <- "Gal(b1-3)GlcNAc"
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("empty motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  expected <- logical(0)
  names(expected) <- character(0)
  expect_equal(has_motifs(glycan, character(0)), expected)
})


test_that("motif names are kept", {
  glycan <- "Gal(b1-3)GlcNAc"
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(names(result), names(motifs))
})


test_that("missing motifs uses default motifs", {
  glycan <- "Gal(b1-3)GlcNAc"
  result <- has_motifs(glycan)
  expect_snapshot(result[1:10])
})


# ----- have_motif() -----
test_that("`have_motif()` works for IUPAC-condensed", {
  motif <- "Gal"
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  result <- have_motif(glycans, motif)
  expected <- c("Gal(b1-3)GlcNAc" = TRUE, "Man(b1-4)GlcNAc" = FALSE, "GlcNAc" = FALSE)
  expect_equal(result, expected)
})


test_that("`have_motif()` works for glycan graphs", {
  motif <- glyparse::parse_iupac_condensed("Gal")
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  glycans <- purrr::map(glycans, glyparse::parse_iupac_condensed)
  expect_equal(have_motif(glycans, motif), c(TRUE, FALSE, FALSE))
})


test_that("`have_motif()` works for motif name", {
  glycans <- c(
    G1 = "Gal(b1-3)Gal(b1-3)GalNAc(a1-",
    G2 = "Gal(b1-3)GalNAc(a1-?)GlcNAc"
  )
  motif <- "O-Glycan core 1"
  result <- have_motif(glycans, motif)
  expected <- c(G1 = TRUE, G2 = FALSE)
  expect_equal(result, expected)
})


test_that("`have_motif()` with motif name and custom alignment", {
  glycans <- c(
    G1 = "Gal(b1-3)Gal(b1-3)GalNAc(a1-",
    G2 = "Gal(b1-3)GalNAc(a1-?)GlcNAc"
  )
  motif <- "O-Glycan core 1"
  expect_snapshot(have_motif(glycans, motif, alignment = "whole"))
})


test_that("some bad glycans", {
  motif <- "Gal"
  glycans <- c("Gal", "bad1", "bad2")
  expect_snapshot(have_motif(glycans, motif), error = TRUE)
})


test_that("empty glycans", {
  motif <- "Gal"
  expected <- logical(0)
  names(expected) <- character(0)
  expect_equal(have_motif(character(0), motif), expected)
})


test_that("glycan names are kept", {
  motif <- "Gal"
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  result <- have_motif(glycans, motif)
  expect_equal(names(result), names(glycans))
})
