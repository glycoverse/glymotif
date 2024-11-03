# ----- has_motifs() -----
test_that("`has_motifs()` works for IUPAC-condensed", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(Gal = TRUE, GlcNAc = TRUE, GalNAc = FALSE))
})


test_that("`has_motifs()` works for motif names", {
  motifs <- c("GM1", "GD1a", "GD1alpha")
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(GM1 = FALSE, GD1a = FALSE, GD1alpha = FALSE))
})


test_that("`has_motifs()` works for glycan graphs", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  motifs <- purrr::map(motifs, glyparse::parse_iupac_condensed)
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(result, c(TRUE, TRUE, FALSE))
})


test_that("some bad motifs", {
  motifs <- c("Gal", "bad1", "bad2")
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("empty motifs", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  expected <- logical(0)
  names(expected) <- character(0)
  expect_equal(has_motifs(glycan, character(0)), expected)
})


test_that("motif names are kept", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  motifs <- c(M1 = "Gal", M2 = "GlcNAc", M3 = "GalNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(names(result), names(motifs))
})


test_that("missing motifs uses default motifs", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
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
