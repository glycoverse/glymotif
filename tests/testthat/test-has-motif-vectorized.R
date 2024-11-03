# ----- has_motifs() -----
test_that("`has_motifs()` works", {
  motifs <- c("Gal", "GlcNAc", "GalNAc")
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  expect_equal(has_motifs(glycan, motifs), c(TRUE, TRUE, FALSE))
})


test_that("some bad motifs", {
  motifs <- c("Gal", "bad1", "bad2")
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  expect_snapshot(has_motifs(glycan, motifs), error = TRUE)
})


test_that("empty motifs", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  expect_equal(has_motifs(glycan, character(0)), logical(0))
})


test_that("motif names are kept", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GlcNAc")
  motifs <- c(Gal = "Gal", GlcNAc = "GlcNAc", GalNAc = "GalNAc")
  result <- has_motifs(glycan, motifs)
  expect_equal(names(result), names(motifs))
})


# ----- have_motif() -----
test_that("`have_motif()` works", {
  motif <- "Gal"
  glycans <- c("Gal(b1-3)GlcNAc", "Man(b1-4)GlcNAc", "GlcNAc")
  expect_equal(have_motif(glycans, motif), c(TRUE, FALSE, FALSE))
})


test_that("some bad glycans", {
  motif <- "Gal"
  glycans <- c("Gal", "bad1", "bad2")
  expect_snapshot(have_motif(glycans, motif), error = TRUE)
})


test_that("empty glycans", {
  motif <- "Gal"
  expect_equal(have_motif(character(0), motif), logical(0))
})


test_that("glycan names are kept", {
  motif <- "Gal"
  glycans <- c(G1 = "Gal(b1-3)GlcNAc", G2 = "Man(b1-4)GlcNAc", G3 = "GlcNAc")
  result <- have_motif(glycans, motif)
  expect_equal(names(result), names(glycans))
})
