test_that("count motifs in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-")
  expect_equal(count_motif(glycan, motif), 2L)
})


test_that("count motifs with branching", {
  glycan <- glyparse::parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc")
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count symmetrical motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal")
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count 0 motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc")
  expect_equal(count_motif(glycan, motif), 0L)
})
