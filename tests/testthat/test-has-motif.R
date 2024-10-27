test_that("simple positive case for NE glycan", {
  glycan <- glyrepr::o_glycan_core_2(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "ne")
  expect_true(has_motif(glycan, motif))
})


test_that("simple negative case for NE glycan", {
  glycan <- glyrepr::o_glycan_core_2(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc", mode = "ne")
  expect_false(has_motif(glycan, motif))
})


test_that("complex positive case for NE glycan", {
  glycan <- glyrepr::n_glycan_core(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc", mode = "ne")
  expect_true(has_motif(glycan, motif))
})


test_that("complex negative case for NE glycan", {
  glycan <- glyrepr::n_glycan_core(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc", mode = "ne")
  expect_false(has_motif(glycan, motif))
})


test_that("simple positive case for DN glycan", {
  glycan <- glyrepr::o_glycan_core_2(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "dn")
  expect_true(has_motif(glycan, motif))
})


test_that("simple negative case for DN glycan", {
  glycan <- glyrepr::o_glycan_core_2(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc", mode = "dn")
  expect_false(has_motif(glycan, motif))
})


test_that("complex positive case for DN glycan", {
  glycan <- glyrepr::n_glycan_core(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc", mode = "dn")
  expect_true(has_motif(glycan, motif))
})


test_that("complex negative case for DN glycan", {
  glycan <- glyrepr::n_glycan_core(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc", mode = "dn")
  expect_false(has_motif(glycan, motif))
})
