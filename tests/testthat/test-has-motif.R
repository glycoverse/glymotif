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


test_that("motif larger than glycan returns FALSE", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc", mode = "ne")
  expect_false(has_motif(glycan, motif))
})


test_that("multiple instances of motif in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc", mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "ne")
  expect_true(has_motif(glycan, motif))
})



test_that("NE glycan ignore linkages", {
  glycan <- glyrepr::o_glycan_core_2(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc", mode = "ne")
  expect_true(has_motif(glycan, motif, ignore_linkages = TRUE))
})


test_that("DN glycan ignore linkages", {
  glycan <- glyrepr::o_glycan_core_2(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc", mode = "dn")
  expect_true(has_motif(glycan, motif, ignore_linkages = TRUE))
})


test_that("wrong glycan types", {
  glycan <- igraph::make_empty_graph()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc", mode = "ne")
  expect_error(has_motif(glycan, motif), "`glycan` and `motif` must be 'glycan_graph' objects.")
})


test_that("wrong motif types", {
  glycan <- glyrepr::o_glycan_core_2(mode = "ne")
  motif <- igraph::make_empty_graph()
  expect_error(has_motif(glycan, motif), "`glycan` and `motif` must be 'glycan_graph' objects.")
})


test_that("NULL inputs", {
  expect_error(has_motif(NULL, NULL), "`glycan` and `motif` must be 'glycan_graph' objects.")
})


test_that("invalid data types", {
  expect_error(has_motif(123, "motif"), "`glycan` and `motif` must be 'glycan_graph' objects.")
})


test_that("ND glycan and NE motif", {
  glycan <- glyrepr::o_glycan_core_2(mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "ne")
  expect_true(has_motif(glycan, motif))
})


test_that("NE glycan and DN motif", {
  glycan <- glyrepr::o_glycan_core_2(mode = "ne")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc", mode = "dn")
  expect_true(has_motif(glycan, motif))
})
