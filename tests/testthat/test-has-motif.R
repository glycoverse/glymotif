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


test_that("concrete glycan and generic motif", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyrepr::convert_glycan_mono_type(motif, "generic")
  expect_true(has_motif(glycan, motif))
})


test_that("concrete glycan and simple motif", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyrepr::convert_glycan_mono_type(motif, "simple")
  expect_true(has_motif(glycan, motif))
})


test_that("generic glycan and simple motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "generic")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyrepr::convert_glycan_mono_type(motif, "simple")
  expect_true(has_motif(glycan, motif))
})


test_that("simple glycan and concrete motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "simple")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_snapshot(has_motif(glycan, motif), error = TRUE)
})


test_that("simple glycan and generic motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "simple")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyrepr::convert_glycan_mono_type(motif, "generic")
  expect_snapshot(has_motif(glycan, motif), error = TRUE)
})


test_that("generic glycan and concrete motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "generic")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_snapshot(has_motif(glycan, motif), error = TRUE)
})


test_that("obscure linkages in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_false(has_motif(glycan, motif))
})


test_that("obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc")
  expect_true(has_motif(glycan, motif))
})


test_that("many obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)Gal(b?-3)GalNAc")
  expect_true(has_motif(glycan, motif))
})


test_that("obscure linkage in motif but false", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b2-?)GalNAc")
  expect_false(has_motif(glycan, motif))
})


test_that("obscure linkage in both motif and glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc")
  expect_true(has_motif(glycan, motif))
})


test_that("obscure linkage in DN motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc", mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc", mode = "dn")
  expect_true(has_motif(glycan, motif))
})


test_that("obscure linkage in DN motif but false", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc", mode = "dn")
  motif <- glyparse::parse_iupac_condensed("Gal(b2-?)GalNAc", mode = "dn")
  expect_false(has_motif(glycan, motif))
})


test_that("whole alignment simple positive case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyrepr::o_glycan_core_2()
  expect_true(has_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment simple negative case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyrepr::o_glycan_core_1()
  expect_false(has_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyrepr::n_glycan_core()
  expect_true(has_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_false(has_motif(glycan, motif, alignment = "whole"))
})


test_that("core alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal")
  expect_true(has_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal")
  expect_false(has_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_true(has_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_false(has_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-4)GlcNAc")
  expect_true(has_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_false(has_motif(glycan, motif, alignment = "core"))
})


test_that("terminal alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal")
  expect_true(has_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal")
  expect_false(has_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_true(has_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_false(has_motif(glycan, motif, alignment = "terminal"))
})


test_that("termimal alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
  expect_true(has_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)[Fuc(a1-3)]GlcNAc")
  expect_false(has_motif(glycan, motif, alignment = "terminal"))
})
