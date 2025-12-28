# ========== Input Data Types ==========
test_that("wrong glycan types", {
  glycan <- igraph::make_empty_graph()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc(b1-")
  expect_error(have_motif(glycan, motif), "`glycans` must be a 'glyrepr_structure' object")
})


test_that("wrong motif types", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- igraph::make_empty_graph()
  expect_error(have_motif(glycan, motif), "The `motif` argument must be a scalar vector")
})


test_that("IUPAC-condensed used as input", {
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(b1-"
  motif <- "Gal(b1-3)GalNAc(b1-"
  expect_true(have_motif(glycan, motif))
})


test_that("auto_parse supports multiple glycan structure formats", {
  # Test IUPAC-condensed format
  glycan_iupac <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(b1-"
  motif_iupac <- "Gal(b1-3)GalNAc(b1-"
  expect_true(have_motif(glycan_iupac, motif_iupac))

  # Test IUPAC-short format
  glycan_short <- "Galb3(GlcNAcb6)GalNAca-"
  motif_short <- "Galb3GalNAca-"
  expect_true(have_motif(glycan_short, motif_short))
})


test_that("motif name used as input", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- "O-Glycan core 1"
  expect_true(have_motif(glycan, motif))
})


test_that("unkown motif name used as input", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- "unknown motif name"
  expect_error(have_motif(glycan, motif), "Some motifs are neither valid glycan structures nor known motif names")
})


test_that("bad glycan structure", {
  glycan <- "bad structure"
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(b1-")
  expect_error(have_motif(glycan, motif), "Some glycans could not be parsed as valid glycan structures")
})


test_that("warning when user-provided alignment is different from database", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-3)GalNAc(a1-")
  motif <- "O-Glycan core 1"
  expect_warning(
    result <- have_motif(glycan, motif, alignment = "terminal"),
    "The provided alignment type.*is different from.*the motif's alignment type"
  )
  expect_false(result)
})


# ========== Monosaccharide Types ==========
test_that("concrete glycan and generic motif", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(?1-")
  motif <- glyrepr::convert_to_generic(motif)
  expect_true(have_motif(glycan, motif))
})


test_that("generic glycan and concrete motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "generic")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(?1-")
  expect_false(have_motif(glycan, motif))
})


# ========== Basic Topologies ==========
test_that("simple positive case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(?1-")
  expect_true(have_motif(glycan, motif))
})


test_that("simple negative case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc(?1-")
  expect_false(have_motif(glycan, motif))
})


test_that("complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif))
})


test_that("complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-")
  expect_false(have_motif(glycan, motif))
})


test_that("motif larger than glycan returns FALSE", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  expect_false(have_motif(glycan, motif))
})


test_that("multiple instances of motif in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  expect_true(have_motif(glycan, motif))
})


# ========== Linkages =========
test_that("ignoring linkages", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc(?1-")
  expect_true(have_motif(glycan, motif, ignore_linkages = TRUE))
})


patrick::with_parameters_test_that("obscure linkages in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  expect_false(have_motif(glycan, motif))
})


patrick::with_parameters_test_that("obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  motif_iupac <- stringr::str_glue("Gal({linkage})GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed(motif_iupac)
  expect_true(have_motif(glycan, motif))
}, linkage = c("b1-?", "b?-3", "b?-?", "??-?", "b1-3/6", "b1-6/3"))


test_that("many obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)Gal(b?-3)GalNAc(a1-")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkage in motif but false", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b2-?)GalNAc(a1-")
  expect_false(have_motif(glycan, motif))
})


test_that("obscure linkage in both motif and glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc(a1-")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkages on the same edge", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc(a1-")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkages on the same edge, motif obscurer", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(?1-?)GalNAc(a1-")
  expect_true(have_motif(glycan, motif))
})


# ========== Alignment ==========
test_that("whole alignment simple positive case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyrepr::o_glycan_core_2()
  expect_true(have_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment simple negative case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyrepr::o_glycan_core_1()
  expect_false(have_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyrepr::n_glycan_core()
  expect_true(have_motif(glycan, motif, alignment = "whole"))
})


test_that("whole alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(a1-")
  expect_false(have_motif(glycan, motif, alignment = "whole"))
})


test_that("core alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal(a1-")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(a1-")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(a1-")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc(?1-")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc(?1-")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("terminal alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(?1-")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(?1-")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(?1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal(?1-")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("termimal alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)[Fuc(a1-3)]GlcNAc(?1-")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("exact alignment simple positive case", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-6)GalNAc(a1-")
  motif <- glyrepr::as_glycan_structure("Gal(b1-")
  expect_true(have_motif(glycan, motif, alignment = "exact"))
})


test_that("exact alignment simple negative case 1", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc(?1-")
  motif <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_false(have_motif(glycan, motif, alignment = "exact"))
})


test_that("exact alignment simple negative case 2", {
  glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  motif <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_false(have_motif(glycan, motif, alignment = "exact"))
})


test_that("exact alignment complex positive case", {
  glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc(?1-")
  motif <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_true(have_motif(glycan, motif, alignment = "exact"))
})


test_that("custom alignment is same as database", {
  expect_false(have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc(?1-", "O-Glycan core 1", alignment = "core"))
  expect_true(have_motif("Gal(b1-3)GalNAc(a1-", "O-Glycan core 1", alignment = "core"))
})


test_that("custom alignment is different from database", {
  expect_warning(
    result <- have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc(?1-", "O-Glycan core 1", alignment = "substructure"),
    "The provided alignment type.*is different from.*the motif's alignment type"
  )
  expect_true(result)
})


# ========== Substituents ==========
test_that("substituents are considered", {
  glycan1 <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  glycan2 <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")

  expect_true(have_motif(glycan1, glycan1))
  expect_false(have_motif(glycan1, glycan2))
  expect_false(have_motif(glycan2, glycan1))
  expect_true(have_motif(glycan2, glycan2))
})


test_that("obscure substituent linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac?Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure substituent linkages in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  expect_false(have_motif(glycan, motif))
})


test_that("multiple substituents in motif", {
  glycan1 <- glyparse::parse_iupac_condensed("Glc3Me6S(?1-")
  glycan2 <- glyparse::parse_iupac_condensed("Glc?Me6S(?1-")
  glycan3 <- glyparse::parse_iupac_condensed("Glc3Me(?1-")

  motif1 <- glyparse::parse_iupac_condensed("Glc3Me6S(?1-")
  motif2 <- glyparse::parse_iupac_condensed("Glc?Me6S(?1-")
  motif3 <- glyparse::parse_iupac_condensed("Glc3Me?S(?1-")
  motif4 <- glyparse::parse_iupac_condensed("Glc3Me(?1-")
  motif5 <- glyparse::parse_iupac_condensed("Glc(?1-")

  expect_true(have_motif(glycan1, motif1))
  expect_true(have_motif(glycan1, motif2))
  expect_true(have_motif(glycan1, motif3))
  expect_false(have_motif(glycan1, motif4))
  expect_false(have_motif(glycan1, motif5))
  expect_false(have_motif(glycan2, motif1))
  expect_true(have_motif(glycan2, motif2))
  expect_false(have_motif(glycan3, motif1))
})


test_that("substituents are not mandatory if strict_sub is FALSE", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")

  expect_true(have_motif(glycan, motif, strict_sub = FALSE))
  expect_false(have_motif(glycan, motif, strict_sub = TRUE))
})


# ========== Anomers ==========
test_that("anomer right for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-")
  expect_true(have_motif(glycan, motif))
})


patrick::with_parameters_test_that(
  "obscure anomer right for inside motifs",
  {
    glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
    motif <- glyparse::parse_iupac_condensed(motif_iupac)
    expect_true(have_motif(glycan, motif))
  },
  motif_iupac = c("Gal(?1-", "Gal(b?-", "Gal(??-")
)


test_that("anomer wrong for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Gal(a1-")
  expect_false(have_motif(glycan, motif))
})


test_that("obscure anomer wrong for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif <- glyparse::parse_iupac_condensed("Gal(a?-")
  expect_false(have_motif(glycan, motif))
})


test_that("anomer right for core motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-")
  expect_true(have_motif(glycan, motif))
})


test_that("anomer wrong for core motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(a1-")
  expect_false(have_motif(glycan, motif))
})


# ========== have_motifs ==========
test_that("have_motifs works with multiple motifs", {
  glycan1 <- as.character(glyrepr::o_glycan_core_2())
  glycan2 <- "Gal(b1-?)[GlcNAc(b1-6)]GalNAc(?1-"
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("core2", "test_glycan")
  
  motifs <- c("Gal(b1-3)GalNAc(?1-", "Gal(b1-4)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), c("Gal(b1-3)GalNAc(?1-", "Gal(b1-4)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-"))
  expect_equal(rownames(result), c("core2", "test_glycan"))
  
  # Check specific results
  expect_true(result[1, "Gal(b1-3)GalNAc(?1-"])
  expect_false(result[1, "Gal(b1-4)GalNAc(?1-"])
  expect_true(result[1, "GlcNAc(b1-6)GalNAc(?1-"])
  
  expect_false(result[2, "Gal(b1-3)GalNAc(?1-"])
  expect_false(result[2, "Gal(b1-4)GalNAc(?1-"])
  expect_true(result[2, "GlcNAc(b1-6)GalNAc(?1-"])
})


test_that("have_motifs works without glycan names", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc(?1-")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-"))
  
  # Check that row names contain IUPAC strings
  expect_true(grepl("Gal\\(b1-3\\)", rownames(result)[1]))
  expect_true(grepl("Gal\\(b1-\\?\\)", rownames(result)[2]))
})


test_that("have_motifs works with different alignments", {
  glycan <- glyparse::parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(?1-")
  glycans <- c(glycan, glycan)
  
  motifs <- c("Gal(a1-3)Gal(a1-4)Gal(?1-", "Gal(a1-4)Gal(a1-6)Gal(?1-")
  alignments <- c("core", "terminal")
  
  result <- have_motifs(glycans, motifs, alignments = alignments)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("have_motifs handles empty motifs", {
  glycan <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan)
  motifs <- character(0)
  
  expect_error(have_motifs(glycans, motifs), "`motifs` cannot be empty")
})


test_that("have_motifs handles invalid motifs argument", {
  glycan <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan)
  motifs <- 123

  expect_error(have_motifs(glycans, motifs), "`motifs` must be a 'glyrepr_structure' object")
})


test_that("have_motifs handles mismatched alignments length", {
  glycan <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan)
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-", "Gal(b1-4)GalNAc(?1-")
  alignments <- c("substructure", "core")  # length 2, motifs length 3
  
  expect_error(have_motifs(glycans, motifs, alignments = alignments), 
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})


test_that("have_motifs works with single alignment for all motifs", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc(?1-")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-")
  
  result <- have_motifs(glycans, motifs, alignments = "substructure")
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("have_motifs works with ignore_linkages", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1)
  
  motifs <- c("Gal(b1-3)GalNAc(?1-", "Gal(b1-4)GalNAc(?1-")
  
  result_normal <- have_motifs(glycans, motifs, ignore_linkages = FALSE)
  result_ignore <- have_motifs(glycans, motifs, ignore_linkages = TRUE)
  
  expect_false(result_normal[1, "Gal(b1-4)GalNAc(?1-"])
  expect_true(result_ignore[1, "Gal(b1-4)GalNAc(?1-"])
})


# ========== Monosaccharide Type Matching ==========
test_that("have_motif: concrete glycan matches generic motif", {
  # Man (concrete) should match Hex (generic) after conversion
  expect_true(have_motif("Man(?1-", "Hex(?1-"))
  expect_true(have_motif("Gal(?1-", "Hex(?1-"))
  expect_true(have_motif("Glc(?1-", "Hex(?1-"))
})

test_that("have_motif: generic glycan does not match concrete motif", {
  # Hex (generic) should not match Man (concrete) - names don't match
  expect_false(have_motif("Hex(?1-", "Man(?1-"))
  expect_false(have_motif("Hex(?1-", "Gal(?1-"))
  expect_false(have_motif("Hex(?1-", "Glc(?1-"))
})

test_that("have_motif: same type matches", {
  # Concrete vs concrete - exact match
  expect_true(have_motif("Man(?1-", "Man(?1-"))
  expect_true(have_motif("Gal(?1-", "Gal(?1-"))
  
  # Generic vs generic - exact match
  expect_true(have_motif("Hex(?1-", "Hex(?1-"))
  expect_true(have_motif("HexNAc(?1-", "HexNAc(?1-"))
})

test_that("have_motif: complex structures with type conversion", {
  # Complex structure: concrete glycan should match when motif is generic
  concrete_glycan <- "Gal(b1-3)GalNAc(?1-"
  generic_motif <- "Hex(b1-3)HexNAc(?1-"
  expect_true(have_motif(concrete_glycan, generic_motif))
  
  # Reverse should not match
  generic_glycan <- "Hex(b1-3)HexNAc(?1-"
  concrete_motif <- "Gal(b1-3)GalNAc(?1-"
  expect_false(have_motif(generic_glycan, concrete_motif))
})

test_that("have_motifs: matrix with mixed types", {
  # Test the key example from the user's original question
  glycans <- c("Hex(?1-", "Man(?1-")
  motifs <- c("Hex(?1-", "Man(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  # Check structure
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_equal(rownames(result), c("Hex(?1-", "Man(?1-"))
  expect_equal(colnames(result), c("Hex(?1-", "Man(?1-"))
  
  # Check expected results
  expect_true(result["Hex(?1-", "Hex(?1-"])    # generic vs generic: match
  expect_false(result["Hex(?1-", "Man(?1-"])   # generic vs concrete: no match
  expect_true(result["Man(?1-", "Hex(?1-"])    # concrete vs generic: match (after conversion)
  expect_true(result["Man(?1-", "Man(?1-"])    # concrete vs concrete: match
})

test_that("have_motifs: only concrete motifs", {
  # Test case that used to throw error
  glycans <- c("Hex(?1-", "Man(?1-")
  motifs <- c("Man(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 1))
  expect_equal(rownames(result), c("Hex(?1-", "Man(?1-"))
  expect_equal(colnames(result), c("Man(?1-"))
  
  # Only Man glycan should match Man motif
  expect_false(result["Hex(?1-", "Man(?1-"])   # generic glycan vs concrete motif: no match
  expect_true(result["Man(?1-", "Man(?1-"])    # concrete glycan vs concrete motif: match
})

test_that("have_motifs: only generic motifs", {
  glycans <- c("Hex(?1-", "Man(?1-", "Gal(?1-")
  motifs <- c("Hex(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 1))
  expect_equal(rownames(result), c("Hex(?1-", "Man(?1-", "Gal(?1-"))
  expect_equal(colnames(result), c("Hex(?1-"))
  
  # All should match the generic motif
  expect_true(result["Hex(?1-", "Hex(?1-"])    # generic vs generic: match
  expect_true(result["Man(?1-", "Hex(?1-"])    # concrete vs generic: match (after conversion)
  expect_true(result["Gal(?1-", "Hex(?1-"])    # concrete vs generic: match (after conversion)
})

test_that("have_motifs: complex mixed type scenarios", {
  # More complex test with various type combinations
  glycans <- c("Hex(?1-", "Man(?1-", "Gal(?1-", "HexNAc(?1-", "GalNAc(?1-")
  motifs <- c("Hex(?1-", "Man(?1-", "HexNAc(?1-", "Fuc(?1-")
  
  result <- have_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(5, 4))
  
  # Test some key combinations
  expect_true(result["Hex(?1-", "Hex(?1-"])      # generic vs generic
  expect_false(result["Hex(?1-", "Man(?1-"])     # generic vs concrete
  expect_true(result["Man(?1-", "Hex(?1-"])      # concrete vs generic (Man converts to Hex)
  expect_true(result["Man(?1-", "Man(?1-"])      # concrete vs concrete (same)
  expect_false(result["Man(?1-", "HexNAc(?1-"])  # concrete vs concrete (different)
  expect_true(result["GalNAc(?1-", "HexNAc(?1-"]) # concrete vs generic (GalNAc converts to HexNAc)
})

test_that("have_motif: vectorized glycans with single motif", {
  # Test vectorized behavior with mixed types
  glycans <- c("Hex(?1-", "Man(?1-", "Gal(?1-")
  motif <- "Hex(?1-"  # generic motif
  
  result <- have_motif(glycans, motif)
  
  expect_length(result, 3)
  expect_true(all(result))  # All should match the generic motif
  
  # Test with concrete motif
  motif <- "Man(?1-"  # concrete motif
  result <- have_motif(glycans, motif)
  
  expect_length(result, 3)
  expect_false(result[1])  # Hex should not match Man
  expect_true(result[2])   # Man should match Man
  expect_false(result[3])  # Gal should not match Man
})

test_that("have_motif works for repeated glycans", {
  # This is to test the internal `fast_convert_to_generic` function.
  glycans <- c("Man(??-", "Man(??-", "Gal(??-", "Gal(??-", "GlcNAc(??-")
  motif <- "Hex(??-"

  result <- have_motif(glycans, motif)
  expect_equal(result, c(TRUE, TRUE, TRUE, TRUE, FALSE))
})
