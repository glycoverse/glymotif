# ========== Input Data Types ==========
test_that("wrong glycan types", {
  glycan <- igraph::make_empty_graph()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc")
  expect_snapshot(have_motif(glycan, motif), error = TRUE)
})


test_that("wrong motif types", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- igraph::make_empty_graph()
  expect_snapshot(have_motif(glycan, motif), error = TRUE)
})


test_that("IUPAC-condensed used as input", {
  glycan <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc"
  motif <- "Gal(b1-3)GalNAc"
  expect_true(have_motif(glycan, motif))
})


test_that("motif name used as input", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- "O-Glycan core 1"
  expect_true(have_motif(glycan, motif))
})


test_that("unkown motif name used as input", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- "unknown motif name"
  expect_snapshot(have_motif(glycan, motif), error = TRUE)
})


test_that("bad glycan IUPAC", {
  glycan <- "bad IUPAC"
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_snapshot(have_motif(glycan, motif), error = TRUE)
})


test_that("warning when user-provided alignment is different from database", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-3)GalNAc(a1-")
  motif <- "O-Glycan core 1"
  expect_snapshot(have_motif(glycan, motif, alignment = "terminal"))
})


# ========== Monosaccharide Types ==========
test_that("concrete glycan and generic motif", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyrepr::convert_mono_type(motif, "generic")
  expect_true(have_motif(glycan, motif))
})


test_that("generic glycan and concrete motif", {
  glycan <- glyrepr::o_glycan_core_2(mono_type = "generic")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_snapshot(have_motif(glycan, motif), error = TRUE)
})


# ========== Basic Topologies ==========
test_that("simple positive case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("simple negative case", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc")
  expect_false(have_motif(glycan, motif))
})


test_that("complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc")
  expect_false(have_motif(glycan, motif))
})


test_that("motif larger than glycan returns FALSE", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  expect_false(have_motif(glycan, motif))
})


test_that("multiple instances of motif in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_true(have_motif(glycan, motif))
})


# ========== Linkages =========
test_that("ignoring linkages", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc")
  expect_true(have_motif(glycan, motif, ignore_linkages = TRUE))
})


patrick::with_parameters_test_that("obscure linkages in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  expect_false(have_motif(glycan, motif))
})


patrick::with_parameters_test_that("obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif_iupac <- stringr::str_glue("Gal({linkage})GalNAc")
  motif <- glyparse::parse_iupac_condensed(motif_iupac)
  expect_true(have_motif(glycan, motif))
}, linkage = c("b1-?", "b?-3", "b?-?", "??-?", "b1-3/6", "b1-6/3"))


test_that("many obscure linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)Gal(b?-3)GalNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkage in motif but false", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b2-?)GalNAc")
  expect_false(have_motif(glycan, motif))
})


test_that("obscure linkage in both motif and glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkages on the same edge", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-?)GalNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure linkages on the same edge, motif obscurer", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-?)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(?1-?)GalNAc")
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
  motif <- glyparse::parse_iupac_condensed("Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_false(have_motif(glycan, motif, alignment = "whole"))
})


test_that("core alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif, alignment = "core"))
})


test_that("core alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(b1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_false(have_motif(glycan, motif, alignment = "core"))
})


test_that("terminal alignment simple positive case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment simple negative case", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex positive case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-6)Man(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment complex negative case", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(b1-4)GlcNAc")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("termimal alignment more complex positive case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("terminal alignment more complex negative case", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)[Fuc(a1-3)]GlcNAc")
  expect_false(have_motif(glycan, motif, alignment = "terminal"))
})


test_that("custom alignment is same as database", {
  expect_false(have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc", "O-Glycan core 1", alignment = "core"))
  expect_true(have_motif("Gal(b1-3)GalNAc(a1-", "O-Glycan core 1", alignment = "core"))
})


test_that("custom alignment is different from database", {
  expect_snapshot(have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc", "O-Glycan core 1", alignment = "substructure"))
})


# ========== Substituents ==========
test_that("substituents are considered", {
  glycan1 <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc")
  glycan2 <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")

  expect_true(have_motif(glycan1, glycan1))
  expect_false(have_motif(glycan1, glycan2))
  expect_false(have_motif(glycan2, glycan1))
  expect_true(have_motif(glycan2, glycan2))
})


test_that("obscure substituent linkages in motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac?Ac(a2-3)Gal(b1-4)GlcNAc")
  expect_true(have_motif(glycan, motif))
})


test_that("obscure substituent linkages in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
  expect_false(have_motif(glycan, motif))
})


# ========== Anomers ==========
test_that("anomer right for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-")
  expect_true(have_motif(glycan, motif))
})


patrick::with_parameters_test_that(
  "obscure anomer right for inside motifs",
  {
    glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
    motif <- glyparse::parse_iupac_condensed(motif_iupac)
    expect_true(have_motif(glycan, motif))
  },
  motif_iupac = c("Gal(?1-", "Gal(b?-", "Gal(??-")
)


test_that("anomer wrong for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(a1-")
  expect_false(have_motif(glycan, motif))
})


test_that("obscure anomer wrong for inside motif", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")
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
  glycan1 <- glyrepr::o_glycan_core_2()
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("core2", "test_glycan")
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc", "GlcNAc(b1-6)GalNAc")
  
  result <- have_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
  expect_equal(names(result), c("glycan", "Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc", "GlcNAc(b1-6)GalNAc"))
  expect_equal(result$glycan, c("core2", "test_glycan"))
  
  # Check specific results
  expect_true(result[1, "Gal(b1-3)GalNAc"][[1]])
  expect_false(result[1, "Gal(b1-4)GalNAc"][[1]])
  expect_true(result[1, "GlcNAc(b1-6)GalNAc"][[1]])
  
  expect_false(result[2, "Gal(b1-3)GalNAc"][[1]])
  expect_false(result[2, "Gal(b1-4)GalNAc"][[1]])
  expect_true(result[2, "GlcNAc(b1-6)GalNAc"][[1]])
})


test_that("have_motifs works without glycan names", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc", "GlcNAc(b1-6)GalNAc")
  
  result <- have_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("glycan", "Gal(b1-3)GalNAc", "GlcNAc(b1-6)GalNAc"))
  
  # Check that glycan column contains IUPAC strings
  expect_true(grepl("Gal\\(b1-3\\)", result$glycan[1]))
  expect_true(grepl("Gal\\(b1-\\?\\)", result$glycan[2]))
})


test_that("have_motifs works with different alignments", {
  glycan <- glyparse::parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
  glycans <- c(glycan, glycan)
  names(glycans) <- c("glycan1", "glycan2")
  
  motifs <- c("Gal(a1-3)Gal(a1-4)Gal", "Gal(a1-4)Gal(a1-6)Gal")
  alignments <- c("core", "terminal")
  
  result <- have_motifs(glycans, motifs, alignments = alignments)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
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
  
  expect_error(have_motifs(glycans, motifs), "`motifs` must be a character vector")
})


test_that("have_motifs handles mismatched alignments length", {
  glycan <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan)
  motifs <- c("Gal(b1-3)GalNAc", "GlcNAc(b1-6)GalNAc", "Gal(b1-4)GalNAc")
  alignments <- c("substructure", "core")  # length 2, motifs length 3
  
  expect_error(have_motifs(glycans, motifs, alignments = alignments), 
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})


test_that("have_motifs works with single alignment for all motifs", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-?)[GlcNAc(b1-6)]GalNAc")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc", "GlcNAc(b1-6)GalNAc")
  
  result <- have_motifs(glycans, motifs, alignments = "substructure")
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})


test_that("have_motifs works with ignore_linkages", {
  glycan1 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1)
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc")
  
  result_normal <- have_motifs(glycans, motifs, ignore_linkages = FALSE)
  result_ignore <- have_motifs(glycans, motifs, ignore_linkages = TRUE)
  
  expect_false(result_normal[1, "Gal(b1-4)GalNAc"][[1]])
  expect_true(result_ignore[1, "Gal(b1-4)GalNAc"][[1]])
})
