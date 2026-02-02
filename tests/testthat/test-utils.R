test_that("prepare_struc_names preserves glycan_structure names", {
  glycan1 <- glyrepr::o_glycan_core_1()
  glycan2 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("core1", "core2")

  result <- prepare_struc_names(glycans, glycans)
  expect_equal(result, c("core1", "core2"))
})

test_that("prepare_struc_names returns IUPAC when no names on glycan_structure", {
  glycan1 <- glyrepr::o_glycan_core_1()
  glycan2 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1, glycan2)

  result <- prepare_struc_names(glycans, glycans)
  expect_equal(result, as.character(glycans))
})

test_that("prepare_struc_names preserves character vector names", {
  chars <- c("a", "b")
  names(chars) <- c("name_a", "name_b")

  result <- prepare_struc_names(chars, chars)
  expect_equal(result, c("name_a", "name_b"))
})

test_that("prepare_struc_names returns values when character vector has no names", {
  chars <- c("a", "b")

  result <- prepare_struc_names(chars, chars)
  expect_equal(result, c("a", "b"))
})

test_that("prepare_motif_args handles dynamic_motifs_spec", {
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  spec <- dynamic_motifs(max_size = 2)

  result <- prepare_motif_args(glycans, spec, alignments = NULL, ignore_linkages = FALSE, match_degree = NULL, single_motif = FALSE, strict_sub = TRUE)

  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_equal(unique(result$alignments), "substructure")
  expect_null(result$match_degree)
})

test_that("prepare_motif_args handles branch_motifs_spec", {
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  spec <- branch_motifs()

  result <- prepare_motif_args(glycans, spec, alignments = NULL, ignore_linkages = FALSE, match_degree = NULL, single_motif = FALSE, strict_sub = TRUE)

  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_true(length(result$motifs) > 0)
  expect_type(result$match_degree, "list")
})

test_that("prepare_motif_args errors on alignments conflict with spec", {
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  spec <- dynamic_motifs()

  expect_error(
    prepare_motif_args(glycans, spec, alignments = "core", ignore_linkages = FALSE, match_degree = NULL, single_motif = FALSE, strict_sub = TRUE),
    "Cannot specify.*alignments"
  )
})
