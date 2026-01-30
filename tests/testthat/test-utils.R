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

test_that("has_duplicate_motifs detects duplicate character motifs", {
  motifs <- c("Hex(b1-", "Hex(b1-")
  expect_true(has_duplicate_motifs(motifs))
})

test_that("has_duplicate_motifs returns FALSE for unique character motifs", {
  motifs <- c("Gal(b1-", "GlcNAc(b1-")
  expect_false(has_duplicate_motifs(motifs))
})

test_that("has_duplicate_motifs detects duplicate glyrepr_structure motifs", {
  skip_if_not_installed("glyparse")
  motifs <- glyparse::parse_iupac_condensed(c("Hex(b1-", "Hex(b1-"))
  expect_true(has_duplicate_motifs(motifs))
})

test_that("has_duplicate_motifs returns FALSE for unique glyrepr_structure motifs", {
  skip_if_not_installed("glyparse")
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-", "GlcNAc(b1-"))
  expect_false(has_duplicate_motifs(motifs))
})

test_that("has_duplicate_motifs handles single motif", {
  expect_false(has_duplicate_motifs("Hex(b1-"))
})

test_that("has_duplicate_motifs handles empty motifs", {
  expect_false(has_duplicate_motifs(character(0)))
})
