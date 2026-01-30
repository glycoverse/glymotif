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
