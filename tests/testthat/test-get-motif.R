test_that("get all motifs", {
  expect_snapshot(db_motifs())
})


test_that("check if motif is known", {
  expect_true(is_known_motif("N-Glycan core basic"))
  expect_false(is_known_motif("bad motif"))
})


test_that("get known motif", {
  expect_error(get_motif_structure("bad motif"), 'Unknown motif: "bad motif".')
})


test_that("getting motif structure works", {
  result <- get_motif_structure("N-Glycan core basic")
  expect_snapshot(result)  # A `glyrepr_structure`
})


test_that("getting many motif graphs", {
  # The order in the result should be the same as the order in the input.
  # There was a bug that the order was not preserved.
  result1 <- get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
  result1 <- as.character(result1)
  expected1 <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  expect_identical(result1, expected1)

  result2 <- get_motif_structure(c("O-Glycan core 2", "O-Glycan core 1"))
  result2 <- as.character(result2)
  expected2 <- c("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-", "Gal(b1-3)GalNAc(a1-")
  expect_identical(result2, expected2)
})


test_that("getting motif alignment works", {
  result <- get_motif_alignment("N-Glycan core basic")
  expect_identical(result, "core")
})


test_that("getting many motif alignments", {
  result1 <- get_motif_alignment(c("O-Glycan core 1", "Lewis x"))
  result1 <- as.character(result1)
  expected1 <- c("core", "substructure")
  expect_identical(result1, expected1)

  result2 <- get_motif_alignment(c("Lewis x", "O-Glycan core 1"))
  result2 <- as.character(result2)
  expected2 <- c("substructure", "core")
  expect_identical(result2, expected2)
})


test_that("getting motif aglycon works", {
  result <- get_motif_aglycon("N-Glycan core basic")
  expect_identical(result, "Asn")
})


test_that("getting many motif aglycons", {
  result1 <- get_motif_aglycon(c("O-Glycan core 1", "O-Glycan core 2"))
  result1 <- as.character(result1)
  expected1 <- c("Ser/Thr", "Ser/Thr")
  expect_identical(result1, expected1)

  result2 <- get_motif_aglycon(c("O-Glycan core 2", "O-Glycan core 1"))
  result2 <- as.character(result2)
  expected2 <- c("Ser/Thr", "Ser/Thr")
  expect_identical(result2, expected2)
})
