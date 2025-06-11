test_that("get available motifs", {
  expect_snapshot(available_motifs())
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
  result <- get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
  expect_snapshot(result)  # A list of `glyrepr_structure`
})


test_that("getting motif alignment works", {
  result <- get_motif_alignment("N-Glycan core basic")
  expect_identical(result, "core")
})


test_that("getting many motif alignments", {
  result <- get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
  expect_identical(result, c("O-Glycan core 1" = "core", "O-Glycan core 2" = "core"))
})


test_that("getting motif aglycon works", {
  result <- get_motif_aglycon("N-Glycan core basic")
  expect_identical(result, "Asn")
})


test_that("getting many motif aglycons", {
  result <- get_motif_aglycon(c("O-Glycan core 1", "O-Glycan core 2"))
  expect_identical(result, c("O-Glycan core 1" = "Ser/Thr", "O-Glycan core 2" = "Ser/Thr"))
})
