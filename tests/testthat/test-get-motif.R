test_that("get available motifs", {
  expect_snapshot(available_motifs())
})


test_that("check if motif is known", {
  expect_true(is_known_motif("N-Glycan core basic"))
  expect_false(is_known_motif("bad motif"))
})


test_that("get known motif", {
  expect_error(get_motif_graph("bad motif"), 'Unknown motif: "bad motif".')
})


test_that("getting motif graph works", {
  result <- get_motif_graph("N-Glycan core basic")
  expect_snapshot(print(result))
})


test_that("getting motif alignment works", {
  result <- get_motif_alignment("N-Glycan core basic")
  expect_identical(result, "core")
})


test_that("getting motif aglycon works", {
  result <- get_motif_aglycon("N-Glycan core basic")
  expect_identical(result, "Asn")
})
