test_that("get available motifs", {
  expect_snapshot(available_motifs())
})


test_that("check if motif is known", {
  expect_true(is_known_motif("N-Glycan core basic"))
  expect_false(is_known_motif("bad motif"))
})


test_that("get known motif", {
  expect_error(get_motif("bad motif"), 'Unknown motif: "bad motif".')
})


test_that("getting motif works", {
  result <- get_motif("N-Glycan core basic")
  expect_length(result, 3)
  expect_s3_class(result$graph, "glycan_graph")
  expect_identical(result$aglycon, "Asn")
  expect_identical(result$alignment, "core")
})
