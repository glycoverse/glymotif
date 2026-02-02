test_that("dynamic_motifs() creates a dynamic_motifs_spec object", {
  spec <- dynamic_motifs()
  expect_s3_class(spec, "dynamic_motifs_spec")
  expect_equal(spec$max_size, 3)
})

test_that("dynamic_motifs() accepts custom max_size", {
  spec <- dynamic_motifs(max_size = 5)
  expect_equal(spec$max_size, 5)
})

test_that("dynamic_motifs() validates max_size", {
  expect_error(dynamic_motifs(max_size = "invalid"), "Must be")
  expect_error(dynamic_motifs(max_size = 0), ">=")
})

test_that("branch_motifs() creates a branch_motifs_spec object", {
  spec <- branch_motifs()
  expect_s3_class(spec, "branch_motifs_spec")
})
