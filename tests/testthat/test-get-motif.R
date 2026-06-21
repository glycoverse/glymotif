test_that("db_motifs() creates a db_motifs_spec object", {
  spec <- db_motifs()
  expect_s3_class(spec, "db_motifs_spec")
  expect_snapshot(db_motifs())
})


test_that("db_motif_info() returns database motif metadata", {
  info <- db_motif_info()

  expect_s3_class(info, "tbl_df")
  expect_named(
    info,
    c(
      "source",
      "source_id",
      "accession",
      "name",
      "alignment",
      "glycan_structure"
    )
  )
  expect_true(all(info$source == "GlyGen"))
  expect_true(all(info$source_id == "GGM"))
  expect_s3_class(info$glycan_structure, "glyrepr_structure")
  expect_equal(length(info$glycan_structure), nrow(info))
})
