test_that("db_motifs() creates a db_motifs_spec object", {
  spec <- db_motifs()
  expect_s3_class(spec, "db_motifs_spec")
  expect_identical(spec$source_id, "GGM")
  expect_snapshot(db_motifs())
})


test_that("db_motifs() validates source IDs", {
  spec <- db_motifs(source_id = c("GGM", "CCRC"))

  expect_s3_class(spec, "db_motifs_spec")
  expect_identical(spec$source_id, c("GGM", "CCRC"))
  expect_error(
    db_motifs(source_id = "UNKNOWN"),
    "Unknown motif source"
  )
})


test_that("db_motif_info() returns database motif metadata", {
  info <- db_motif_info()
  expected_sources <- c(
    GGM = "GlyGen Motifs",
    CCRC = "CCRC Motifs",
    GTC = "GlyTouCan Motifs",
    GE = "GlycoEpitope Epitopes",
    GD = "Glydin",
    GDB = "Glydin - BiOligo",
    GDV = "Glydin - Cermav",
    GDC = "Glydin - Cummings",
    GDH = "Glydin - Hayes",
    GDSB = "Glydin - SugarBind",
    UCM = "UniCarbKB Motifs",
    GM = "All Motifs"
  )

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
  expect_setequal(info$source_id, names(expected_sources))
  purrr::iwalk(expected_sources, function(source, source_id) {
    collection_info <- info[info$source_id == source_id, ]
    expect_gt(nrow(collection_info), 0)
    expect_equal(unique(collection_info$source), source)
  })
  expect_s3_class(info$glycan_structure, "glyrepr_structure")
  expect_equal(length(info$glycan_structure), nrow(info))
})


test_that("non-GGM motifs use IUPAC and alignment as names", {
  info <- db_motif_info()
  non_ggm_info <- info[info$source_id != "GGM", ]

  expect_false(any(is.na(non_ggm_info$name)))
  expect_identical(
    non_ggm_info$name,
    paste0(
      as.character(non_ggm_info$glycan_structure),
      "#",
      non_ggm_info$alignment
    )
  )
})


test_that("database motifs are unique within source by structure and alignment", {
  info <- db_motif_info()
  motif_keys <- paste(
    info$source_id,
    as.character(info$glycan_structure),
    info$alignment,
    sep = "\r"
  )

  expect_false(any(duplicated(motif_keys)))
})
