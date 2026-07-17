test_that("lenient fuzzy motifs bypass generic residue pruning", {
  motif <- glyparse::parse_iupac_condensed(
    "Gal?NAc(b1-3)Neu5Ac(a2-"
  )
  glycan <- glyrepr::convert_to_generic(
    glyparse::parse_iupac_condensed("GalNAc(b1-3)Neu5Ac(a2-")
  )

  expect_true(have_motif(glycan, motif, mode = "lenient"))
  expect_identical(count_motif(glycan, motif, mode = "lenient"), 1L)
  expect_equal(
    match_motif(glycan, motif, mode = "lenient"),
    list(list(c(1L, 2L)))
  )

  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)
  expect_true(.g_have_motif(glycan_graph, motif_graph, mode = "lenient"))
  expect_identical(
    .g_count_motif(glycan_graph, motif_graph, mode = "lenient"),
    1L
  )
  expect_equal(
    .g_match_motif(glycan_graph, motif_graph, mode = "lenient"),
    list(c(1L, 2L))
  )
})

test_that("plural lenient matching uses graph-local generic residue keys", {
  glycans <- glyrepr::convert_to_generic(
    glyparse::parse_iupac_condensed(c(
      "GalNAc(b1-3)Neu5Ac(a2-",
      "Gal(b1-3)Neu5Ac(a2-",
      "Glc(b1-4)GlcNAc(a1-"
    ))
  )
  motifs <- glyparse::parse_iupac_condensed(c(
    "Gal?NAc(b1-3)Neu5Ac(a2-",
    "Gal(b1-3)Neu5Ac(a2-",
    "Man(b1-4)GlcNAc(a1-"
  ))

  expected_have <- diag(TRUE, 3)
  expected_count <- diag(1L, 3)
  expect_identical(
    unname(have_motifs(glycans, motifs, mode = "lenient")),
    expected_have
  )
  expect_identical(
    unname(count_motifs(glycans, motifs, mode = "lenient")),
    expected_count
  )
  expect_equal(
    match_motifs(glycans, motifs, mode = "lenient"),
    list(
      list(list(c(1L, 2L)), list(), list()),
      list(list(), list(c(1L, 2L)), list()),
      list(list(), list(), list(c(1L, 2L)))
    )
  )

  batch <- prepare_match_batch(glycans, motifs, mode = "lenient")
  expect_identical(
    purrr::map_chr(batch$motif_profiles, "key_mode"),
    c("none", "generic", "generic")
  )
})
