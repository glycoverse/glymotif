test_that(".g_* motif functions match high-level scalar behavior", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)GalNAc(b1-",
    "Gal(b1-4)GalNAc(b1-"
  ))
  motif <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(b1-")
  glycan_graphs <- glyrepr::get_structure_graphs(glycans)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  expect_identical(
    purrr::map_lgl(glycan_graphs, .g_have_motif, motif_graph = motif_graph),
    unname(have_motif(glycans, motif))
  )
  expect_identical(
    purrr::map_int(glycan_graphs, .g_count_motif, motif_graph = motif_graph),
    unname(count_motif(glycans, motif))
  )
  expect_equal(
    purrr::map(glycan_graphs, .g_match_motif, motif_graph = motif_graph),
    unname(match_motif(glycans, motif))
  )
})

test_that(".g_* motif functions support matching options", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-4)GalNAc(b1-")
  motif <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(b1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  expect_false(.g_have_motif(glycan_graph, motif_graph))
  expect_true(.g_have_motif(glycan_graph, motif_graph, ignore_linkages = TRUE))
  expect_identical(.g_count_motif(glycan_graph, motif_graph), 0L)
  expect_identical(
    .g_count_motif(glycan_graph, motif_graph, ignore_linkages = TRUE),
    1L
  )
  expect_equal(
    .g_match_motif(glycan_graph, motif_graph, ignore_linkages = TRUE),
    list(c(1L, 2L))
  )
})

test_that("optional .g_* motif arguments must be named", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  motif <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  expect_error(
    .g_have_motif(glycan_graph, motif_graph, "whole"),
    "must be empty"
  )
  expect_error(
    .g_count_motif(glycan_graph, motif_graph, "whole"),
    "must be empty"
  )
  expect_error(
    .g_match_motif(glycan_graph, motif_graph, "whole"),
    "must be empty"
  )

  expect_true(.g_have_motif(glycan_graph, motif_graph, alignment = "whole"))
  expect_identical(
    .g_count_motif(glycan_graph, motif_graph, alignment = "whole"),
    1L
  )
  expect_equal(
    .g_match_motif(glycan_graph, motif_graph, alignment = "whole"),
    list(c(1L, 2L))
  )
})

test_that(".g_* motif functions keep mono-type compatibility as caller contract", {
  glycan <- glyrepr::as_glycan_structure("Man(?1-")
  generic_motif <- glyrepr::as_glycan_structure("Hex(?1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  generic_motif_graph <- glyrepr::get_structure_graphs(generic_motif)

  expect_false(.g_have_motif(glycan_graph, generic_motif_graph))
  expect_identical(.g_count_motif(glycan_graph, generic_motif_graph), 0L)
  expect_equal(.g_match_motif(glycan_graph, generic_motif_graph), list())

  igraph::V(glycan_graph)$mono <- glyrepr::convert_to_generic(
    igraph::V(glycan_graph)$mono
  )
  expect_true(.g_have_motif(glycan_graph, generic_motif_graph))
})

test_that(".g_* motif functions recycle scalar match_degree like high-level API", {
  glycan <- glyparse::parse_iupac_condensed(
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  )
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  expect_identical(
    .g_have_motif(glycan_graph, motif_graph, match_degree = TRUE),
    unname(have_motif(glycan, motif, match_degree = TRUE))
  )
  expect_identical(
    .g_count_motif(glycan_graph, motif_graph, match_degree = TRUE),
    unname(count_motif(glycan, motif, match_degree = TRUE))
  )
  expect_equal(
    .g_match_motif(glycan_graph, motif_graph, match_degree = TRUE),
    unname(match_motif(glycan, motif, match_degree = TRUE)[[1]])
  )
})
