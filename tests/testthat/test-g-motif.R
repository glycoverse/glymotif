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

test_that(".g_* motif functions support generic and mixed residue rules", {
  concrete <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  generic <- glyrepr::as_glycan_structure("Hex(b1-3)HexNAc(a1-")
  mixed_motif <- glyrepr::as_glycan_structure("Hex(b1-3)GalNAc(a1-")
  concrete_graph <- glyrepr::get_structure_graphs(concrete)
  generic_graph <- glyrepr::get_structure_graphs(generic)
  mixed_motif_graph <- glyrepr::get_structure_graphs(mixed_motif)

  expect_identical(.g_have_motif(concrete_graph, mixed_motif_graph), TRUE)
  expect_identical(.g_count_motif(concrete_graph, mixed_motif_graph), 1L)
  expect_equal(
    .g_match_motif(concrete_graph, mixed_motif_graph),
    list(c(1L, 2L))
  )
  expect_identical(.g_have_motif(generic_graph, mixed_motif_graph), FALSE)
  expect_identical(
    .g_have_motif(generic_graph, mixed_motif_graph, mode = "lenient"),
    TRUE
  )
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
