test_that("native VF2 prunes mixed residues and incompatible linkages", {
  glycan <- glyparse::parse_iupac_condensed(
    "Gal(b1-3)GalNAc(a1-3)GlcNAc(?1-"
  )
  mixed_motif <- glyparse::parse_iupac_condensed(
    "Hex(b1-3)GalNAc(a1-"
  )
  wrong_linkage <- glyparse::parse_iupac_condensed(
    "Hex(b1-4)GalNAc(a1-"
  )
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  mixed_graph <- glyrepr::get_structure_graphs(mixed_motif)
  wrong_graph <- glyrepr::get_structure_graphs(wrong_linkage)

  expect_identical(
    perform_compatible_vf2(
      glycan_graph,
      mixed_graph,
      alignment = "substructure",
      ignore_linkages = FALSE,
      strict_sub = TRUE,
      match_degree = NULL,
      mode = "strict"
    ),
    list(c(1L, 2L))
  )
  expect_identical(
    perform_compatible_vf2(
      glycan_graph,
      wrong_graph,
      alignment = "substructure",
      ignore_linkages = FALSE,
      strict_sub = TRUE,
      match_degree = NULL,
      mode = "strict"
    ),
    list()
  )
})


test_that("native VF2 can stop after the first compatible mapping", {
  glycan <- glyparse::parse_iupac_condensed(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(?1-"
  )
  motif <- glyparse::parse_iupac_condensed("Man(?1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)
  match_all <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = NULL,
    mode = "strict"
  )
  match_first <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = NULL,
    mode = "strict",
    first_only = TRUE
  )

  expect_length(match_all, 3L)
  expect_length(match_first, 1L)
  expect_in(match_first[[1]], match_all)
})


test_that("native VF2 deduplicates motif automorphisms by mapped vertices", {
  glycan <- glyparse::parse_iupac_condensed(
    "Man(a1-3)[Man(a1-6)]Man(b1-"
  )
  motif <- glyparse::parse_iupac_condensed(
    "Hex(?1-?)[Hex(?1-?)]Hex(?1-"
  )
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  result <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "whole",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = NULL,
    mode = "strict"
  )

  expect_identical(length(result), 1L)
  expect_setequal(result[[1]], seq_len(igraph::vcount(glycan_graph)))
})


test_that("native VF2 preserves legacy mapping order", {
  glycan <- glyparse::parse_iupac_condensed(
    "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(?1-"
  )
  motif <- glyparse::parse_iupac_condensed("Hex(?1-?)Hex(?1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)
  colors <- colorize_graphs(glycan_graph, motif_graph)
  candidates <- perform_vf2(
    glycan_graph,
    motif_graph,
    colors$glycan_colors,
    colors$motif_colors
  )
  context <- prepare_validation_context(
    glycan_graph,
    motif_graph,
    alignment = "substructure",
    ignore_linkages = FALSE
  )
  valid_mask <- vapply(
    candidates,
    function(candidate) {
      is_valid_result(
        candidate,
        glycan = glycan_graph,
        motif = motif_graph,
        alignment = "substructure",
        ignore_linkages = FALSE,
        strict_sub = TRUE,
        mode = "strict",
        context = context
      )
    },
    logical(1L)
  )
  valid <- candidates[valid_mask]
  expected <- unique_vf2_res(valid)

  result <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "substructure",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = NULL,
    mode = "strict"
  )

  expect_identical(result, expected)
})


test_that("native VF2 ignores alignment when degree matching is supplied", {
  glycan <- glyparse::parse_iupac_condensed(
    "GlcNAc(b1-6)Gal(b1-6)Gal(a1-"
  )
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)Gal(a1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)

  terminal <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "terminal",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = NULL,
    mode = "strict"
  )
  degree_controlled <- perform_compatible_vf2(
    glycan_graph,
    motif_graph,
    alignment = "terminal",
    ignore_linkages = FALSE,
    strict_sub = TRUE,
    match_degree = c(FALSE, FALSE),
    mode = "strict"
  )

  expect_length(terminal, 0L)
  expect_identical(degree_controlled, list(c(2L, 3L)))
})
