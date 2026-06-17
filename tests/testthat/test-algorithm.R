test_that("core alignment root filter rejects incompatible cores", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)GlcNAc(a1-")
  generic_glycan <- glyrepr::convert_to_generic(glycan)
  generic_motif <- glyrepr::convert_to_generic(
    glyparse::parse_iupac_condensed("Gal(b1-6)Gal(a1-")
  )

  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)
  generic_glycan_graph <- glyrepr::get_structure_graphs(generic_glycan)
  generic_motif_graph <- glyrepr::get_structure_graphs(generic_motif)

  expect_false(core_alignment_root_can_match(
    glycan_graph,
    motif_graph,
    alignment = "core",
    strict_sub = TRUE
  ))
  expect_true(core_alignment_root_can_match(
    generic_glycan_graph,
    generic_motif_graph,
    alignment = "core",
    strict_sub = TRUE
  ))
  expect_true(core_alignment_root_can_match(
    glycan_graph,
    motif_graph,
    alignment = "substructure",
    strict_sub = TRUE
  ))
})

test_that("core alignment root mismatch bypasses VF2", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-6)GlcNAc(a1-")

  testthat::local_mocked_bindings(
    perform_vf2 = function(...) {
      stop("VF2 should not run")
    },
    .package = "glymotif"
  )

  expect_false(have_motif(glycan, motif, alignment = "core"))
  expect_identical(count_motif(glycan, motif, alignment = "core"), 0L)
  expect_equal(match_motif(glycan, motif, alignment = "core"), list(list()))
})
