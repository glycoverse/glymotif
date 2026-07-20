test_that("integer graph helpers preserve igraph results", {
  structure <- glyparse::auto_parse(
    "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
  )
  graph <- glyrepr::get_structure_graphs(structure)

  expect_identical(
    graph_vertex_attr(graph, "mono"),
    igraph::vertex_attr(
      graph,
      "mono",
      index = seq_len(igraph::vcount(graph))
    )
  )
  expect_identical(
    graph_edge_attr(graph, "linkage"),
    igraph::edge_attr(
      graph,
      "linkage",
      index = seq_len(igraph::ecount(graph))
    )
  )
  expect_identical(
    graph_degree(graph, mode = "out"),
    unname(igraph::degree(graph, mode = "out"))
  )
  expect_identical(
    graph_neighbor_ids(graph, 3L, mode = "out"),
    as.integer(igraph::neighbors(graph, 3L, mode = "out"))
  )
  expect_identical(
    graph_subcomponent_ids(graph, 3L, mode = "out"),
    as.integer(igraph::subcomponent(graph, 3L, mode = "out"))
  )
})

test_that("integer VF2 mappings preserve mapping order", {
  glycan <- glyparse::auto_parse(
    "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
  )
  motif <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")
  glycan_graph <- glyrepr::get_structure_graphs(glycan)
  motif_graph <- glyrepr::get_structure_graphs(motif)
  colors <- colorize_graphs(glycan_graph, motif_graph)

  expected <- igraph::graph.get.subisomorphisms.vf2(
    glycan_graph,
    motif_graph,
    vertex.color1 = colors$glycan_colors,
    vertex.color2 = colors$motif_colors
  )
  expected <- lapply(expected, as.integer)

  expect_identical(
    perform_vf2(
      glycan_graph,
      motif_graph,
      colors$glycan_colors,
      colors$motif_colors
    ),
    expected
  )
})

test_that("structure levels preserve glyrepr classifications", {
  structures <- list(
    intact = glyparse::auto_parse("Gal(b1-3)GalNAc(a1-"),
    partial = glyparse::auto_parse("Gal(b1-?)GalNAc(a1-"),
    topological = glyparse::auto_parse("Gal(??-?)GalNAc(??-"),
    basic = glyparse::auto_parse("Hex(??-?)HexNAc(??-")
  )

  expected <- vapply(
    structures,
    \(structure) suppressWarnings(glyrepr::get_structure_level(structure)),
    character(1L)
  )
  actual <- vapply(structures, structure_level, character(1L))

  expect_identical(actual, expected)
  expect_identical(
    structure_level(c(structures$intact, structures$topological)),
    suppressWarnings(
      glyrepr::get_structure_level(
        c(structures$intact, structures$topological)
      )
    )
  )
  expect_identical(
    structure_level(glyrepr::glycan_structure()),
    NA_character_
  )
})

test_that("repeated matching does not retain graph weak references", {
  glycans <- glyparse::auto_parse(
    c(
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-",
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
    )
  )
  motifs <- glyparse::auto_parse(
    c(
      "Gal(b1-4)GlcNAc(b1-",
      "Man(a1-6)Man(b1-"
    )
  )
  repeat_matching <- function(n) {
    for (i in seq_len(n)) {
      have_motif(glycans, motifs[[1]])
      have_motifs(glycans, motifs)
      count_motif(glycans, motifs[[1]])
      match_motif(glycans, motifs[[1]])
    }
  }

  repeat_matching(20L)
  before <- gc(full = TRUE)[, "used"]
  repeat_matching(200L)
  after <- gc(full = TRUE)[, "used"]
  retained <- after - before

  expect_lt(unname(retained[["Ncells"]]), 1000)
  expect_lt(unname(retained[["Vcells"]]), 4000)
})
