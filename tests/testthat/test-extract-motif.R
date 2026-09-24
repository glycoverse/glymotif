test_that("extract_motif works for intact glycans", {
  glycan <- "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
  result <- extract_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Gal(b1-",
    "GlcNAc(a1-",
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(a1-6)GalNAc(a1-",
    "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
  ))
  expect_setequal(as.character(result), as.character(expected))
})

test_that("extract_motif works for topological glycans", {
  glycan <- "Hex(??-?)[HexNAc(??-?)]HexNAc(??-"
  result <- extract_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Hex(??-",
    "HexNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "HexNAc(??-?)HexNAc(??-",
    "Hex(??-?)[HexNAc(??-?)]HexNAc(??-"
  ))
  expect_setequal(as.character(result), as.character(expected))
})

test_that("extract_motif works for multiple glycans", {
  glycan <- c(
    "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  )
  result <- extract_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Gal(b1-",
    "GlcNAc(a1-",
    "GalNAc(a1-",
    "Neu5Ac(a2-",
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(a1-6)GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-",
    "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  ))
  expect_setequal(as.character(result), as.character(expected))
})

test_that("extract_motif preserves distinct residue and linkage details", {
  glycans <- c(
    "Gal(a1-3)GlcNAc(a1-",
    "Gal(b1-3)GlcNAc(a1-",
    "Gal3S(b1-3)GlcNAc(a1-",
    "Gal6S(b1-3)GlcNAc(a1-",
    "Gal(b1-4)GlcNAc(a1-"
  )

  result <- extract_motif(glycans, max_size = 2)

  expect_length(result, 10)
})

test_that("extract_motif keeps linkages paired with their branches", {
  glycans <- c(
    "Gal(b1-3)[Fuc(a1-2)]GlcNAc(a1-",
    "Gal(a1-2)[Fuc(b1-3)]GlcNAc(a1-"
  )

  result <- extract_motif(glycans, max_size = 3)

  expect_length(result, 11)
})

test_that("extract_motif's max_size works", {
  glycan <- "Glc(a1-3)Glc(a1-3)Glc(a1-3)Glc(a1-"
  result <- extract_motif(glycan, max_size = 3)
  expected <- glyrepr::as_glycan_structure(c(
    "Glc(a1-",
    "Glc(a1-3)Glc(a1-",
    "Glc(a1-3)Glc(a1-3)Glc(a1-"
  ))
  expect_setequal(as.character(result), as.character(expected))
})

test_that("optional extract_motif arguments must be named", {
  glycan <- "Gal(b1-3)GalNAc(a1-"

  expect_error(
    extract_motif(glycan, 2),
    "must be empty"
  )
  expect_no_error(extract_motif(glycan, max_size = 2))
})

test_that("extract_branch_motif works for intact glycan structures", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  res <- extract_branch_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif works for topological glycan structures", {
  glycans <- "NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-"
  res <- extract_branch_motif(glycans)
  expected <- glyrepr::as_glycan_structure(c(
    "NeuAc(??-?)Hex(??-?)HexNAc(??-",
    "Hex(??-?)HexNAc(??-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif works for multiple glycans", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  res <- extract_branch_motif(glycans)
  expected <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif works for glycans with bisecting GlcNAc", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  res <- extract_branch_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif works for tri-antennary glycans", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  res <- extract_branch_motif(glycan)
  expected <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif works for glycans with complex branching patterns", {
  glycan <- "Neu5Ac(a2-3)[Fuc(a1-2)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  res <- extract_branch_motif(glycan)
  expected <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)[Fuc(a1-2)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif warns for glycans other than N-glycans", {
  paucimannose <- "Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_snapshot_warning(extract_branch_motif(paucimannose))
  expect_length(suppressWarnings(extract_branch_motif(paucimannose)), 0)
})

test_that("optional extract_branch_motif arguments must be named", {
  glycan <- glyrepr::n_glycan_core()

  expect_error(
    extract_branch_motif(glycan, TRUE),
    "must be empty"
  )
  expect_no_error(extract_branch_motif(glycan, including_core = TRUE))
})

test_that("extract_branch_motif with including_core works for topological glycans", {
  glycan <- "HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-"
  res <- extract_branch_motif(glycan, including_core = TRUE)
  expected <- glyrepr::as_glycan_structure(c(
    "HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Hex(??-?)HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif with including_core works for intact glycans", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  res <- extract_branch_motif(glycan, including_core = TRUE)
  expected <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "Gal(b1-4)GlcNAc(b1-?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif with including_core works for glycans without linkages", {
  glycan <- "GlcNAc(??-?)Man(??-?)[Gal(??-?)GlcNAc(??-?)Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  res <- extract_branch_motif(glycan, including_core = TRUE)
  expected <- glyrepr::as_glycan_structure(c(
    "GlcNAc(??-?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "Gal(??-?)GlcNAc(??-?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  ))
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif appends cores element-wise for mixed vectors", {
  concrete <- glyrepr::as_glycan_structure(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  generic <- glyrepr::convert_to_generic(concrete)
  mixed <- glyrepr::as_glycan_structure(
    "Hex(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )

  result <- extract_branch_motif(
    c(concrete = concrete, generic = generic, mixed = mixed),
    including_core = TRUE
  )
  expected <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-?)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "Hex(b1-4)HexNAc(b1-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Hex(b1-4)GlcNAc(b1-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"
  ))

  expect_identical(as.character(result), as.character(expected))
  expect_identical(
    glyrepr::get_mono_type(result),
    c("concrete", "generic", "mixed")
  )
})

test_that("native extraction preserves ordered R enumeration and input graphs", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)[Fuc(a1-2)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-",
    "Gal3S(b1-3)[Gal6S(b1-6)]GlcNAc(a1-",
    "Gal(a1-3)[Gal(b1-6)]GlcNAc(b1-",
    "Hex(??-?)[GlcNAc(b1-6)]HexNAc(??-"
  ))
  graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
  before <- serialize(graphs, NULL)
  for (size in c(0, 1, 2, 2.5, 3, 5, Inf)) {
    seen <- new.env(hash = TRUE, parent = emptyenv())
    expected <- unlist(
      lapply(graphs, function(g) {
        .extract_motifs_from_graph(g, size, seen)
      }),
      recursive = FALSE
    )
    expected <- unique(glyrepr::as_glycan_structure(expected))
    actual <- extract_motif(glycans, max_size = size)
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(names(actual), names(expected))
  }
  expect_identical(serialize(graphs, NULL), before)
})

test_that("native branch extraction deduplicates complete subtrees in order", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal3S(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
  matches <- match_motif(
    glycans,
    glyrepr::as_glycan_structure(
      "HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"
    ),
    alignment = "core"
  )
  expected <- unique(glyrepr::as_glycan_structure(
    .extract_branch_subtrees_r(graphs, matches)
  ))
  expect_identical(
    as.character(extract_branch_motif(glycans)),
    as.character(expected)
  )
  expect_length(expected, 2L)
})

test_that("native extraction handles empty inputs and singleton structures", {
  expect_length(extract_motif(glyrepr::glycan_structure()), 0L)
  expect_length(extract_branch_motif(glyrepr::glycan_structure()), 0L)
  expect_identical(as.character(extract_motif("Gal(b1-")), "Gal(b1-")
})

test_that("native enumeration prunes wide branches before forming products", {
  g <- igraph::make_star(41, mode = "out")
  g <- igraph::set_vertex_attr(g, "mono", value = c("GlcNAc", rep("Gal", 40)))
  g <- igraph::set_vertex_attr(g, "sub", value = rep("", 41))
  g <- igraph::set_edge_attr(g, "linkage", value = rep("??-?", 40))
  g$anomer <- "??"
  motifs <- .extract_motif_candidates_native(list(g), list(1L), 3)
  expect_identical(vapply(motifs, igraph::vcount, numeric(1)), c(1, 2, 3))
})

test_that("floating branch attributes retain the reference extraction path", {
  glycan <- glyrepr::as_glycan_structure(
    "{Fuc(a1-?)}Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  graphs <- glyrepr::get_structure_graphs(glycan, return_list = TRUE)
  expect_identical(.has_floating_motif_parts(graphs[[1]]), TRUE)
  expect_snapshot(error = TRUE, extract_branch_motif(glycan))
})
