test_that("native matching agrees with graph matching and leaves inputs intact", {
  glycans <- glyrepr::as_glycan_structure(c(
    "{Gal(a1-3)|2,3}Man(a1-3)[Glc(a1-6)]GlcNAc(b1-",
    "{6S}Gal(a1-3)Glc(a1-",
    "Gal(b1-3)GalNAc(a1-"
  ))
  motifs <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)HexNAc(??-",
    "Gal6S(a1-"
  ))
  gg <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
  mg <- glyrepr::get_structure_graphs(motifs, return_list = TRUE)
  for (mode in c("strict", "lenient")) {
    for (sf in c(TRUE, FALSE)) {
      expected <- lapply(mg, function(m) {
        vapply(
          gg,
          function(g) {
            .g_count_motif(g, m, mode = mode, strict_floating = sf)
          },
          integer(1L)
        )
      })
      expect_identical(
        unname(count_motifs(
          glycans,
          motifs,
          mode = mode,
          strict_floating = sf
        )),
        unname(do.call(cbind, expected))
      )
    }
    expected <- lapply(mg, function(m) {
      unname(lapply(gg, function(g) .g_match_motif(g, m, mode = mode)))
    })
    expect_identical(
      unname(lapply(match_motifs(glycans, motifs, mode = mode), unname)),
      unname(expected)
    )
  }
  before <- serialize(list(glycans, motifs), NULL)
  invisible(match_motifs(glycans, motifs))
  expect_identical(serialize(list(glycans, motifs), NULL), before)
})

test_that("native calls retain missing glycans and zero-row matrix shapes", {
  glycans <- glyrepr::glycan_structure(NA_character_)
  expect_identical(have_motif(glycans, "Man"), NA)
  expect_identical(count_motif(glycans, "Man"), NA_integer_)
  expect_identical(
    match_motif(glycans, glyrepr::as_glycan_structure("Man")),
    list(NULL)
  )
  expect_identical(dim(have_motifs(glycans[integer()], "Man")), c(0L, 1L))
})

test_that("floating combination limits remain enforced by native matching", {
  glycan <- native_structure_input(glyrepr::as_glycan_structure("Gal(a1-3)Glc"))
  motif <- native_structure_input(glyrepr::as_glycan_structure("Gal"))
  run <- function(n) {
    glycan$graphs[[1]]$attributes$floating_substituents <- rep(
      list(list(parents = 1:2, substituent = "?S")),
      n
    )
    cpp_match_structures(
      glycan,
      motif,
      native_monosaccharide_dictionary(),
      "substructure",
      FALSE,
      FALSE,
      FALSE,
      list(NULL),
      "have",
      TRUE,
      256L
    )
  }
  expect_identical(run(8L), matrix(TRUE, 1L, 1L))
  expect_snapshot(error = TRUE, run(9L))
})

test_that("alditol and ordinary residues match only their own form in every alignment", {
  forms <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)GalNAc-ol(a1-")
  graphs <- glyrepr::get_structure_graphs(
    glyrepr::as_glycan_structure(forms),
    return_list = TRUE
  )
  # Legacy graphs omit the alditol attribute.
  graphs[[1]] <- igraph::delete_graph_attr(graphs[[1]], "alditol")
  for (alignment in c("substructure", "core", "terminal", "whole")) {
    expect_identical(
      unname(have_motifs(forms, forms, alignments = alignment)),
      diag(2) == 1
    )
    for (g in seq_along(graphs)) {
      for (m in seq_along(graphs)) {
        expect_identical(
          .g_have_motif(graphs[[g]], graphs[[m]], alignment = alignment),
          g == m
        )
      }
    }
  }
})

test_that("only the main reducing end matches alditol motifs", {
  glycan <- glyrepr::as_glycan_structure(
    "{GalNAc(a1-3)|2,3}Man(a1-3)[Glc(a1-6)]GalNAc-ol(b1-"
  )
  motifs <- glyrepr::as_glycan_structure(c("HexNAc", "HexNAc-ol"))
  graph <- glyrepr::get_structure_graphs(glycan)
  motif_graphs <- glyrepr::get_structure_graphs(motifs, return_list = TRUE)
  # The floating GalNAc (node 1) can attach to either branch. Node 4 is the root.
  for (i in seq_along(motifs)) {
    node <- c(1L, 4L)[[i]]
    args <- list(
      ignore_linkages = TRUE,
      strict_sub = FALSE,
      mode = "lenient",
      match_degree = FALSE
    )
    expect_identical(
      do.call(match_motif, c(list(glycan, motifs[i]), args)),
      list(list(node))
    )
    expect_identical(
      do.call(.g_match_motif, c(list(graph, motif_graphs[[i]]), args)),
      list(node)
    )
    expect_identical(count_motif(glycan, motifs[i]), 1L)
  }
})
