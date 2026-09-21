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
