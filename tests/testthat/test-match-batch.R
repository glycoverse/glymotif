test_that("batch motif functions restore repeated and missing glycans", {
  glycan_13 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  glycan_14 <- glyparse::parse_iupac_condensed("Man(b1-4)GlcNAc(a1-")
  glycans <- c(
    glycan_13,
    glycan_13,
    glyrepr::glycan_structure(NA),
    glycan_14
  )
  names(glycans) <- c("g1", "g1-copy", "missing", "g2")
  motifs <- glyparse::parse_iupac_condensed(c(
    m13 = "Hex(b1-3)HexNAc(?1-",
    m14 = "Hex(b1-4)HexNAc(?1-"
  ))

  expected_have <- matrix(
    c(TRUE, TRUE, NA, FALSE, FALSE, FALSE, NA, TRUE),
    ncol = 2,
    dimnames = list(names(glycans), names(motifs))
  )
  expected_count <- matrix(
    c(1L, 1L, NA_integer_, 0L, 0L, 0L, NA_integer_, 1L),
    ncol = 2,
    dimnames = list(names(glycans), names(motifs))
  )
  expected_match <- list(
    m13 = list(
      g1 = list(c(1L, 2L)),
      `g1-copy` = list(c(1L, 2L)),
      missing = NULL,
      g2 = list()
    ),
    m14 = list(
      g1 = list(),
      `g1-copy` = list(),
      missing = NULL,
      g2 = list(c(1L, 2L))
    )
  )

  expect_identical(have_motifs(glycans, motifs), expected_have)
  expect_identical(count_motifs(glycans, motifs), expected_count)
  expect_identical(match_motifs(glycans, motifs), expected_match)
})

test_that("batch motif functions preserve empty glycan shapes", {
  glycans <- glyrepr::glycan_structure()
  motifs <- glyparse::parse_iupac_condensed(c(
    m1 = "Gal(b1-",
    m2 = "GlcNAc(b1-"
  ))

  expect_identical(dim(have_motifs(glycans, motifs)), c(0L, 2L))
  expect_identical(dim(count_motifs(glycans, motifs)), c(0L, 2L))
  expect_identical(
    match_motifs(glycans, motifs),
    list(m1 = list(), m2 = list())
  )
})

test_that("batch motif functions keep per-motif degree constraints", {
  glycans <- glyparse::parse_iupac_condensed(
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  )
  motifs <- glyparse::parse_iupac_condensed(c(
    gal = "Gal(b1-3)GalNAc(a1-",
    glcnac = "GlcNAc(b1-6)GalNAc(a1-"
  ))
  match_degree <- list(c(FALSE, FALSE), c(FALSE, TRUE))

  expect_identical(
    unname(have_motifs(glycans, motifs, match_degree = match_degree)),
    matrix(c(TRUE, FALSE), nrow = 1)
  )
  expect_identical(
    unname(count_motifs(glycans, motifs, match_degree = match_degree)),
    matrix(c(1L, 0L), nrow = 1)
  )
})
