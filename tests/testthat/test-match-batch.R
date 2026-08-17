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

test_that("batch motif functions match heterogeneous vectors element-wise", {
  glycans <- glyrepr::as_glycan_structure(c(
    concrete = "Gal(b1-3)GalNAc(a1-",
    mixed = "Hex(b1-3)GalNAc(a1-",
    generic = "Hex(b1-3)HexNAc(a1-",
    missing = NA_character_,
    concrete_copy = "Gal(b1-3)GalNAc(a1-"
  ))
  motifs <- glyrepr::as_glycan_structure(c(
    generic = "Hex(b1-3)HexNAc(a1-",
    concrete = "Gal(b1-3)GalNAc(a1-",
    mixed = "Hex(b1-3)GalNAc(a1-"
  ))
  expected_have <- rbind(
    concrete = c(generic = TRUE, concrete = TRUE, mixed = TRUE),
    mixed = c(generic = TRUE, concrete = FALSE, mixed = TRUE),
    generic = c(generic = TRUE, concrete = FALSE, mixed = FALSE),
    missing = c(generic = NA, concrete = NA, mixed = NA),
    concrete_copy = c(generic = TRUE, concrete = TRUE, mixed = TRUE)
  )
  expected_count <- expected_have
  storage.mode(expected_count) <- "integer"

  matches <- match_motifs(glycans, motifs)
  match_counts <- vapply(
    matches,
    function(motif_matches) {
      vapply(
        motif_matches,
        \(x) if (is.null(x)) NA_integer_ else length(x),
        integer(1L)
      )
    },
    integer(length(glycans))
  )
  dimnames(match_counts) <- dimnames(expected_count)

  expect_identical(have_motifs(glycans, motifs), expected_have)
  expect_identical(count_motifs(glycans, motifs), expected_count)
  expect_identical(match_counts, expected_count)
  expect_identical(names(matches), names(motifs))
  expect_identical(
    unname(lapply(matches, names)),
    rep(list(names(glycans)), length(motifs))
  )

  reversed <- motifs[3:1]
  expect_identical(
    have_motifs(glycans, reversed),
    expected_have[, 3:1, drop = FALSE]
  )
  expect_identical(
    count_motifs(glycans, reversed),
    expected_count[, 3:1, drop = FALSE]
  )
  expect_identical(match_motifs(glycans, reversed), matches[3:1])
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

test_that("batch profiles keep strict fuzzy composition keys", {
  glycans <- glyparse::parse_iupac_condensed(c(
    match = "Glc3Me6S(a1-",
    miss = "Glc4Me3S(a1-"
  ))
  motifs <- glyparse::parse_iupac_condensed(c(
    fuzzy = "Glc?Me6S(a1-",
    exact = "Glc3Me6S(a1-"
  ))

  batch <- prepare_match_batch(glycans, motifs)

  expect_identical(
    unname(purrr::map_chr(
      batch$motif_profiles,
      \(profile) profile$composition_profile$key_mode
    )),
    c("base", "exact")
  )
  expect_identical(
    have_motifs(glycans, motifs),
    matrix(
      c(TRUE, FALSE, TRUE, FALSE),
      nrow = 2,
      dimnames = list(names(glycans), names(motifs))
    )
  )
})
