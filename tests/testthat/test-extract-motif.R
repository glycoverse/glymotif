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