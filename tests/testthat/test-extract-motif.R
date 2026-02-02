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
  expected <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)[Fuc(a1-2)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_setequal(as.character(res), as.character(expected))
})

test_that("extract_branch_motif rejects glycans other than N-glycans", {
  glycan <- "Gal(b1-3)GalNAc(a1-"
  expect_error(extract_branch_motif(glycan), "must be N-glycans")
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
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "Gal(b1-4)GlcNAc(b1-2)Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-"
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