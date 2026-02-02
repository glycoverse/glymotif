# Load required packages
library(glyparse)

# Check if two match results are the same, ignoring orders.
expect_match_motif_equal <- function(result, expected) {
  # Process all match results of one glycan
  process_motif_glycan_match <- function(res) {
    sort(purrr::map_chr(res, ~ paste(sort(.x), collapse = ",")))
  }
  # Process all match results of all glycans (return value of `match_motif()`)
  process_match_res <- function(res) {
    purrr::map(res, process_motif_glycan_match)
  }
  processed_result <- process_match_res(result)
  processed_expected <- process_match_res(expected)
  expect_equal(processed_result, processed_expected)
}

expect_match_motifs_equal <- function(result, expected) {
  # Process all match results of one motif to one glycan
  process_motif_glycan_match <- function(res) {
    sort(purrr::map_chr(res, ~ paste(sort(.x), collapse = ",")))
  }
  # Process all match results of one motif to all glycans
  process_motif_match_res <- function(res) {
    purrr::map(res, process_motif_glycan_match)
  }
  # Process all match results of all motifs to all glycans (return value of `match_motifs()`)
  process_motifs_match_res <- function(res) {
    purrr::map(res, process_motif_match_res)
  }
  processed_result <- process_motifs_match_res(result)
  processed_expected <- process_motifs_match_res(expected)
  expect_equal(processed_result, processed_expected)
}

# ========== `match_motif()` ==========
test_that("match_motif works", {
  glycan <- glyrepr::n_glycan_core()
  motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
  result <- match_motif(glycan, motif)
  expect_match_motif_equal(result, list(list(c(1, 2, 3))))
})

test_that("match_motif rejects non-glyrepr_structure objects", {
  glycan <- glyrepr::n_glycan_core()
  motif <- "Man(a1-3)[Man(a1-6)]Man(b1-"
  expect_error(match_motif(glycan, motif), "must be a 'glyrepr_structure' object")

  glycan <- "Man(a1-3)[Man(a1-6)]Man(b1-"
  motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
  expect_error(match_motif(glycan, motif), "must be a 'glyrepr_structure' object")
})

test_that("match_motif return empty list when motif not found", {
  glycan <- glyrepr::n_glycan_core()
  motif <- parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  result <- match_motif(glycan, motif)
  expect_match_motif_equal(result, list(list()))
})

# ========== `match_motifs()` ==========
test_that("match_motifs works", {
  glycan <- glyrepr::n_glycan_core()
  motifs <- c(
    "Man(a1-3)[Man(a1-6)]Man(b1-",
    "Man(a1-?)Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  )
  motifs <- parse_iupac_condensed(motifs)
  result <- match_motifs(glycan, motifs)
  expected <- list(
    # Motif 1
    list(list(c(1, 2, 3))),
    # Motif 2
    list(list(c(1, 3, 4, 5), c(2, 3, 4, 5)))
  )
  expect_match_motifs_equal(result, expected)
})

test_that("match_motifs accepts character vectors and parses them", {
  glycan <- glyrepr::n_glycan_core()
  motifs <- c("Man(a1-3)[Man(a1-6)]Man(b1-", "GlcNAc(b1-4)GlcNAc(?1-")
  # Character vectors should be auto-parsed
  result <- match_motifs(glycan, motifs)
  expect_type(result, "list")
  expect_length(result, 2)

  glycan_str <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  motifs <- parse_iupac_condensed(c("Man(a1-3)[Man(a1-6)]Man(b1-", "GlcNAc(b1-4)GlcNAc(?1-"))
  # Character glycans should be auto-parsed
  result <- match_motifs(glycan_str, motifs)
  expect_type(result, "list")
  expect_length(result, 2)
})

test_that("match_motifs returns empty lists when motifs not found", {
  glycan <- glyrepr::n_glycan_core()
  motifs <- c("Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-")
  motifs <- parse_iupac_condensed(motifs)
  result <- match_motifs(glycan, motifs)
  expected <- list(
    # Motif 1
    list(list()),
    # Motif 2
    list(list())
  )
  expect_match_motifs_equal(result, expected)
})

# ========== Multiple glycans ==========
test_that("match_motif works with multiple glycans", {
  glycans <- c(glyrepr::n_glycan_core(), glyrepr::o_glycan_core_2())
  motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
  result <- match_motif(glycans, motif)
  expected <- list(
    # Glycan 1 (N-glycan core)
    list(c(1, 2, 3)),
    # Glycan 2 (O-glycan core 2) - no match
    list()
  )
  expect_match_motif_equal(result, expected)
})

test_that("match_motifs works with multiple glycans", {
  glycans <- c(glyrepr::n_glycan_core(), glyrepr::o_glycan_core_2())
  motifs <- c("Man(a1-3)[Man(a1-6)]Man(b1-", "Gal(b1-3)GalNAc(?1-")
  motifs <- parse_iupac_condensed(motifs)
  result <- match_motifs(glycans, motifs)
  expected <- list(
    # Motif 1
    list(
      list(c(1, 2, 3)),  # Glycan 1 match
      list()             # Glycan 2 no match
    ),
    # Motif 2
    list(
      list(),            # Glycan 1 no match
      list(c(1, 3))      # Glycan 2 match (Gal=node1, GalNAc=node3)
    )
  )
  expect_match_motifs_equal(result, expected)
})

# ========== Alignment parameter ==========
test_that("match_motif works with whole alignment", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- glyrepr::o_glycan_core_2()
  result <- match_motif(glycan, motif, alignment = "whole")
  expected <- list(list(c(1, 2, 3)))
  expect_match_motif_equal(result, expected)
})

test_that("match_motif returns empty with whole alignment when not identical", {
  glycan <- glyrepr::o_glycan_core_2()
  motif <- parse_iupac_condensed("Gal(b1-3)GalNAc(?1-")
  result <- match_motif(glycan, motif, alignment = "whole")
  expected <- list(list())
  expect_match_motif_equal(result, expected)
})

test_that("match_motif works with core alignment", {
  glycan <- parse_iupac_condensed("GlcNAc(b1-6)Gal(b1-6)Gal(a1-")
  motif <- parse_iupac_condensed("Gal(b1-6)Gal(a1-")
  result <- match_motif(glycan, motif, alignment = "core")
  expected <- list(list(c(2, 3)))  # Gal=node2, Gal=node3 (core is node1)
  expect_match_motif_equal(result, expected)
})

test_that("match_motif works with terminal alignment", {
  glycan <- parse_iupac_condensed("Gal(b1-3)Gal(b1-4)GlcNAc(b1-")
  motif <- parse_iupac_condensed("Gal(b1-3)Gal(b1-")
  result <- match_motif(glycan, motif, alignment = "terminal")
  expected <- list(list(c(1, 2)))  # Terminal Gal=node1, Gal=node2
  expect_match_motif_equal(result, expected)
})

test_that("match_motifs works with different alignments", {
  glycan <- glyrepr::o_glycan_core_2()
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-")
  motifs <- parse_iupac_condensed(motifs)
  alignments <- c("substructure", "whole")
  result <- match_motifs(glycan, motifs, alignments = alignments)
  expected <- list(
    # Motif 1 with substructure alignment
    list(list(c(1, 3))),  # Gal=node1, GalNAc=node3
    # Motif 2 with whole alignment - no match (not identical)
    list(list())
  )
  expect_match_motifs_equal(result, expected)
})

# ========== ignore_linkages parameter ==========
test_that("match_motif works with ignore_linkages = TRUE", {
  glycan <- glyrepr::o_glycan_core_2()  # Gal(b1-3)GalNAc and GlcNAc(b1-6)GalNAc
  motif <- parse_iupac_condensed("Gal(b1-4)GalNAc(?1-")  # Wrong linkage

  # Should not match with default ignore_linkages = FALSE
  result_false <- match_motif(glycan, motif, ignore_linkages = FALSE)
  expect_match_motif_equal(result_false, list(list()))

  # Should match with ignore_linkages = TRUE
  result_true <- match_motif(glycan, motif, ignore_linkages = TRUE)
  expected <- list(list(c(1, 3)))  # Gal=node1, GalNAc=node3
  expect_match_motif_equal(result_true, expected)
})

test_that("match_motifs works with ignore_linkages = TRUE", {
  glycan <- glyrepr::o_glycan_core_2()
  motifs <- c("Gal(b1-4)GalNAc(?1-", "GlcNAc(b1-3)GalNAc(?1-")  # Wrong linkages
  motifs <- parse_iupac_condensed(motifs)

  result <- match_motifs(glycan, motifs, ignore_linkages = TRUE)
  expected <- list(
    # Motif 1
    list(list(c(1, 3))),  # Gal=node1, GalNAc=node3
    # Motif 2
    list(list(c(2, 3)))   # GlcNAc=node2, GalNAc=node3
  )
  expect_match_motifs_equal(result, expected)
})

# ========== Multiple matches ==========
test_that("match_motif finds multiple matches in same glycan", {
  glycan <- parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(a1-")
  motif <- parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  result <- match_motif(glycan, motif)
  expected <- list(list(c(2, 3)))  # Gal=node2, GalNAc=node3 (terminal match)
  expect_match_motif_equal(result, expected)
})

test_that("match_motif finds multiple instances of simple motif", {
  glycan <- parse_iupac_condensed("Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  motif <- parse_iupac_condensed("Man(a1-")
  result <- match_motif(glycan, motif)
  # Should find Man residues at nodes 4 and 5 (terminal Man residues)
  expected <- list(list(c(1), c(2)))  # Terminal Man nodes
  expect_match_motif_equal(result, expected)
})

# ========== Substituents ==========
test_that("match_motif considers substituents", {
  glycan1 <- parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  glycan2 <- parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(?1-")
  motif1 <- parse_iupac_condensed("Neu5Ac9Ac(a2-3)Gal(b1-")
  motif2 <- parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-")

  # Exact match
  result1 <- match_motif(glycan1, motif1)
  expect_match_motif_equal(result1, list(list(c(1, 2))))  # Neu5Ac=node1, Gal=node2

  # Mismatch due to substituent
  result2 <- match_motif(glycan1, motif2)
  expect_match_motif_equal(result2, list(list()))

  # Match without substituent
  result3 <- match_motif(glycan2, motif2)
  expect_match_motif_equal(result3, list(list(c(1, 2))))  # Neu5Ac=node1, Gal=node2
})

# ========== Edge cases ==========
test_that("match_motif works with single monosaccharide", {
  glycan <- parse_iupac_condensed("Man(?1-")
  motif <- parse_iupac_condensed("Man(?1-")
  result <- match_motif(glycan, motif)
  expected <- list(list(c(1)))
  expect_match_motif_equal(result, expected)
})

test_that("match_motif works with complex branched structure", {
  glycan <- glyrepr::n_glycan_core()
  motif <- parse_iupac_condensed("GlcNAc(b1-4)GlcNAc(?1-")
  result <- match_motif(glycan, motif)
  expected <- list(list(c(4, 5)))
  expect_match_motif_equal(result, expected)
})

# ========== Linkage flexibility ==========
test_that("match_motif works with flexible linkages", {
  glycan <- parse_iupac_condensed("Fuc(a1-6)GlcNAc(b1-")
  motif <- parse_iupac_condensed("Fuc(a1-?)GlcNAc(b1-")  # Flexible position
  result <- match_motif(glycan, motif)
  expected <- list(list(c(1, 2)))
  expect_match_motif_equal(result, expected)
})

test_that("match_motif works with flexible anomer", {
  glycan <- parse_iupac_condensed("Gal(a1-3)GalNAc(b1-")
  motif <- parse_iupac_condensed("Gal(?1-3)GalNAc(b1-")  # Flexible anomer
  result <- match_motif(glycan, motif)
  expected <- list(list(c(1, 2)))
  expect_match_motif_equal(result, expected)
})

# ========== Error handling ==========
test_that("match_motifs handles mismatched alignments length", {
  glycan <- glyrepr::o_glycan_core_2()
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-", "Man(a1-")
  motifs <- parse_iupac_condensed(motifs)
  alignments <- c("substructure", "core")  # length 2, motifs length 3

  expect_error(match_motifs(glycan, motifs, alignments = alignments),
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})

test_that("match_motifs works with single alignment for all motifs", {
  glycan <- glyrepr::o_glycan_core_2()
  motifs <- c("Gal(b1-3)GalNAc(?1-", "GlcNAc(b1-6)GalNAc(?1-")
  motifs <- parse_iupac_condensed(motifs)

  result <- match_motifs(glycan, motifs, alignments = "substructure")
  expected <- list(
    # Motif 1
    list(list(c(1, 3))),  # Gal=node1, GalNAc=node3
    # Motif 2
    list(list(c(2, 3)))   # GlcNAc=node2, GalNAc=node3
  )
  expect_match_motifs_equal(result, expected)
})

# ========== Return value structure ==========
test_that("match_motif returns correct nested list structure", {
  glycan <- glyrepr::n_glycan_core()
  motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
  result <- match_motif(glycan, motif)

  # Check structure: list of glycans, each containing list of matches
  expect_type(result, "list")
  expect_length(result, 1)  # One glycan
  expect_type(result[[1]], "list")
  expect_length(result[[1]], 1)  # One match
  expect_type(result[[1]][[1]], "integer")
  expect_length(result[[1]][[1]], 3)  # Three nodes in motif
})

test_that("match_motifs returns correct nested list structure", {
  glycan <- glyrepr::n_glycan_core()
  motifs <- c("Man(a1-3)[Man(a1-6)]Man(b1-", "GlcNAc(b1-4)GlcNAc(?1-")
  motifs <- parse_iupac_condensed(motifs)
  result <- match_motifs(glycan, motifs)

  # Check structure: list of motifs, each containing list of glycans, each containing list of matches
  expect_type(result, "list")
  expect_length(result, 2)  # Two motifs
  expect_type(result[[1]], "list")
  expect_length(result[[1]], 1)  # One glycan
  expect_type(result[[1]][[1]], "list")
  expect_length(result[[1]][[1]], 1)  # One match for first motif
  expect_type(result[[1]][[1]][[1]], "integer")
})

# ========== Empty results ==========
test_that("match_motif handles empty results correctly", {
  glycan <- parse_iupac_condensed("Man(?1-")
  motif <- parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  result <- match_motif(glycan, motif)

  expect_type(result, "list")
  expect_length(result, 1)  # One glycan
  expect_type(result[[1]], "list")
  expect_length(result[[1]], 0)  # No matches
})

test_that("match_motifs handles mixed empty and non-empty results", {
  glycan <- glyrepr::o_glycan_core_2()
  motifs <- c("Gal(b1-3)GalNAc(?1-", "Neu5Ac(a2-3)Gal(b1-")  # First matches, second doesn't
  motifs <- parse_iupac_condensed(motifs)
  result <- match_motifs(glycan, motifs)

  expect_length(result, 2)  # Two motifs
  expect_length(result[[1]][[1]], 1)  # First motif has one match
  expect_length(result[[2]][[1]], 0)  # Second motif has no matches
})

# ========== Name preservation ==========
test_that("match_motif preserves names from glycan_structure input", {
  glycan1 <- glyrepr::o_glycan_core_1()
  glycan2 <- glyrepr::o_glycan_core_2()
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("core1", "core2")

  result <- match_motif(glycans, glyrepr::o_glycan_core_1())
  expect_equal(names(result), c("core1", "core2"))
})

test_that("match_motif returns no names when input has no names", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())

  result <- match_motif(glycans, glyrepr::o_glycan_core_1())
  expect_null(names(result))
})

# ========== match_degree ==========
test_that("match_motif respects match_degree", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  result_default <- match_motif(glycan, motif)
  expect_true(length(result_default[[1]]) > 0)
  expect_match_motif_equal(match_motif(glycan, motif, match_degree = c(FALSE, TRUE)), list(list()))
})

test_that("match_motifs validates match_degree list", {
  glycans <- glyrepr::o_glycan_core_2()
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"))

  expect_error(match_motifs(glycans, motifs, match_degree = list(c(TRUE, FALSE))), "match_degree.*same length")
})



# ========== match_motifs naming ==========
test_that("match_motifs names outer list by motifs and inner by glycans", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  names(glycans) <- c("core1", "core2")

  motifs <- c(m1 = glyrepr::o_glycan_core_1(), m2 = glyrepr::o_glycan_core_2())

  result <- match_motifs(glycans, motifs)
  expect_equal(names(result), c("m1", "m2"))
  expect_equal(names(result[[1]]), c("core1", "core2"))
})

# ========== Duplicate motifs ==========
test_that("match_motifs raises error for duplicate motifs", {
  glycan <- glyrepr::n_glycan_core()
  glycans <- c(glycan)
  motifs <- glyparse::parse_iupac_condensed(c("Man(b1-", "Man(b1-"))

  expect_error(
    match_motifs(glycans, motifs),
    "cannot have duplications"
  )
})

# ========== Integration tests with motif specs ==========
test_that("match_motifs works with dynamic_motifs()", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))

  result <- match_motifs(glycans, dynamic_motifs(max_size = 2))

  expect_type(result, "list")
  expect_equal(length(result), length(extract_motif(glycans, max_size = 2)))
})

test_that("match_motifs works with branch_motifs()", {
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )

  result <- match_motifs(glycans, branch_motifs())

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

# ========== IUPAC Names for dynamic_motifs and branch_motifs ==========
test_that("match_motifs returns IUPAC strings as outer list names for dynamic_motifs", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))
  result <- match_motifs(glycans, dynamic_motifs(max_size = 2))
  
  # Outermost list names should be the IUPAC strings of extracted motifs
  expect_type(names(result), "character")
  expect_equal(length(names(result)), length(result))
  # All names should be valid IUPAC strings (non-empty)
  expect_true(all(nchar(names(result)) > 0))
})

test_that("match_motifs returns trimmed IUPAC strings as outer list names for branch_motifs", {
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  result <- match_motifs(glycans, branch_motifs())
  
  # Names should be present
  expect_type(names(result), "character")
  expect_equal(length(names(result)), length(result))
  expect_true(all(nchar(names(result)) > 0))
  
  # Names should NOT contain the core suffix
  expect_false(any(grepl(")Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc", names(result), fixed = TRUE)))
  expect_false(any(grepl(")Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc", names(result), fixed = TRUE)))
  
  # Names should end with the branch root linkage pattern
  expect_true(all(grepl("\\([a-z]1-.$", names(result))))
})
