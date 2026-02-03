test_that("dynamic_motifs() creates a dynamic_motifs_spec object", {
  spec <- dynamic_motifs()
  expect_s3_class(spec, "dynamic_motifs_spec")
  expect_equal(spec$max_size, 3)
})

test_that("dynamic_motifs() accepts custom max_size", {
  spec <- dynamic_motifs(max_size = 5)
  expect_equal(spec$max_size, 5)
})

test_that("dynamic_motifs() validates max_size", {
  expect_error(dynamic_motifs(max_size = "invalid"), "Must be")
  expect_error(dynamic_motifs(max_size = 0), ">=")
})

test_that("branch_motifs() creates a branch_motifs_spec object", {
  spec <- branch_motifs()
  expect_s3_class(spec, "branch_motifs_spec")
})

test_that("resolve_motif_spec errors on alignments conflict for dynamic_motifs", {
  spec <- dynamic_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_error(
    resolve_motif_spec(glycans, spec, alignments = "substructure", match_degree = NULL),
    "Cannot specify.*alignments"
  )
})

test_that("resolve_motif_spec errors on match_degree conflict for dynamic_motifs", {
  spec <- dynamic_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_error(
    resolve_motif_spec(glycans, spec, alignments = NULL, match_degree = list(c(TRUE, TRUE))),
    "Cannot specify.*match_degree"
  )
})

test_that("resolve_motif_spec extracts motifs for dynamic_motifs", {
  spec <- dynamic_motifs(max_size = 2)
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))
  result <- resolve_motif_spec(glycans, spec, alignments = NULL, match_degree = NULL)

  expect_type(result, "list")
  expect_named(result, c("motifs", "alignments", "match_degree"))
  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_equal(unique(result$alignments), "substructure")
  expect_null(result$match_degree)
})

test_that("resolve_motif_spec extracts motifs for branch_motifs", {
  spec <- branch_motifs()
  glycans <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
  )
  result <- resolve_motif_spec(glycans, spec, alignments = NULL, match_degree = NULL)

  expect_type(result, "list")
  expect_named(result, c("motifs", "alignments", "match_degree"))
  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_true(length(result$motifs) > 0)
  expect_type(result$match_degree, "list")
  expect_equal(length(result$match_degree), length(result$motifs))
})

test_that("strict_sub is locked for dynamic_motifs", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  expect_error(
    have_motifs(glycans, dynamic_motifs(), strict_sub = FALSE),
    "Cannot specify.*strict_sub"
  )
})

test_that("strict_sub is locked for branch_motifs", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_error(
    have_motifs(glycans, branch_motifs(), strict_sub = FALSE),
    "Cannot specify.*strict_sub"
  )
})

test_that("ignore_linkages is locked for dynamic_motifs", {
  glycans <- c(glyrepr::o_glycan_core_1(), glyrepr::o_glycan_core_2())
  expect_error(
    have_motifs(glycans, dynamic_motifs(), ignore_linkages = TRUE),
    "Cannot specify.*ignore_linkages"
  )
})

test_that("ignore_linkages is locked for branch_motifs", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_error(
    have_motifs(glycans, branch_motifs(), ignore_linkages = TRUE),
    "Cannot specify.*ignore_linkages"
  )
})
