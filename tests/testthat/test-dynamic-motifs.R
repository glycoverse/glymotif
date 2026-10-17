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
  expect_identical(spec$include_empty_branch, FALSE)
  expect_identical(spec$include_mannose_branch, FALSE)
})

test_that("branch_motifs counts special branches without changing defaults", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-3)[GlcNAc(b1-4)][Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
  ))
  baseline <- count_motifs(glycans, branch_motifs())
  expect_identical(colnames(baseline), "GlcNAc(b1-")
  expect_identical(unname(baseline[, 1]), c(0L, 0L, 1L, 2L, 0L))

  spec <- branch_motifs(
    include_empty_branch = TRUE,
    include_mannose_branch = TRUE
  )
  counts <- count_motifs(glycans, spec)
  expect_setequal(colnames(counts), c("GlcNAc(b1-", "Man(a1-", "empty"))
  expect_identical(unname(counts[, "empty"]), c(2L, 1L, 0L, 0L, 2L))
  expect_identical(unname(counts[, "Man(a1-"]), c(0L, 1L, 1L, 0L, 0L))
  expect_identical(counts[, "GlcNAc(b1-"], baseline[, 1])
  expect_identical(have_motifs(glycans, spec), counts > 0L)
  expect_setequal(
    colnames(count_motifs(glycans, branch_motifs(include_empty_branch = TRUE))),
    c("GlcNAc(b1-", "empty")
  )
  expect_setequal(
    colnames(count_motifs(
      glycans,
      branch_motifs(include_mannose_branch = TRUE)
    )),
    c("GlcNAc(b1-", "Man(a1-")
  )
})

test_that("branch_motifs matches full mannose subtrees and returns arm indices", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  spec <- branch_motifs(
    include_empty_branch = TRUE,
    include_mannose_branch = TRUE
  )
  counts <- count_motifs(glycans, spec)
  expect_identical(unname(counts[, "Man(a1-2)Man(a1-"]), c(1L, 0L))
  expect_identical(unname(counts[, "Man(a1-"]), c(0L, 1L))

  matches <- match_motifs(glycans, spec)
  expect_equal(
    matches[["Man(a1-2)Man(a1-"]],
    list(list(c(1, 2)), list())
  )
  expect_equal(matches[["Man(a1-"]], list(list(), list(1)))
  expect_equal(matches$empty, list(list(4), list(3)))
})

test_that("special branch matching supports generic and mixed topologies", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Hex(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "HexNAc(??-?)Man(??-?)[HexNAc(??-?)Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  ))
  spec <- branch_motifs(
    include_empty_branch = TRUE,
    include_mannose_branch = TRUE
  )
  counts <- count_motifs(glycans, spec)
  expect_identical(unname(counts[, "empty (generic)"]), c(2L, 1L, 0L))
  expect_identical(unname(counts[, "Hex(??-"]), c(0L, 1L, 0L))
  expect_identical(unname(counts[, "HexNAc(??-"]), c(0L, 0L, 2L))
  expect_identical(have_motifs(glycans, spec), counts > 0L)
})

test_that("special branch options require logical scalars", {
  glycan <- glyrepr::n_glycan_core()
  expect_snapshot(error = TRUE, branch_motifs(include_empty_branch = NA))
  expect_snapshot(error = TRUE, branch_motifs(include_mannose_branch = 1))
  expect_snapshot(
    error = TRUE,
    extract_branch_motif(glycan, include_empty_branch = NA)
  )
  expect_snapshot(
    error = TRUE,
    extract_branch_motif(glycan, include_mannose_branch = 1)
  )
})

test_that("mannose motifs match concrete subtrees on generic and mixed cores", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Man(??-?)Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "Man(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Man(??-?)Man(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-"
  ))
  spec <- branch_motifs(include_mannose_branch = TRUE)
  counts <- count_motifs(glycans, spec)
  expect_identical(anyDuplicated(colnames(counts)), 0L)
  expect_equal(
    unname(counts),
    rbind(c(1L, 1L, 1L), c(0L, 1L, 0L), c(0L, 1L, 1L))
  )
  expect_identical(have_motifs(glycans, spec), counts > 0L)
  matches <- match_motifs(glycans, spec)
  for (i in seq_along(glycans)) {
    expect_equal(matches[[i]][[i]], list(1L))
  }
  expect_identical(
    unname(count_motifs(
      glyrepr::n_glycan_core(),
      branch_motifs(include_empty_branch = TRUE, include_mannose_branch = TRUE)
    )),
    matrix(2L, nrow = 1)
  )
})

test_that("empty branch matching counts substituted arms individually", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Man6S(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-3)[Man6S(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man6S(a1-3)[Man6S(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  spec <- branch_motifs(include_empty_branch = TRUE)
  counts <- count_motifs(glycans, spec)
  expect_setequal(colnames(counts), c("empty", "empty (Man6S)"))
  expect_identical(unname(counts[, "empty"]), c(1L, 1L, 0L))
  expect_identical(unname(counts[, "empty (Man6S)"]), c(1L, 1L, 2L))
  expect_identical(have_motifs(glycans, spec), counts > 0L)
  matches <- match_motifs(glycans, spec)
  expect_equal(
    matches[["empty (Man6S)"]],
    list(list(1L), list(2L), list(1L, 2L))
  )
})

test_that("substituted empty core paths support generic and mixed identities", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Hex6S(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Man6S(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc6S(b1-4)GlcNAc(b1-"
  ))
  for (i in seq_along(glycans)) {
    glycan <- glycans[i]
    spec <- branch_motifs(include_empty_branch = TRUE)
    counts <- count_motifs(glycan, spec)
    expect_identical(sum(counts), 2L)
    expect_identical(have_motifs(glycan, spec), counts > 0L)
    expect_setequal(
      unlist(match_motifs(glycan, spec), use.names = FALSE),
      c(1L, 2L)
    )
  }
})

test_that("mannose branches retain substituents on generic source cores", {
  glycan <- glyrepr::as_glycan_structure(
    "Man6S(a1-2)Hex6S(??-?)[Hex(??-?)]Hex(??-?)HexNAc6S(??-?)HexNAc(??-"
  )
  spec <- branch_motifs(include_mannose_branch = TRUE)
  counts <- count_motifs(glycan, spec)
  expect_identical(colnames(counts), "Man6S(a1-")
  expect_identical(unname(counts), matrix(1L, nrow = 1))
  expect_identical(have_motifs(glycan, spec), counts > 0L)
  expect_equal(match_motifs(glycan, spec)[[1]][[1]], list(1L))
})

test_that("db_motifs() creates a db_motifs_spec object", {
  spec <- db_motifs()
  expect_s3_class(spec, "db_motifs_spec")
})

test_that("resolve_motif_spec errors on alignments conflict for dynamic_motifs", {
  spec <- dynamic_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_error(
    resolve_motif_spec(
      glycans,
      spec,
      alignments = "substructure",
      match_degree = NULL
    ),
    "Cannot specify.*alignments"
  )
})

test_that("resolve_motif_spec errors on match_degree conflict for dynamic_motifs", {
  spec <- dynamic_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_error(
    resolve_motif_spec(
      glycans,
      spec,
      alignments = NULL,
      match_degree = list(c(TRUE, TRUE))
    ),
    "Cannot specify.*match_degree"
  )
})

test_that("resolve_motif_spec extracts motifs for dynamic_motifs", {
  spec <- dynamic_motifs(max_size = 2)
  glycans <- glyrepr::as_glycan_structure(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Man(b1-4)GlcNAc(b1-"
  ))
  result <- resolve_motif_spec(
    glycans,
    spec,
    alignments = NULL,
    match_degree = NULL
  )

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
  result <- resolve_motif_spec(
    glycans,
    spec,
    alignments = NULL,
    match_degree = NULL
  )

  expect_type(result, "list")
  expect_named(result, c("motifs", "alignments", "match_degree"))
  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_true(length(result$motifs) > 0)
  expect_type(result$match_degree, "list")
  expect_equal(length(result$match_degree), length(result$motifs))
})

test_that("resolve_motif_spec extracts motifs for db_motifs", {
  spec <- db_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  result <- resolve_motif_spec(
    glycans,
    spec,
    alignments = NULL,
    match_degree = NULL
  )
  info <- db_motif_info()
  info <- info[info$source_id == "GGM", ]

  expect_type(result, "list")
  expect_named(
    result,
    c("motifs", "alignments", "match_degree", "allow_duplicate_motifs")
  )
  expect_s3_class(result$motifs, "glyrepr_structure")
  expect_identical(
    unname(as.character(result$motifs)),
    as.character(info$glycan_structure)
  )
  expect_identical(result$alignments, info$alignment)
  expect_null(result$match_degree)
  expect_true(result$allow_duplicate_motifs)
  expect_gt(sum(duplicated(as.character(result$motifs))), 0)
})

test_that("resolve_motif_spec errors on alignments conflict for db_motifs", {
  spec <- db_motifs()
  glycans <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-")
  expect_error(
    resolve_motif_spec(
      glycans,
      spec,
      alignments = "substructure",
      match_degree = NULL
    ),
    "Cannot specify.*alignments"
  )
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
