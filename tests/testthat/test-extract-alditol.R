test_that("extracted alditol status belongs to the original reducing end", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc-ol(a1-")
  result <- extract_motif(glycan, max_size = 2)
  expect_setequal(as.character(result), c(
    "Gal(b1-", "GlcNAc-ol(a1-", "Gal(b1-4)GlcNAc-ol(a1-"
  ))

  graph <- glyrepr::get_structure_graphs(glycan)
  reference <- .extract_motifs_from_graph(
    graph, 2, new.env(hash = TRUE, parent = emptyenv())
  )
  expect_identical(
    as.character(result),
    as.character(glyrepr::as_glycan_structure(reference))
  )
  expect_identical(glyrepr::get_alditol(glycan), TRUE)
})

test_that("deduplication distinguishes ordinary and reduced motifs", {
  input <- c("Glc(a1-", "Glc-ol(a1-")
  for (glycans in list(input, rev(input))) {
    expect_identical(as.character(extract_motif(glycans)), glycans)
    graphs <- as.list(glyrepr::as_glycan_structure(glycans))
    seen <- new.env(hash = TRUE, parent = emptyenv())
    reference <- unlist(lapply(graphs, function(g) {
      .extract_motifs_from_graph(g, 1, seen)
    }), recursive = FALSE)
    expect_identical(
      as.character(glyrepr::as_glycan_structure(reference)), glycans
    )
  }
})

test_that("branch motifs do not inherit a reduced N-glycan core", {
  glycan <- paste0(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc-ol(a1-"
  )
  result <- extract_branch_motif(glycan)
  expect_identical(as.character(result), "Gal(b1-4)GlcNAc(b1-")
  expect_identical(glyrepr::get_alditol(result), FALSE)
  expect_false(any(glyrepr::get_alditol(
    extract_branch_motif(glycan, including_core = TRUE)
  )))
})
