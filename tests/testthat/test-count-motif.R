test_that("count motifs in glycan", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-"
  expect_equal(count_motif(glycan, motif), 2L)
})


test_that("count motifs with branching", {
  glycan <- "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-"
  motif <- "Man(b1-?)[Man(b1-?)]GalNAc(b1-"
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count symmetrical motif", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-3)Gal(b1-"
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count 0 motif", {
  glycan <- "Gal(b1-3)Gal(b1-3)GalNAc(b1-"
  motif <- "Gal(b1-4)GalNAc(b1-"
  expect_equal(count_motif(glycan, motif), 0L)
})


# ========== count_motifs ==========
test_that("count_motifs works with multiple motifs", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  names(glycans) <- c("double_gal", "single_gal")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Gal(b1-3)GalNAc(b1-", "Gal(b1-"))
  expect_equal(rownames(result), c("double_gal", "single_gal"))
  
  # Check specific counts
  expect_equal(result[1, "Gal(b1-3)GalNAc(b1-"], 1L)
  expect_equal(result[1, "Gal(b1-"], 2L)
  expect_equal(result[2, "Gal(b1-3)GalNAc(b1-"], 1L)
  expect_equal(result[2, "Gal(b1-"], 1L)
})


test_that("count_motifs works without glycan names", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("Gal(b1-3)GalNAc(b1-", "Gal(b1-"))
  
  # Check that row names contain IUPAC strings
  expect_true(grepl("Gal\\(b1-3\\)", rownames(result)[1]))
  expect_true(grepl("Gal\\(b1-3\\)", rownames(result)[2]))
})


test_that("count_motifs works with complex branching motifs", {
  glycans <- c("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-", "Gal(b1-3)Gal(b1-")
  names(glycans) <- c("complex", "simple")
  
  motifs <- c("Man(b1-?)[Man(b1-?)]GalNAc(b1-", "Man(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  
  # Check complex structure counts
  expect_equal(result[1, "Man(b1-?)[Man(b1-?)]GalNAc(b1-"], 1L)
  expect_equal(result[1, "Man(b1-"], 2L)
  expect_equal(result[1, "Gal(b1-"], 0L)
  
  # Check simple structure counts
  expect_equal(result[2, "Man(b1-?)[Man(b1-?)]GalNAc(b1-"], 0L)
  expect_equal(result[2, "Man(b1-"], 0L)
  expect_equal(result[2, "Gal(b1-"], 2L)
})


test_that("count_motifs works with different alignments", {
  glycans <- c("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(a1-", "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(a1-")
  names(glycans) <- c("glycan1", "glycan2")
  
  motifs <- c("Gal(a1-3)Gal(a1-4)Gal(a1-", "Gal(a1-4)Gal(a1-6)Gal(a1-")
  alignments <- c("core", "terminal")
  
  result <- count_motifs(glycans, motifs, alignments = alignments)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("count_motifs handles zero counts", {
  glycans <- c("Gal(b1-3)GalNAc(b1-", "GlcNAc(b1-4)GlcNAc(b1-")
  
  motifs <- c("Man(b1-", "Fuc(a1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_true(is.matrix(result))
  expect_true(all(result[, "Man(b1-"] == 0L))
  expect_true(all(result[, "Fuc(a1-"] == 0L))
})


test_that("count_motifs handles empty motifs", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- character(0)
  
  expect_error(count_motifs(glycans, motifs), "`motifs` cannot be empty")
})


test_that("count_motifs handles invalid motifs argument", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- 123

  expect_error(count_motifs(glycans, motifs), "`motifs` must be a 'glyrepr_structure' object")
})


test_that("count_motifs handles mismatched alignments length", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-", "GalNAc(b1-")
  alignments <- c("substructure", "core")  # length 2, motifs length 3
  
  expect_error(count_motifs(glycans, motifs, alignments = alignments), 
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})


test_that("count_motifs works with single alignment for all motifs", {
  glycans <- c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-")
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs, alignments = "substructure")
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})


test_that("count_motifs works with ignore_linkages", {
  glycans <- "Gal(b1-3)GalNAc(b1-"
  
  motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-4)GalNAc(b1-")
  
  result_normal <- count_motifs(glycans, motifs, ignore_linkages = FALSE)
  result_ignore <- count_motifs(glycans, motifs, ignore_linkages = TRUE)
  
  expect_equal(result_normal[1, "Gal(b1-4)GalNAc(b1-"], 0L)
  expect_equal(result_ignore[1, "Gal(b1-4)GalNAc(b1-"], 1L)
})
