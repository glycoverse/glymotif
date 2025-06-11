test_that("count motifs in glycan", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-")
  expect_equal(count_motif(glycan, motif), 2L)
})


test_that("count motifs with branching", {
  glycan <- glyparse::parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc")
  motif <- glyparse::parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc")
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count symmetrical motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal")
  expect_equal(count_motif(glycan, motif), 1L)
})


test_that("count 0 motif", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  motif <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc")
  expect_equal(count_motif(glycan, motif), 0L)
})


# ========== count_motifs ==========
test_that("count_motifs works with multiple motifs", {
  glycan1 <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("double_gal", "single_gal")
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("glycan", "Gal(b1-3)GalNAc", "Gal(b1-"))
  expect_equal(result$glycan, c("double_gal", "single_gal"))
  
  # Check specific counts
  expect_equal(unname(result[1, "Gal(b1-3)GalNAc"][[1]]), 1L)
  expect_equal(unname(result[1, "Gal(b1-"][[1]]), 2L)
  expect_equal(unname(result[2, "Gal(b1-3)GalNAc"][[1]]), 1L)
  expect_equal(unname(result[2, "Gal(b1-"][[1]]), 1L)
})


test_that("count_motifs works without glycan names", {
  glycan1 <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("glycan", "Gal(b1-3)GalNAc", "Gal(b1-"))
  
  # Check that glycan column contains IUPAC strings
  expect_true(grepl("Gal\\(b1-3\\)", result$glycan[1]))
  expect_true(grepl("Gal\\(b1-3\\)", result$glycan[2]))
})


test_that("count_motifs works with complex branching motifs", {
  glycan1 <- glyparse::parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc")
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal")
  glycans <- c(glycan1, glycan2)
  names(glycans) <- c("complex", "simple")
  
  motifs <- c("Man(b1-?)[Man(b1-?)]GalNAc", "Man(b1-", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
  
  # Check complex structure counts
  expect_equal(unname(result[1, "Man(b1-?)[Man(b1-?)]GalNAc"][[1]]), 1L)
  expect_equal(unname(result[1, "Man(b1-"][[1]]), 2L)
  expect_equal(unname(result[1, "Gal(b1-"][[1]]), 0L)
  
  # Check simple structure counts
  expect_equal(unname(result[2, "Man(b1-?)[Man(b1-?)]GalNAc"][[1]]), 0L)
  expect_equal(unname(result[2, "Man(b1-"][[1]]), 0L)
  expect_equal(unname(result[2, "Gal(b1-"][[1]]), 1L)
})


test_that("count_motifs works with different alignments", {
  glycan <- glyparse::parse_iupac_condensed("Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal")
  glycans <- c(glycan, glycan)
  names(glycans) <- c("glycan1", "glycan2")
  
  motifs <- c("Gal(a1-3)Gal(a1-4)Gal", "Gal(a1-4)Gal(a1-6)Gal")
  alignments <- c("core", "terminal")
  
  result <- count_motifs(glycans, motifs, alignments = alignments)
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})


test_that("count_motifs handles zero counts", {
  glycan1 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycan2 <- glyparse::parse_iupac_condensed("GlcNAc(b1-4)GlcNAc")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Man(b1-", "Fuc(a1-")
  
  result <- count_motifs(glycans, motifs)
  
  expect_s3_class(result, "tbl_df")
  expect_true(all(result[["Man(b1-"]] == 0L))
  expect_true(all(result[["Fuc(a1-"]] == 0L))
})


test_that("count_motifs handles empty motifs", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan)
  motifs <- character(0)
  
  expect_error(count_motifs(glycans, motifs), "`motifs` cannot be empty")
})


test_that("count_motifs handles invalid motifs argument", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan)
  motifs <- 123
  
  expect_error(count_motifs(glycans, motifs), "`motifs` must be either")
})


test_that("count_motifs handles mismatched alignments length", {
  glycan <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan)
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-", "GalNAc")
  alignments <- c("substructure", "core")  # length 2, motifs length 3
  
  expect_error(count_motifs(glycans, motifs, alignments = alignments), 
               "`alignments` must be NULL, a single value, or have the same length as `motifs`")
})


test_that("count_motifs works with single alignment for all motifs", {
  glycan1 <- glyparse::parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc")
  glycan2 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan1, glycan2)
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-")
  
  result <- count_motifs(glycans, motifs, alignments = "substructure")
  
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})


test_that("count_motifs works with ignore_linkages", {
  glycan1 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc")
  glycans <- c(glycan1)
  
  motifs <- c("Gal(b1-3)GalNAc", "Gal(b1-4)GalNAc")
  
  result_normal <- count_motifs(glycans, motifs, ignore_linkages = FALSE)
  result_ignore <- count_motifs(glycans, motifs, ignore_linkages = TRUE)
  
  expect_equal(unname(result_normal[1, "Gal(b1-4)GalNAc"][[1]]), 0L)
  expect_equal(unname(result_ignore[1, "Gal(b1-4)GalNAc"][[1]]), 1L)
})
