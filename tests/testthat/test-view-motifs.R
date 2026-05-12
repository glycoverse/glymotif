test_that("view_motif returns a glydraw plot for matched residues", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")

  result <- suppressMessages(view_motif(glycan, motif))

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
})

test_that("view_motif accepts structure strings", {
  result <- suppressMessages(view_motif(
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(b1-",
    "Gal(b1-3)GalNAc(b1-"
  ))

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
})

test_that("view_motif returns an unhighlighted plot when no match is found", {
  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  expect_snapshot(result <- view_motif(glycan, motif))

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
})

test_that("view_motif rejects multiple glycans or motifs", {
  glycan <- c(glyrepr::n_glycan_core(), glyrepr::o_glycan_core_2())
  motif <- glyparse::parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")

  expect_error(
    view_motif(glycan, motif),
    "Only one glycan can be visualized at a time"
  )

  glycan <- glyrepr::n_glycan_core()
  motif <- glyparse::parse_iupac_condensed(c(
    "Man(a1-3)[Man(a1-6)]Man(b1-",
    "GlcNAc(b1-4)GlcNAc(?1-"
  ))

  expect_error(
    view_motif(glycan, motif),
    "Only one motif can be visualized at a time"
  )
})
