# Run from the glymotif repository root: Rscript benchmarks/rcpp-matcher/run.R
suppressPackageStartupMessages(library(glyrepr))
suppressPackageStartupMessages(library(glydb))
repo <- normalizePath(".")
pkgload::load_all(repo, quiet = TRUE)
out <- "benchmarks/rcpp-matcher"
Rcpp::sourceCpp(
  file.path(out, "matcher.cpp"),
  cacheDir = file.path(tempdir(), "matcher-cache")
)
stopifnot(packageVersion("igraph") == "2.3.3")
dictionary <- get("monosaccharides", asNamespace("glyrepr"))
native <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  output = "have"
) {
  if (is.null(alignments)) {
    alignments <- rep("substructure", length(motifs))
  }
  if (length(alignments) == 1L) {
    alignments <- rep(alignments, length(motifs))
  }
  if (is.null(match_degree)) {
    match_degree <- rep(list(NULL), length(motifs))
  }
  cpp_structure_match(
    glycans,
    motifs,
    dictionary,
    alignments,
    ignore_linkages,
    strict_sub,
    mode == "lenient",
    match_degree,
    output
  )
}
# Keep public argument preparation in the timed path for the conservative result.
validated <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  output = "have"
) {
  p <- glymotif:::prepare_motif_args(
    glycans,
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = mode
  )
  native(
    p$glycans,
    p$motifs,
    p$alignments,
    p$ignore_linkages,
    p$strict_sub,
    p$match_degree,
    p$mode,
    output
  )
}
oracle <- function(glycans, motifs, ..., output = "have") {
  f <- switch(
    output,
    have = glymotif::have_motifs,
    count = glymotif::count_motifs,
    match = glymotif::match_motifs
  )
  suppressWarnings(f(glycans, motifs, ...))
}
