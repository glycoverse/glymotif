suppressPackageStartupMessages(library(glyrepr))
suppressPackageStartupMessages(library(glydb))
pkgload::load_all(".", quiet = TRUE)
out <- "benchmarks/rcpp-floating-matcher"
stopifnot(packageVersion("igraph") == "2.3.3")
Rcpp::sourceCpp(
  file.path(out, "matcher.cpp"),
  cacheDir = file.path(tempdir(), "floating-cache")
)
dictionary <- get("monosaccharides", asNamespace("glyrepr"))
native <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  strict_floating = TRUE,
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
  cpp_floating_match(
    glycans,
    motifs,
    dictionary,
    alignments,
    ignore_linkages,
    strict_sub,
    mode == "lenient",
    match_degree,
    output,
    strict_floating
  )
}
validated <- function(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = "strict",
  strict_floating = TRUE,
  output = "have"
) {
  p <- glymotif:::prepare_motif_args(
    glycans,
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = mode,
    strict_floating = strict_floating
  )
  native(
    p$glycans,
    p$motifs,
    p$alignments,
    p$ignore_linkages,
    p$strict_sub,
    p$match_degree,
    p$mode,
    p$strict_floating,
    output
  )
}
oracle <- function(
  glycans,
  motifs,
  ...,
  strict_floating = TRUE,
  output = "have"
) {
  f <- switch(
    output,
    have = glymotif::have_motifs,
    count = glymotif::count_motifs,
    match = glymotif::match_motifs
  )
  if (output == "match") {
    suppressWarnings(f(glycans, motifs, ...))
  } else {
    suppressWarnings(f(glycans, motifs, ..., strict_floating = strict_floating))
  }
}
