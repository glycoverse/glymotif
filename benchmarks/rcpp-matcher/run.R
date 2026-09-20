source("benchmarks/rcpp-matcher/setup.R")
checks <- list()
verify <- function(label, g, m, ...) {
  for (kind in c("have", "count", "match")) {
    expected <- oracle(g, m, ..., output = kind)
    actual <- native(g, m, ..., output = kind)
    if (!identical(actual, expected)) {
      saveRDS(
        list(
          label = label,
          g = g,
          m = m,
          args = list(...),
          kind = kind,
          expected = expected,
          actual = actual
        ),
        file.path(out, "failure.rds")
      )
      stop(
        label,
        " / ",
        kind,
        ": ",
        paste(all.equal(actual, expected), collapse = "; ")
      )
    }
    checks[[length(checks) + 1L]] <<- data.frame(
      case = label,
      output = kind,
      pairs = length(g) * length(m),
      identical = TRUE
    )
  }
}
fixture_text <- c(
  "Gal(b1-4)GlcNAc",
  "Gal(b1-3)GlcNAc",
  "Gal(b1-3/4)GlcNAc",
  "Gal(?1-?)GlcNAc",
  "Gal(??-?)GlcNAc(??-",
  "Hex(??-?)HexNAc(??-",
  "Hex(b1-4)GlcNAc",
  "Gal",
  "Man",
  "Hex",
  "HexNAc",
  "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
  "Hex(?1-?)[Hex(?1-?)]Hex(?1-",
  "Glc3Me6S(a1-",
  "Glc?Me6S(a1-",
  "Glc3Me?S(a1-",
  "Glc?Me?Me(a1-",
  "Glc3Me6Me(a1-",
  "Glc(a1-",
  "GlcNAc",
  "Glc?NAc",
  "Neu5Ac",
  "Neu?Ac",
  "Gal(a1-3)Gal(a1-4)Gal(a1-6)Gal(?1-",
  "Galf(b1-4)GlcfNAc",
  "D-Fuc",
  "Glc-ol"
)
fixtures <- unique(as_glycan_structure(fixture_text))
names(fixtures) <- paste0("fixture", seq_along(fixtures))
fg <- c(fixtures, fixtures[c(1, 1, 3)], as_glycan_structure(NA_character_))
cat("Fixture structures:", length(fixtures), "\n")
for (mode in c("strict", "lenient")) {
  for (alignment in c("substructure", "core", "terminal", "whole")) {
    for (strict in c(TRUE, FALSE)) {
      for (ignore in c(TRUE, FALSE)) {
        label <- paste(mode, alignment, strict, ignore, sep = "/")
        verify(
          label,
          fg,
          fixtures,
          mode = mode,
          alignments = rep(alignment, length(fixtures)),
          strict_sub = strict,
          ignore_linkages = ignore
        )
      }
    }
  }
}
cat("Fixture parameter grid passed\n")
for (degree in c(FALSE, TRUE)) {
  masks <- lapply(attr(fixtures, "graphs"), function(g) {
    rep(degree, igraph::vcount(g))
  })
  masks <- unname(masks[match(
    as.character(fixtures),
    names(attr(fixtures, "graphs"))
  )])
  # The masks above deliberately follow the structure keys, not dictionary order.
  verify(paste0("degree-", degree), fg, fixtures, match_degree = masks)
}
verify("empty-glycans", fixtures[integer()], fixtures)
all_g <- glydb::glydb_structures()
floating <- has_floating_parts(all_g) | has_floating_substituents(all_g)
ordinary <- all_g[!is.na(all_g) & !floating]
set.seed(20260920)
sample_ids <- sample(seq_along(ordinary), min(512L, length(ordinary)))
g <- ordinary[sample_ids]
m <- unique(c(fixtures[1:13], g[1:12]))
# Strip names only from motifs to check both output naming modes.
names(m) <- NULL
cat(
  "Corpus:",
  length(all_g),
  "ordinary:",
  length(ordinary),
  "floating:",
  sum(floating),
  "\n"
)
for (mode in c("strict", "lenient")) {
  for (alignment in c("substructure", "core", "terminal", "whole")) {
    verify(
      paste("corpus", mode, alignment),
      g,
      m,
      mode = mode,
      alignments = rep(alignment, length(m))
    )
    cat("Corpus parity:", mode, alignment, "\n")
  }
}
write.csv(
  do.call(rbind, checks),
  file.path(out, "parity.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    seed = 20260920,
    sample_ids = sample_ids,
    glycan_keys = as.character(g),
    motif_keys = as.character(m),
    corpus_keys = as.character(ordinary)
  ),
  file.path(out, "inputs.rds")
)
cat(
  "Parity complete:",
  sum(vapply(checks, function(x) x$pairs, numeric(1))),
  "pair-output checks\n"
)

# Each measurement starts with objects and includes extraction, preparation,
# C++ graph construction, matching, restoration and returned R result allocation.
# Parsing input strings, compilation, and static residue dictionary setup excluded.
source(file.path(out, "timing.R"))
benchmark("scalar", g[1], m[1], iterations = 100L)
benchmark("512x25", g, m)
benchmark("512x25", g, m, kind = "count")
benchmark("512x25", g, m, kind = "match")
benchmark("duplicates-4096x25", rep(g[1:64], 64), m)
benchmark("full-corpus", ordinary, m[1:13], reps = 5L)
raw <- do.call(rbind, rows)
summary <- aggregate(
  cbind(elapsed_s, cpu_s) ~ workload + output + engine,
  raw,
  median
)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("DONE\n")
