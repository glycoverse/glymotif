source("benchmarks/rcpp-floating-matcher/setup.R")
text <- c(
  "{Gal(a1-3)|2,3}Man(a1-3)[Glc(a1-6)]GlcNAc(b1-",
  "{6S}Gal(a1-3)Glc(a1-",
  "{Neu5Ac(a2-3)|2,3}Gal(??-?)[Gal(??-?)]GlcNAc(??-",
  "{Gal(a1-3)|3,4}{Gal(a1-4)|3,4}Man(a1-3)[Man(a1-6)]GlcNAc(b1-",
  "{6S|2,3}{Fuc(a1-6)|2,3}Gal(a1-3)Glc(a1-",
  "{Gal(??-?)GlcNAc(??-?)}Gal(??-?)GlcNAc(??-?)Gal(??-",
  "{Gal(a1-3)|3,4}{Gal(a1-3)|3,4}Man(a1-3)[Man(a1-6)]GlcNAc(b1-",
  "{Man(a1-?)}{Man(a1-?)}Man(a1-",
  "{Neu5Ac(a2-3/6)|2,3}Gal(b1-4)GalNAc(a1-",
  "{?S|1,2}Gal(a1-3)Glc(a1-",
  "{3/6S|1,2}Gal(a1-3)Glc(a1-",
  "{6S}{6S}Gal(a1-3)Glc(a1-",
  "Gal6S(a1-3)Glc(a1-"
)
g <- as_glycan_structure(text)
g <- c(g, g[1], as_glycan_structure(NA_character_))
names(g) <- paste0("g", seq_along(g))
m <- as_glycan_structure(c(
  "Gal(a1-3)Man(a1-",
  "Gal6S(a1-",
  "Neu5Ac(??-?)Gal(??-",
  "Gal(a1-3)[Gal(a1-4)]Man(a1-",
  "Gal(??-?)GlcNAc(??-?)Gal(??-",
  "Gal(??-?)GlcNAc(??-",
  "Man(a1-?)Man(a1-",
  "Hex(??-?)Hex(??-",
  "Gal?S(a1-",
  "Gal3/6S(a1-",
  "Glc(a1-",
  "Fuc(a1-6)Gal(a1-"
))
names(m) <- paste0("m", seq_along(m))
rows <- list()
for (i in seq_len(length(text))) {
  graph <- attr(g, "graphs")[[as.character(g[i])]]
  a <- cpp_localizations(graph)
  b <- glyrepr::enumerate_floating_graph_localizations(graph)
  stopifnot(length(a$variants) == nrow(b))
  for (j in seq_len(nrow(b))) {
    bg <- b$graph[[j]]
    expected <- list(
      edges = matrix(
        as.integer(igraph::as_edgelist(bg, names = FALSE)),
        ncol = 2
      ),
      linkage = unname(igraph::edge_attr(bg, "linkage")),
      sub = unname(igraph::vertex_attr(bg, "sub")),
      parents = b$assignments[[j]]$parent_node
    )
    stopifnot(identical(a$variants[[j]], expected))
  }
}
cat("Fixture localization graphs identical\n")
verify <- function(label, ...) {
  for (output in c("have", "count", "match")) {
    actual <- native(g, m, ..., output = output)
    expected <- oracle(g, m, ..., output = output)
    if (!identical(actual, expected)) {
      saveRDS(
        list(
          label = label,
          output = output,
          actual = actual,
          expected = expected
        ),
        file.path(out, "fixture-failure.rds")
      )
      stop(
        label,
        " / ",
        output,
        ": ",
        paste(all.equal(actual, expected), collapse = "; ")
      )
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      case = label,
      output = output,
      pairs = length(g) * length(m),
      identical = TRUE
    )
  }
  write.csv(
    do.call(rbind, rows),
    file.path(out, "fixture-parity.csv"),
    row.names = FALSE
  )
  cat(label, "passed\n")
}
for (mode in c("strict", "lenient")) {
  for (alignment in c("substructure", "core", "terminal", "whole")) {
    for (sf in c(TRUE, FALSE)) {
      verify(
        paste(mode, alignment, sf),
        mode = mode,
        alignments = alignment,
        strict_floating = sf
      )
    }
  }
}
for (ignore in c(TRUE, FALSE)) {
  for (ss in c(TRUE, FALSE)) {
    verify(
      paste("options", ignore, ss),
      ignore_linkages = ignore,
      strict_sub = ss,
      strict_floating = FALSE
    )
  }
}
for (degree in c(TRUE, FALSE)) {
  masks <- lapply(as.character(m), function(key) {
    rep(degree, igraph::vcount(attr(m, "graphs")[[key]]))
  })
  verify(paste("degree", degree), match_degree = masks, strict_floating = FALSE)
}
# Same public rejection of unresolved floating motifs.
stopifnot(
  inherits(tryCatch(native(g, g[1]), error = identity), "error"),
  inherits(tryCatch(oracle(g, g[1]), error = identity), "error")
)
saveRDS(list(g = g, m = m), file.path(out, "fixture-inputs.rds"))
cat("FIXTURES DONE\n")
