source("benchmarks/rcpp-floating-matcher/setup.R")
audit <- read.csv(file.path(out, "localization-audit.csv"))
g <- glydb_structures()[audit$id[audit$result == "identical"]]
m <- as_glycan_structure(c(
  "Gal(??-?)GlcNAc(??-",
  "Neu5Ac(??-?)Gal(??-",
  "Fuc(??-?)GlcNAc(??-",
  "Hex(??-?)HexNAc(??-",
  "Gal6S(??-",
  "Man(??-?)Man(??-"
))
ng <- length(g)
nm <- length(m)
graphs <- attr(g, "graphs")[as.character(g)]
mg <- attr(m, "graphs")[as.character(m)]
profiles <- lapply(mg, glymotif:::new_motif_composition_profile)
linked <- vapply(mg, glymotif:::graph_has_linkages, logical(1))
low <- high <- matrix(
  NA_integer_,
  ng,
  nm,
  dimnames = list(as.character(g), NULL)
)
maps <- replicate(nm, vector("list", ng), simplify = FALSE)
# Cache R localization/mapping results only for the correctness oracle. Timed
# public API calls in benchmark.R recompute everything from input objects.
for (i in seq_len(ng)) {
  variants <- glymotif:::floating_localization_graphs(graphs[[i]])
  for (j in seq_len(nm)) {
    result <- lapply(variants, function(graph) {
      glymotif:::.match_motif_single(
        graph,
        mg[[j]],
        linked[[j]],
        profiles[[j]],
        alignment = "substructure"
      )
    })
    low[i, j] <- min(lengths(result))
    high[i, j] <- max(lengths(result))
    maps[[j]][[i]] <- glymotif:::unique_floating_matches(result)
  }
  if (i %% 10L == 0L) cat("Reference matching", i, "/", ng, "\n")
}
expected <- list(
  have_true = low > 0L,
  have_false = high > 0L,
  count_true = low,
  count_false = high,
  match = maps
)
for (sf in c(TRUE, FALSE)) {
  for (kind in c("have", "count", "match")) {
    ref <- if (kind == "match") {
      expected$match
    } else {
      expected[[paste(kind, tolower(sf), sep = "_")]]
    }
    actual <- native(g, m, strict_floating = sf, output = kind)
    if (!identical(actual, ref)) {
      saveRDS(
        list(kind = kind, sf = sf, actual = actual, ref = ref),
        file.path(out, "corpus-failure.rds")
      )
      stop(kind, "/", sf, ": ", paste(all.equal(actual, ref), collapse = "; "))
    }
    cat("Full corpus identical", kind, sf, "\n")
  }
}
# Cross-check the assembled reference against the public API on a fixed subset.
for (sf in c(TRUE, FALSE)) {
  for (kind in c("have", "count", "match")) {
    stopifnot(identical(
      oracle(g[1:12], m, strict_floating = sf, output = kind),
      native(g[1:12], m, strict_floating = sf, output = kind)
    ))
  }
}
saveRDS(
  list(g = g, m = m, expected = expected),
  file.path(out, "corpus-results.rds")
)
write.csv(
  data.frame(
    output = c("have", "count", "match"),
    pairs = ng * nm,
    floating_modes = 2L,
    identical = TRUE
  ),
  file.path(out, "corpus-parity.csv"),
  row.names = FALSE
)
cat("CORPUS MATCH AUDIT DONE\n")
