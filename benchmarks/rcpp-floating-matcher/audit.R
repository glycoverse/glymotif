source("benchmarks/rcpp-floating-matcher/setup.R")
all <- glydb_structures()
ids <- which(has_floating_parts(all) | has_floating_substituents(all))
g <- all[ids]
graphs <- attr(g, "graphs")[as.character(g)]
error_text <- function(f) {
  tryCatch(
    {
      f()
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}
rows <- list()
for (i in seq_along(g)) {
  cpp <- tryCatch(cpp_localizations(graphs[[i]]), error = identity)
  ref <- tryCatch(
    glyrepr::enumerate_floating_graph_localizations(graphs[[i]]),
    error = identity
  )
  if (inherits(ref, "error")) {
    stopifnot(
      inherits(cpp, "error"),
      grepl("max_variants", conditionMessage(ref)),
      grepl("max_variants", conditionMessage(cpp))
    )
    rows[[i]] <- data.frame(
      id = ids[i],
      raw = NA_real_,
      variants = NA_integer_,
      result = "same_limit_error"
    )
  } else {
    stopifnot(!inherits(cpp, "error"), length(cpp$variants) == nrow(ref))
    for (j in seq_len(nrow(ref))) {
      rg <- ref$graph[[j]]
      cg <- cpp$variants[[j]]
      expected <- list(
        edges = matrix(
          as.integer(igraph::as_edgelist(rg, names = FALSE)),
          ncol = 2
        ),
        linkage = unname(igraph::edge_attr(rg, "linkage")),
        sub = unname(igraph::vertex_attr(rg, "sub")),
        parents = ref$assignments[[j]]$parent_node
      )
      if (!identical(cg, expected)) {
        saveRDS(
          list(i = i, j = j, cpp = cg, expected = expected),
          file.path(out, "localization-failure.rds")
        )
        stop("Localization mismatch at ", i, "/", j)
      }
    }
    rows[[i]] <- data.frame(
      id = ids[i],
      raw = cpp$raw,
      variants = length(cpp$variants),
      result = "identical"
    )
  }
  write.csv(
    do.call(rbind, rows),
    file.path(out, "localization-audit.csv"),
    row.names = FALSE
  )
  if (i %% 10L == 0L) cat("Localizations", i, "/", length(g), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
saveRDS(
  list(corpus_ids = ids, keys = as.character(g)),
  file.path(out, "inputs.rds")
)
cat("LOCALIZATION AUDIT DONE\n")
