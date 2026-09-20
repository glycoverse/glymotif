# Additional full-corpus parity, unsupported-input guards and scalar positive timing.
source("benchmarks/rcpp-matcher/setup.R")
source(file.path(out, "timing.R"))
old_timings <- read.csv(file.path(out, "timings.csv"))
rows <- lapply(seq_len(nrow(old_timings)), function(i) {
  old_timings[i, , drop = FALSE]
})
saved <- readRDS(file.path(out, "inputs.rds"))
all_g <- glydb_structures()
floating <- has_floating_parts(all_g) | has_floating_substituents(all_g)
g <- all_g[!is.na(all_g) & !floating]
m <- as_glycan_structure(saved$motif_keys[1:13])
stopifnot(identical(as.character(g), saved$corpus_keys))
extra <- list()
for (kind in c("have", "count", "match")) {
  expected <- oracle(g, m, output = kind)
  actual <- native(g, m, output = kind)
  stopifnot(identical(actual, expected))
  extra[[length(extra) + 1L]] <- data.frame(
    case = "full-corpus",
    output = kind,
    pairs = length(g) * length(m),
    identical = TRUE
  )
  cat("Full corpus identical:", kind, "\n")
}
# All excluded structures must fail explicitly, not return plausible wrong results.
rejected <- vapply(
  which(floating),
  function(i) {
    err <- tryCatch(
      {
        native(all_g[i], m[1])
        NULL
      },
      error = identity
    )
    inherits(err, "error") &&
      grepl("Floating localization is not implemented", conditionMessage(err))
  },
  logical(1)
)
stopifnot(all(rejected))
cat("Floating structures explicitly rejected:", sum(rejected), "\n")
write.csv(
  rbind(read.csv(file.path(out, "parity.csv")), do.call(rbind, extra)),
  file.path(out, "parity.csv"),
  row.names = FALSE
)
stopifnot(isTRUE(unname(oracle(m[1], m[1])[[1]])))
benchmark("scalar-positive", m[1], m[1], iterations = 100L)
raw <- do.call(rbind, rows)
summary <- aggregate(
  cbind(elapsed_s, cpu_s) ~ workload + output + engine,
  raw,
  median
)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(
  paste("floating_rejected", sum(rejected)),
  file.path(out, "guards.txt")
)
unlink(file.path(out, "failure.rds"))
cat("SUPPLEMENT DONE\n")
