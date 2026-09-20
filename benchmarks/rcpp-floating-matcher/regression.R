source("benchmarks/rcpp-floating-matcher/setup.R")
inputs <- readRDS("benchmarks/rcpp-matcher/inputs.rds")
g <- as_glycan_structure(inputs$glycan_keys[1:128])
m <- as_glycan_structure(inputs$motif_keys)
rows <- list()
for (kind in c("have", "count", "match")) {
  before <- serialize(list(g, m), NULL)
  actual <- native(g, m, output = kind)
  stopifnot(identical(before, serialize(list(g, m), NULL)))
  stopifnot(identical(actual, oracle(g, m, output = kind)))
  rows[[length(rows) + 1L]] <- data.frame(
    case = "ordinary-128x25",
    output = kind,
    pairs = length(g) * length(m),
    identical = TRUE
  )
}
z <- readRDS(file.path(out, "fixture-inputs.rds"))
before <- serialize(z, NULL)
invisible(native(z$g, z$m, output = "match"))
stopifnot(identical(before, serialize(z, NULL)))
write.csv(
  do.call(rbind, rows),
  file.path(out, "ordinary-regression.csv"),
  row.names = FALSE
)
cat("ORDINARY REGRESSION AND INPUT IMMUTABILITY PASSED\n")
