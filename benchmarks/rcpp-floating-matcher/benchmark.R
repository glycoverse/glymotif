source("benchmarks/rcpp-floating-matcher/setup.R")
saved <- readRDS(file.path(out, "corpus-results.rds"))
g <- saved$g
m <- saved$m
set.seed(20260920)
rows <- list()
measure <- function(f, n) {
  gc()
  t0 <- proc.time()
  for (i in seq_len(n)) {
    invisible(f())
  }
  dt <- proc.time() - t0
  c(
    elapsed = unname(dt["elapsed"]) / n,
    cpu = unname(dt["user.self"] + dt["sys.self"]) / n
  )
}
benchmark <- function(
  label,
  g,
  m,
  kind = "have",
  sf = TRUE,
  reps = 3L,
  iterations = 1L
) {
  stopifnot(identical(
    oracle(g, m, strict_floating = sf, output = kind),
    validated(g, m, strict_floating = sf, output = kind)
  ))
  funcs <- list(
    current = function() oracle(g, m, strict_floating = sf, output = kind),
    cpp_direct = function() native(g, m, strict_floating = sf, output = kind),
    cpp_validated = function() {
      validated(g, m, strict_floating = sf, output = kind)
    }
  )
  # Parity call above warms current and validated routes; warm the direct route.
  invisible(funcs$cpp_direct())
  for (r in seq_len(reps)) {
    for (engine in sample(names(funcs))) {
      tm <- measure(funcs[[engine]], iterations)
      rows[[length(rows) + 1L]] <<- data.frame(
        workload = label,
        output = kind,
        strict_floating = sf,
        glycans = length(g),
        motifs = length(m),
        engine = engine,
        rep = r,
        iterations = iterations,
        elapsed_s = tm["elapsed"],
        cpu_s = tm["cpu"]
      )
      write.csv(
        do.call(rbind, rows),
        file.path(out, "timings.csv"),
        row.names = FALSE
      )
      cat(label, kind, sf, r, engine, tm["elapsed"], "s\n")
    }
  }
}
# Size-stratified deterministic subset of the 503 supported structures.
a <- read.csv(file.path(out, "localization-audit.csv"))
a <- a[a$result == "identical", ]
idx <- order(a$variants)[unique(round(seq(1, nrow(a), length.out = 32)))]
for (kind in c("have", "count", "match")) {
  benchmark("32x6", g[idx], m, kind = kind)
}
benchmark("32x6-any", g[idx], m, sf = FALSE)
benchmark("503x6", g, m)
summary <- aggregate(
  cbind(elapsed_s, cpu_s) ~ workload + output + strict_floating + engine,
  do.call(rbind, rows),
  median
)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
cat("BENCHMARK DONE\n")
