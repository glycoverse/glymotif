rows <- list()
time_one <- function(f, iterations) {
  gc()
  t0 <- proc.time()
  for (i in seq_len(iterations)) {
    invisible(f())
  }
  dt <- proc.time() - t0
  c(
    elapsed = unname(dt[["elapsed"]]) / iterations,
    cpu = unname(dt[["user.self"]] + dt[["sys.self"]]) / iterations
  )
}
benchmark <- function(label, g, m, kind = "have", iterations = 1L, reps = 7L) {
  stopifnot(identical(
    oracle(g, m, output = kind),
    validated(g, m, output = kind)
  ))
  fns <- list(
    current = function() oracle(g, m, output = kind),
    cpp_direct = function() native(g, m, output = kind),
    cpp_validated = function() validated(g, m, output = kind)
  )
  for (f in fns) {
    invisible(f())
  }
  for (r in seq_len(reps)) {
    order <- sample(names(fns))
    for (engine in order) {
      tm <- time_one(fns[[engine]], iterations)
      rows[[length(rows) + 1L]] <<- data.frame(
        workload = label,
        output = kind,
        glycans = length(g),
        unique_glycans = length(unique(g)),
        motifs = length(m),
        engine = engine,
        rep = r,
        iterations = iterations,
        elapsed_s = tm["elapsed"],
        cpu_s = tm["cpu"]
      )
      cat(label, kind, r, engine, round(tm["elapsed"], 5), "s\n")
      write.csv(
        do.call(rbind, rows),
        file.path(out, "timings.csv"),
        row.names = FALSE
      )
    }
  }
}
