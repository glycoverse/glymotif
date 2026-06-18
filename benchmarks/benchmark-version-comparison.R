#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: benchmark-version-comparison.R <package_path> <version_label> <output_csv>",
    call. = FALSE
  )
}

package_path <- normalizePath(args[[1]], mustWork = TRUE)
version_label <- args[[2]]
output_csv <- args[[3]]

bench_reps <- as.integer(Sys.getenv("GLYMOTIF_BENCH_REPS", "3"))
multi_glycan_n <- as.integer(Sys.getenv("GLYMOTIF_BENCH_MULTI_GLYCAN_N", "200"))

if (is.na(bench_reps) || bench_reps < 1) {
  stop("GLYMOTIF_BENCH_REPS must be a positive integer.", call. = FALSE)
}
if (is.na(multi_glycan_n) || multi_glycan_n < 1) {
  stop(
    "GLYMOTIF_BENCH_MULTI_GLYCAN_N must be a positive integer.",
    call. = FALSE
  )
}

suppressPackageStartupMessages(devtools::load_all(package_path, quiet = TRUE))

time_expr <- function(expr) {
  gc()
  elapsed <- system.time(force(expr))[["elapsed"]]
  unname(elapsed)
}

run_case <- function(fun_name, workload, pair_count, expr_fun) {
  timings <- vapply(
    seq_len(bench_reps),
    function(run) {
      time_expr(suppressMessages(expr_fun()))
    },
    numeric(1)
  )

  data.frame(
    version_label = version_label,
    package_version = as.character(utils::packageVersion("glymotif")),
    package_path = package_path,
    function_name = fun_name,
    workload = workload,
    run = seq_len(bench_reps),
    elapsed_sec = timings,
    pair_count = pair_count,
    stringsAsFactors = FALSE
  )
}

select_spread_indices <- function(n, size) {
  unique(pmax(1L, pmin(n, round(seq(1, n, length.out = size)))))
}

glycans <- glydb::glydb_structures()
motif_names <- db_motifs()
motif_alignments <- get_motif_alignment(motif_names)
alignment_classes <- c("substructure", "core", "terminal", "whole")

single_motif_names <- vapply(
  alignment_classes,
  function(alignment) {
    motif_names[motif_alignments == alignment][[1]]
  },
  character(1)
)

multi_glycans <- glycans[select_spread_indices(length(glycans), multi_glycan_n)]
multi_motif_structures <- get_motif_structure(motif_names)
multi_motif_alignments <- get_motif_alignment(motif_names)
unique_motif_structure_idx <- !duplicated(multi_motif_structures)
unique_motif_structures <- multi_motif_structures[unique_motif_structure_idx]
unique_motif_alignments <- multi_motif_alignments[unique_motif_structure_idx]
single_motif_structures <- get_motif_structure(single_motif_names)
single_motif_alignments <- get_motif_alignment(single_motif_names)

single_workload <- paste0(
  "all_glycans_x_",
  length(single_motif_names),
  "_representative_db_motifs"
)
multi_workload <- paste0(
  length(multi_glycans),
  "_spread_glycans_x_all_db_motifs"
)
unique_multi_workload <- paste0(
  length(multi_glycans),
  "_spread_glycans_x_unique_db_motif_structures"
)

single_pair_count <- length(glycans) * length(single_motif_names)
multi_pair_count <- length(multi_glycans) * length(motif_names)
unique_multi_pair_count <- length(multi_glycans) *
  length(unique_motif_structures)

results <- list(
  run_case(
    "have_motif",
    single_workload,
    single_pair_count,
    function() {
      purrr::walk(single_motif_names, ~ have_motif(glycans, .x))
    }
  ),
  run_case(
    "count_motif",
    single_workload,
    single_pair_count,
    function() {
      purrr::walk(single_motif_names, ~ count_motif(glycans, .x))
    }
  ),
  run_case(
    "match_motif",
    single_workload,
    single_pair_count,
    function() {
      purrr::walk2(
        single_motif_structures,
        single_motif_alignments,
        ~ match_motif(glycans, .x, alignment = .y)
      )
    }
  ),
  run_case(
    "have_motifs",
    multi_workload,
    multi_pair_count,
    function() {
      have_motifs(multi_glycans, motif_names)
    }
  ),
  run_case(
    "count_motifs",
    multi_workload,
    multi_pair_count,
    function() {
      count_motifs(multi_glycans, motif_names)
    }
  ),
  run_case(
    "match_motifs",
    unique_multi_workload,
    unique_multi_pair_count,
    function() {
      match_motifs(
        multi_glycans,
        unique_motif_structures,
        alignments = unique_motif_alignments
      )
    }
  )
)

results <- do.call(rbind, results)
results$glycan_total <- length(glycans)
results$motif_total <- length(motif_names)
results$single_motif_names <- paste(single_motif_names, collapse = " | ")

if (file.exists(output_csv)) {
  existing <- utils::read.csv(output_csv, stringsAsFactors = FALSE)
  results <- rbind(existing, results)
}

utils::write.csv(results, output_csv, row.names = FALSE)
print(results)
