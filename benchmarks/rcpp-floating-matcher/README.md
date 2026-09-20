# Floating structures: C++ object-to-result matching experiment

This extends the earlier [ordinary-structure experiment](../rcpp-matcher/REPORT.md).
It is independent of the production `glymotif` API and retains the earlier files
and measurements unchanged. Entry point: `cpp_floating_match()`.

From the **glymotif repository root**, run in this order:

```sh
Rscript benchmarks/rcpp-floating-matcher/audit.R
Rscript benchmarks/rcpp-floating-matcher/fixtures.R
Rscript benchmarks/rcpp-floating-matcher/guards.R
Rscript benchmarks/rcpp-floating-matcher/corpus.R
Rscript benchmarks/rcpp-floating-matcher/benchmark.R
Rscript benchmarks/rcpp-floating-matcher/regression.R
python3 benchmarks/rcpp-floating-matcher/report.py
```

Rcpp/BH compilation occurs in a temporary cache. `setup.R` loads the current
checkout as the reference and guards igraph version 2.3.3. The C++ prototype
reads igraph's internal arrays; this remains an experimental representation
boundary, not a stable public igraph ABI.

The C++ entry receives the original `glycan_structure` vector, reads cached
forests and floating metadata, enumerates candidate-parent combinations,
checks cycles and joint carbon-slot feasibility, constructs localized native
graphs without renumbering vertices, performs matching, aggregates and returns
R matrices/lists. It makes no R function callbacks. It keeps the current 256
raw-combination limit and rejects unresolved floating motifs as the public API
does. It trusts glyrepr-created input objects, not arbitrary forged S3 lists.

All/any presence and min/max counts follow `strict_floating`. Mapping output
always takes the union across localizations, deduplicating mapping vectors in
first-occurrence order. Presence short-circuits across variants once the result
is determined; current R computes all variants first. Unlike presence/count APIs, the public `match_motifs()`
has no `strict_floating` argument.

`audit.R` compares every localized graph, attachment assignment and error
boundary for all 564 floating corpus structures. `fixtures.R` exercises matching
options; `guards.R` covers cycles, collisions and the combination limit.
`corpus.R` reuses R localization/matching results solely for the correctness
oracle, checks all six have/count/match × floating-mode outputs, and also checks
against public APIs on a fixed subset. `benchmark.R` does **not** use those cached
localizations: each timed public/native call starts with the original objects.

`inputs.rds`, `fixture-inputs.rds` and `corpus-results.rds` freeze the inputs and
reference outputs. CSV files preserve localization checks, matching checks and
raw timings; `REPORT.md` summarizes measured results and limitations.
