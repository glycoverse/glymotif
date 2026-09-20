# Object-to-result Rcpp matcher experiment

This is a standalone experiment stored in the **glymotif** repository. It does not change either package's production API, dependencies or
matching implementation.

Run from the glymotif root:

```r
source("benchmarks/rcpp-matcher/run.R")
source("benchmarks/rcpp-matcher/supplement.R")
```

The scripts load the adjacent glymotif checkout as the oracle, compiles
`matcher.cpp` with Rcpp/BH, runs differential checks, and writes paired timings.
The sourceCpp build is temporary. `inputs.rds` freezes sampled structure keys and
motifs; `parity.csv`, `timings.csv`, `summary.csv` and `sessionInfo.txt` record the
results. `run.log` contains progress. See `REPORT.md` for measured conclusions. Rebuild the report and numeric ranges
with `python3 benchmarks/rcpp-matcher/report.py` after running both R scripts.

## Boundaries

- Native entry: `glyrepr_structure` vectors, with their cached `graphs` attributes.
- C++ reads the object keys and cached igraph arrays, deduplicates glycans, builds
  profiles and compatibility matrices, runs Boost VF2, deduplicates mappings,
  restores duplicates/missing values and builds the R output.
- Matching supports strict/lenient residue types, fuzzy and built-in substituent
  modifications, alternative/unknown linkage positions, root anomers, the four
  alignments, degree masks, presence, counts and full node mappings.
- No R callbacks occur inside `cpp_structure_match()`. The static glyrepr residue
  dictionary is supplied as package data; copying it into C++ is timed each call.
- `cpp_direct` measures the typed-object native route. `cpp_validated` additionally
  retains current `prepare_motif_args()` validation, conversion and warning logic.
  The latter is the conservative estimate for integrating the prototype under the
  public API. All three routes start with the same already-created objects.
- Timings include graph extraction, construction, matching and returned result
  allocation; exclude input parsing, compilation and loading static package data.
  No graph profiles or compatibility matrices are reused between timed calls.
- This prototype explicitly rejects unresolved floating parts/substituents. It
  does not silently drop them or call R as a fallback. Their localization and
  cross-localization aggregation remain a separate, unported algorithm.
- Direct igraph layout access is **experimental and version-dependent**. Tested
  layout: igraph 2.3.3, 10 list slots, attribute slot 9, zero-based endpoints in
  slots 3/4. Production use needs a stable native bridge or owned compact graph
  representation, with version/shape safeguards and upstream compatibility tests.
- The direct entry expects valid glyrepr-created rooted trees. Its guards are not
  a replacement for all public input validation; malformed handcrafted S3 objects,
  character parsing, motif specifications and database-name resolution are outside
  this native experiment. The validated route retains those R-side entry rules.

`matcher.cpp` starts with the existing MIT-licensed glymotif VF2 implementation;
this intentionally keeps the search engine fixed while moving preparation and
batch/result handling into C++. The experiment does not claim a faster new VF2.
