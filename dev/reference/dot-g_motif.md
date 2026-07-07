# Low-Level Motif Matching on Graphs

These functions are low-level variants of
[`have_motif()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
[`count_motif()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
and
[`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md)
for package code that already has compatible igraph objects from
[`glyrepr::get_structure_graphs()`](https://glycoverse.github.io/glyrepr/reference/get_structure_graphs.html).

## Usage

``` r
.g_have_motif(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
)

.g_count_motif(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
)

.g_match_motif(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
)
```

## Arguments

- glycan_graph:

  An igraph glycan graph.

- motif_graph:

  An igraph motif graph.

- alignment:

  A character scalar: `"substructure"`, `"core"`, `"terminal"`, or
  `"whole"`.

- ignore_linkages:

  A logical scalar. If `TRUE`, linkages are ignored.

- strict_sub:

  A logical scalar. If `TRUE`, substituents are matched strictly.

- match_degree:

  A logical vector indicating which motif nodes must match the glycan's
  in- and out-degree exactly. A scalar is recycled to the number of
  motif nodes.

- mode:

  Matching mode. `"strict"` preserves the default behavior; `"lenient"`
  treats glycan-side unknowns as compatible with more specific motif
  fields.

## Value

- `.g_have_motif()` returns a logical scalar.

- `.g_count_motif()` returns an integer scalar.

- `.g_match_motif()` returns a list of integer vectors.

## Details

These functions do no validation, parsing, naming, or
monosaccharide-type conversion. Callers must provide valid,
already-compatible graph objects. In strict mode, when matching a
generic motif graph, callers must pass a glycan graph whose `mono`
vertex attributes have already been converted to the compatible generic
representation.

These functions never call
[`glyrepr::as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.html).

## See also

[`have_motif()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
[`count_motif()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
[`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md)
