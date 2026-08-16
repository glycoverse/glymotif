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
  ...,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient"),
  strict_floating = TRUE
)

.g_count_motif(
  glycan_graph,
  motif_graph,
  ...,
  alignment = "substructure",
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient"),
  strict_floating = TRUE
)

.g_match_motif(
  glycan_graph,
  motif_graph,
  ...,
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

- ...:

  These dots must be empty and are used only to force optional arguments
  to be supplied by name.

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

- strict_floating:

  A logical scalar. For `.g_have_motif()`, `TRUE` requires the motif in
  every possible floating localization and `FALSE` requires it in at
  least one. For `.g_count_motif()`, `TRUE` returns the minimum count
  across localizations and `FALSE` returns the maximum.

## Value

- `.g_have_motif()` returns a logical scalar.

- `.g_count_motif()` returns an integer scalar.

- `.g_match_motif()` returns a list of integer vectors.

## Details

These functions do no validation, parsing, naming, or graph mutation.
Callers must provide valid graph objects. Residue compatibility follows
the high-level matching rules, including generic and mixed motif
residues.

These functions never call
[`glyrepr::as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.html).

Glycan graphs with unresolved floating parts or substituents are matched
across all conflict-free localizations. `.g_match_motif()` returns the
union of mappings from every localization, with node indices referring
to the original unresolved graph.

## See also

[`have_motif()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
[`count_motif()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
[`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md)

## Examples

``` r
library(glyparse)
library(glyrepr)

glycan <- parse_iupac_condensed("Gal(b1-3)GalNAc(b1-")
motif <- parse_iupac_condensed("Gal(b1-")
glycan_graph <- get_structure_graphs(glycan)
motif_graph <- get_structure_graphs(motif)

.g_have_motif(glycan_graph, motif_graph)
#> [1] TRUE
.g_count_motif(glycan_graph, motif_graph)
#> [1] 1
.g_match_motif(glycan_graph, motif_graph)
#> [[1]]
#> [1] 1
#> 
```
