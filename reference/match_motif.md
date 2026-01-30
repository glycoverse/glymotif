# Match Motif(s) in Glycans

These functions find all occurrences of the given `motif`(s) in the
`glycans`. Node-to-node mapping is returned for each match. This
function is NOT useful for most users if you are not interested in the
concrete node mapping. See
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and
[`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
for more information about the matching rules.

- `match_motif()` matches a single motif against multiple glycans

- `match_motifs()` matches multiple motifs against multiple glycans

Different from
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and
[`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
these functions return detailed match information. More specifically,
for each glycan-motif pair, a integer vector is returned, indicating the
node mapping from the motif to the glycan. For example, if the vector is
`c(2, 3, 6)`, it means that the first node in the motif matches the 2nd
node in the glycan, the second node in the motif matches the 3rd node in
the glycan, and the third node in the motif matches the 6th node in the
glycan.

Node indices are only meaningful for
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
so only
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
is supported for `glycans` and `motifs`.

## Usage

``` r
match_motif(
  glycans,
  motif,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE
)

match_motifs(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE
)
```

## Arguments

- glycans:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector. All formats supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)
    are accepted, including IUPAC-condensed, WURCS, GlycoCT, and others.

- motif:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    scalar.

  - A glycan structure string, supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - A known motif name (use
    [`all_motifs()`](https://glycoverse.github.io/glymotif/reference/all_motifs.md)
    to see all available motifs).

- alignment:

  A character string. Possible values are "substructure", "core",
  "terminal", "whole", and "exact". If not provided, the value will be
  decided based on the `motif` argument. If `motif` is a motif name, the
  alignment in the database will be used. Otherwise, "substructure" will
  be used.

- ignore_linkages:

  A logical value. If `TRUE`, linkages will be ignored in the
  comparison. Default is `FALSE`.

- strict_sub:

  A logical value. If `TRUE` (default), substituents will be matched in
  strict mode, which means if the glycan has a substituent in some
  residue, the motif must have the same substituent to be matched.

- motifs:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector, supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - A character vector of motif names (use
    [`all_motifs()`](https://glycoverse.github.io/glymotif/reference/all_motifs.md)
    to see all available motifs).

- alignments:

  A character vector specifying alignment types for each motif. Can be a
  single value (applied to all motifs) or a vector of the same length as
  motifs.

## Value

A nested list of integer vectors.

- `match_motif()`: Two levels of nesting. The outer list corresponds to
  glycans, and the inner list corresponds to matches. Use
  `purrr::pluck(result, glycan_index, match_index)` to access the match
  information. For example, `purrr::pluck(result, 1, 2)` means the 2nd
  match in the 1st glycan.

- `match_motifs()`: Three levels of nesting. The outermost list
  corresponds to motifs, the middle list corresponds to glycans, and the
  innermost list corresponds to matches. Use
  `purrr::pluck(result, motif_index, glycan_index, match_index)` to
  access the match information. For example,
  `purrr::pluck(result, 1, 2, 3)` means the 3rd match in the 2nd glycan
  for the 1st motif. The outermost list is named by `motifs` if they
  have names. The middle list is named by `glycans` if they have names.

## Vertex and Linkage Indices

The indices of vertices and linkages in a glycan correspond directly to
their order in the IUPAC-condensed string, which is printed when you
print a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html).
For example, for the glycan
`Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-)`, the vertices are
"Man", "Man", "Man", "GlcNAc", "GlcNAc", and the linkages are "a1-3",
"a1-6", "b1-4", "b1-4".

Thus, matching the motif "Man(a1-3)Man(b1-4)" to this glycan yields
`c(1, 3)`. This indicates that the first motif vertex (the a1-3 Man)
corresponds to the first vertex in the glycan, and the second motif
vertex (the b1-4 Man) corresponds to the third vertex in the glycan.

## See also

[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)

## Examples

``` r
library(glyparse)
library(glyrepr)

(glycan <- n_glycan_core())
#> <glycan_structure[1]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1

# Let's peek under the hood of the nodes in the glycan
glycan_graph <- get_structure_graphs(glycan)
igraph::V(glycan_graph)$mono # 1, 2, 3, 4, 5
#> [1] "Man"    "Man"    "Man"    "GlcNAc" "GlcNAc"

# Match a single motif against a single glycan
motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")
match_motif(glycan, motif)
#> [[1]]
#> [[1]][[1]]
#> [1] 1 2 3
#> 
#> 

# Match multiple motifs against a single glycan
motifs <- c(
  "Man(a1-3)[Man(a1-6)]Man(b1-",
  "Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
)
motifs <- parse_iupac_condensed(motifs)
match_motifs(glycan, motifs)
#> [[1]]
#> [[1]][[1]]
#> [[1]][[1]][[1]]
#> [1] 1 2 3
#> 
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [[2]][[1]][[1]]
#> [1] 1 3 4 5
#> 
#> 
#> 
```
