# Extract All Substructures (Motifs)

Extract all unique connected subgraphs (motifs) from the input glycans
up to a specified size. This function can be useful combined with
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
or `glydet::quantify_motifs()`. If so, set `alignment` to "substructure"
for these functions.

## Usage

``` r
extract_motif(glycans, max_size = 3)
```

## Arguments

- glycans:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector. All formats supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)
    are accepted.

- max_size:

  The maximum number of monosaccharides in the extracted motifs. Default
  is 3. Note that setting this value very large can be computationally
  expensive. Try the default value first, and increase it progressively
  if needed.

## Value

A
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector containing the unique extracted motifs.

## Examples

``` r
glycan <- "Gal(b1-3)[GlcNAc(a1-6)]GalNAc(a1-"
extract_motif(glycan, max_size = 2)
#> <glycan_structure[5]>
#> [1] Gal(b1-
#> [2] GlcNAc(a1-
#> [3] GalNAc(a1-
#> [4] Gal(b1-3)GalNAc(a1-
#> [5] GlcNAc(a1-6)GalNAc(a1-
#> # Unique structures: 5
```
