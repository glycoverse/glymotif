# Dynamic Motifs Specification

Create a specification for dynamic motif extraction. This should be
passed to the `motifs` argument of
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
[`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
[`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md),
or
[`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md).

## Usage

``` r
dynamic_motifs(max_size = 3)
```

## Arguments

- max_size:

  The maximum number of monosaccharides in the extracted motifs. Default
  is 3. Passed to
  [`extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.md).

## Value

A `dynamic_motifs_spec` object.

## Details

Passing `dynamic_motifs()` to the `motifs` argument of supported
functions will:

1.  Call
    [`extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.md)
    on `glycans` to get all dynamic motifs.

2.  Perform motif matching with `alignments` as "substructure".

In fact, `have_motifs(glycans, dynamic_motifs())` is just a syntatic
sugar of `have_motifs(glycans, extract_motif(glycans))`. This function
exists to align with the
[`db_motifs()`](https://glycoverse.github.io/glymotif/reference/db_motifs.md)
and
[`branch_motifs()`](https://glycoverse.github.io/glymotif/reference/branch_motifs.md)
API.

## See also

[`branch_motifs()`](https://glycoverse.github.io/glymotif/reference/branch_motifs.md),
[`extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.md)

## Examples

``` r
library(glyrepr)
glycans <- c(o_glycan_core_1(), o_glycan_core_2())
have_motifs(glycans, dynamic_motifs())
#>                                   Gal(b1- GalNAc(a1- Gal(b1-3)GalNAc(a1-
#> Gal(b1-3)GalNAc(a1-                  TRUE       TRUE                TRUE
#> Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-    TRUE       TRUE                TRUE
#>                                   GlcNAc(b1- GlcNAc(b1-6)GalNAc(a1-
#> Gal(b1-3)GalNAc(a1-                    FALSE                  FALSE
#> Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-       TRUE                   TRUE
#>                                   Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> Gal(b1-3)GalNAc(a1-                                           FALSE
#> Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-                              TRUE
```
