# Get the Structures or Alignments of Known Motifs

**\[deprecated\]**

`get_motif_structure()` and `get_motif_alignment()` were deprecated in
glymotif 0.16.0. Use
[`db_motif_info()`](https://glycoverse.github.io/glymotif/dev/reference/db_motif_info.md)
to inspect database motifs instead.

## Usage

``` r
get_motif_structure(name)

get_motif_alignment(name)
```

## Arguments

- name:

  A character vector of motif names.

## Value

- `get_motif_structure()`: a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector.

- `get_motif_alignment()`: a character vector of motif alignments.

For `get_motif_alignment()`, if `name` has length greater than 1, the
return value is named with the motif names.

## See also

[`db_motif_info()`](https://glycoverse.github.io/glymotif/dev/reference/db_motif_info.md)

## Examples

``` r
get_motif_structure("LacdiNAc")
#> Warning: `get_motif_structure()` was deprecated in glymotif 0.16.0.
#> ℹ Please use `db_motif_info()` instead.
#> <glycan_structure[1]>
#> [1] GalNAc(b1-4)GlcNAc(?1-
#> # Unique structures: 1
get_motif_alignment("LacdiNAc")
#> Warning: `get_motif_alignment()` was deprecated in glymotif 0.16.0.
#> ℹ Please use `db_motif_info()` instead.
#> [1] "substructure"

get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#> O-Glycan core 1 O-Glycan core 2 
#>          "core"          "core" 
```
