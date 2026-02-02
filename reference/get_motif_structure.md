# Get the Structures, Alignments, or Aglycons of Known Motifs

Given a character vector of motifs names in GlycoMotif GlyGen
Collection, these functions return the structures, alignments, or
aglycons of the motifs.

## Usage

``` r
get_motif_structure(name)

get_motif_alignment(name)

get_motif_aglycon(name)
```

## Arguments

- name:

  A character vector of the motif name.

## Value

`get_motif_structure()` returns a `glycan_graph` when `name` is a
character scalar, and a list of `glycan_graph` when `name` is a
character vector. `get_motif_alignment()` returns a character vector of
motif alignments. `get_motif_aglycon()` returns a character vector of
motif aglycons. For all three functions, if `name` has length greater
than 1, the return value is named with the motif names.

## See also

[`db_motifs()`](https://glycoverse.github.io/glymotif/reference/db_motifs.md),
[`is_known_motif()`](https://glycoverse.github.io/glymotif/reference/is_known_motif.md)

## Examples

``` r
get_motif_structure("N-Glycan core basic")
#> <glycan_structure[1]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-
#> # Unique structures: 1
get_motif_alignment("N-Glycan core basic")
#> [1] "core"
get_motif_aglycon("N-Glycan core basic")
#> [1] "Asn"

get_motif_structure(c("O-Glycan core 1", "O-Glycan core 2"))
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
get_motif_alignment(c("O-Glycan core 1", "O-Glycan core 2"))
#> O-Glycan core 1 O-Glycan core 2 
#>          "core"          "core" 
get_motif_aglycon(c("O-Glycan core 1", "O-Glycan core 2"))
#> O-Glycan core 1 O-Glycan core 2 
#>       "Ser/Thr"       "Ser/Thr" 
```
