# Get Database Motif Information

Returns metadata for all motifs available in the package.

## Usage

``` r
db_motif_info()
```

## Value

A tibble with `source`, `source_id`, `accession`, `name`, `alignment`,
and `glycan_structure` columns.

## Examples

``` r
db_motif_info()
#> # A tibble: 183 × 6
#>    source source_id accession name                    alignment glycan_structure
#>    <chr>  <chr>     <chr>     <chr>                   <chr>     <struct>        
#>  1 GlyGen GGM       000045    Blood group H (type 2)… substruc… Fuc(a1-2)Gal(b1…
#>  2 GlyGen GGM       000004    i antigen               substruc… Gal(b1-4)GlcNAc…
#>  3 GlyGen GGM       000016    LacdiNAc                substruc… GalNAc(b1-4)Glc…
#>  4 GlyGen GGM       000109    GT2                     whole     Neu5Ac(a2-8)Neu…
#>  5 GlyGen GGM       000046    Blood group B (type 1)… substruc… Fuc(a1-2)[Gal(a…
#>  6 GlyGen GGM       000075    LcGg4                   whole     GlcNAc(b1-3)[Ga…
#>  7 GlyGen GGM       000081    Sialosyl paragloboside  whole     Neu5Ac(a2-3)Gal…
#>  8 GlyGen GGM       000022    Sialyl Lewis x          terminal  Neu5Ac(a2-3)Gal…
#>  9 GlyGen GGM       000011    A antigen (type 3)      whole     Fuc(a1-2)[GalNA…
#> 10 GlyGen GGM       000006    Type 1 LN2              substruc… Gal(b1-3)GlcNAc…
#> # ℹ 173 more rows
```
