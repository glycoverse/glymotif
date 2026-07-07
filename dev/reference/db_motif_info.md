# Get Database Motif Information

Returns metadata for all motifs available in the package. You can use
`dplyr::distinct(db_motif_info(), source_id, source)` to get all
available sources.

## Usage

``` r
db_motif_info()
```

## Value

A tibble.

## Details

It contains the following columns:

- `source_id`: the collection identifier of the motif

- `source`: the collection name of the motif

- `accession`: the accession number of the motif

- `name`: the name of the motif

- `alignment`: the alignment of the motif

- `glycan_structure`: the glycan structure (glyrepr::glycan_structure())
  of the motif

## Examples

``` r
db_motif_info()
#> # A tibble: 2,843 × 6
#>    source      source_id accession name               alignment glycan_structure
#>    <chr>       <chr>     <chr>     <chr>              <chr>     <struct>        
#>  1 CCRC Motifs CCRC      000019    Fuc(a1-2)Gal(b1-4… substruc… Fuc(a1-2)Gal(b1…
#>  2 CCRC Motifs CCRC      000075    GlcNAc(b1-3)[GalN… core      GlcNAc(b1-3)[Ga…
#>  3 CCRC Motifs CCRC      000014    Neu5Ac(a2-3)Gal(b… core      Neu5Ac(a2-3)Gal…
#>  4 CCRC Motifs CCRC      000022    Neu5Ac(a2-3)Gal(b… substruc… Neu5Ac(a2-3)Gal…
#>  5 CCRC Motifs CCRC      000078    Fuc(a1-2)Gal(b1-3… core      Fuc(a1-2)Gal(b1…
#>  6 CCRC Motifs CCRC      000099    Gal(b1-3)GalNAc(b… core      Gal(b1-3)GalNAc…
#>  7 CCRC Motifs CCRC      000053    Neu5Ac(a2-3)Gal(b… core      Neu5Ac(a2-3)Gal…
#>  8 CCRC Motifs CCRC      000046    Fuc(a1-2)[Gal(a1-… substruc… Fuc(a1-2)[Gal(a…
#>  9 CCRC Motifs CCRC      000033    Neu5Ac(a2-3)Gal(b… substruc… Neu5Ac(a2-3)Gal…
#> 10 CCRC Motifs CCRC      000109    Neu5Ac(a2-8)Neu5A… core      Neu5Ac(a2-8)Neu…
#> # ℹ 2,833 more rows
```
