# Get All Motifs

This function returns the names of all motifs available in the package.
We use GlycoMotif GlyGen Collection
(https://glycomotif.glyomics.org/glycomotif/GGM) as the source of the
motifs. This function is useful to be integrated with
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md).
For example, use `have_motifs(glycans, all_motifs())` to check against
all motifs.

## Usage

``` r
all_motifs()
```

## Value

A character vector of motif names.

## Examples

``` r
all_motifs()[1:5]
#> [1] "Blood group H (type 2) - Lewis y" "i antigen"                       
#> [3] "LacdiNAc"                         "GT2"                             
#> [5] "Blood group B (type 1) - Lewis b"
```
