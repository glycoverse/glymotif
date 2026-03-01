# Check if a Motif is Known

This function checks if motifs are known motifs in GlycoMotif GlyGen
Collection.

## Usage

``` r
is_known_motif(name)
```

## Arguments

- name:

  A character vector of the motif name.

## Value

A logical vector.

## Examples

``` r
is_known_motif(c("N-Glycan core basic", "O-Glycan core 1", "unknown"))
#> [1]  TRUE  TRUE FALSE
```
