# Check if a Motif is Known

**\[deprecated\]**

`is_known_motif()` was deprecated in glymotif 0.16.0. Use
[`db_motif_info()`](https://glycoverse.github.io/glymotif/reference/db_motif_info.md)
to inspect database motifs instead.

## Usage

``` r
is_known_motif(name)
```

## Arguments

- name:

  A character vector of motif names.

## Value

A logical vector.

## Examples

``` r
is_known_motif(c("N-Glycan core basic", "O-Glycan core 1", "unknown"))
#> Warning: `is_known_motif()` was deprecated in glymotif 0.16.0.
#> ℹ Please use `db_motif_info()` instead.
#> [1]  TRUE  TRUE FALSE
```
