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

When used, motifs will be extracted dynamically from the input glycans
using
[`extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.md),
and motif matching will use `alignment = "substructure"`.

## See also

[`branch_motifs()`](https://glycoverse.github.io/glymotif/reference/branch_motifs.md),
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)

## Examples

``` r
# Use in have_motifs()
# have_motifs(glycans, dynamic_motifs(max_size = 4))
```
