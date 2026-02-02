# Branch Motifs Specification

Create a specification for branch motif extraction. This should be
passed to the `motifs` argument of
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
[`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
[`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md),
or
[`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md).

## Usage

``` r
branch_motifs()
```

## Value

A `branch_motifs_spec` object.

## Details

When used, motifs will be extracted dynamically from the input glycans
using
[`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)
with `including_core = TRUE`. Motif matching will use a custom
`match_degree` where the last 4 nodes of each motif are not required to
match degree exactly (to allow for attachment to the glycan core).

## See also

[`dynamic_motifs()`](https://glycoverse.github.io/glymotif/reference/dynamic_motifs.md),
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)

## Examples

``` r
# Use in have_motifs()
# have_motifs(glycans, branch_motifs())
```
