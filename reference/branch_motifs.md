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

Passing `branch_motifs()` to the `motifs` argument of supported
functions will:

1.  Call
    [`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)
    with `including_core = TRUE` on `glycans` to get all branching
    motifs.

2.  Construct a `match_degree` list based on the motifs.

3.  Perform motif matching using the constructed `match_degree` list.

Specifically, setting `including_core = TRUE` will include an additional
"Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-" suffix to each branching
motif. This suffix helps differentiate branching GlcNAc and bisecting
GlcNAc. Then, the `match_degree` is constructed so that the four
residues in the suffix do not have to match the node degree in the motif
matching process.

Therefore, `have_motifs(glycans, branch_motifs())` doesn't equal to
`have_motifs(glycans, extract_branch_motif(glycans))`. Never use the
results from
[`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)
directly in these functions.

## See also

[`dynamic_motifs()`](https://glycoverse.github.io/glymotif/reference/dynamic_motifs.md),
[`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)

## Examples

``` r
glycans <- c(
  "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
)
have_motifs(glycans, branch_motifs())
#>                                                                             GlcNAc(b1-
#> GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-       TRUE
#> Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-         FALSE
#>                                                                             Gal(b1-4)GlcNAc(b1-
#> GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-               FALSE
#> Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-                   TRUE
```
