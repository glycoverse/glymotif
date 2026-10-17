# Branch Motifs Specification

Create a specification for branch motif extraction. This should be
passed to the `motifs` argument of
[`have_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
or
[`match_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md).

## Usage

``` r
branch_motifs(include_empty_branch = FALSE, include_mannose_branch = FALSE)
```

## Arguments

- include_empty_branch:

  A logical scalar. If `TRUE`, include bare core a3/a6 mannose arms as
  empty branches. Default is `FALSE`.

- include_mannose_branch:

  A logical scalar. If `TRUE`, include full subtrees rooted at Man (or
  generic Hex) directly attached to core a3/a6 mannose. Default is
  `FALSE`.

## Value

A `branch_motifs_spec` object.

## Details

Passing `branch_motifs()` to the `motifs` argument of supported
functions will:

1.  Call
    [`extract_branch_motif()`](https://glycoverse.github.io/glymotif/dev/reference/extract_branch_motif.md)
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

With `include_mannose_branch = TRUE`, Man-rooted subtrees are extracted
and matched in the same way as GlcNAc-rooted subtrees. Identical branch
labels with different source cores receive numeric suffixes to
distinguish the motifs. With `include_empty_branch = TRUE`, bare core
arms are labelled `empty` (or `empty (generic)` or `empty (mixed)` for
other core types). Substituents are preserved and included in the label,
e.g. `empty (Man6S)`. Their arm residue must have no outward
monosaccharides.
[`match_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md)
returns that core arm's node index for each empty-branch match.

Therefore, `have_motifs(glycans, branch_motifs())` doesn't equal to
`have_motifs(glycans, extract_branch_motif(glycans))`. Never use the
results from
[`extract_branch_motif()`](https://glycoverse.github.io/glymotif/dev/reference/extract_branch_motif.md)
directly in these functions.

## See also

[`dynamic_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/dynamic_motifs.md),
[`extract_branch_motif()`](https://glycoverse.github.io/glymotif/dev/reference/extract_branch_motif.md)

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
count_motifs(
  glycans,
  branch_motifs(include_empty_branch = TRUE, include_mannose_branch = TRUE)
)
#>                                                                             GlcNAc(b1-
#> GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-          2
#> Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-             0
#>                                                                             Gal(b1-4)GlcNAc(b1-
#> GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-                   0
#> Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-                      1
#>                                                                             empty
#> GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-     0
#> Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-        1
```
