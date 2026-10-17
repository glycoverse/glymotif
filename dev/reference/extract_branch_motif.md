# Extract Branch Motifs

An N-glycan branching motif is a subtree linked to either the a3- or
a6-core-mannose.

For example:

    Neu5Ac - Gal - GlcNAc - Man
    ~~~~~~~~~~~~~~~~~~~~~      \
      A branching motif         Man - GlcNAc - GlcNAc -
                               /
                            Man

This function returns all unique branching motifs found in the input
glycans.

For branch motif matching with
[`have_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md)
or `glydet::quantify_motifs()`, use the
[`branch_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/branch_motifs.md)
helper, which handles the additional matching constraints.

## Usage

``` r
extract_branch_motif(
  glycans,
  ...,
  including_core = FALSE,
  include_empty_branch = FALSE,
  include_mannose_branch = FALSE
)
```

## Arguments

- glycans:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector. All formats supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)
    are accepted.

- ...:

  These dots must be empty and are used only to force optional arguments
  to be supplied by name.

- including_core:

  A logical scalar. If `TRUE`, the N-glycan core structure
  (`Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-` or
  `Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-`) is appended to each branch
  motif. Concrete branches receive a concrete core; generic and mixed
  branches receive a generic core. Default is `FALSE`. For Man-rooted
  branches, the source core's residue identities and substituents are
  preserved instead.

- include_empty_branch:

  A logical scalar. If `TRUE`, also extract empty branches: core a3/a6
  mannose residues with no outward monosaccharides. These are
  represented by the four-residue core path
  `Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-` (or its generic
  equivalent), even when `including_core = FALSE`. The source residue
  identities and substituents on this path are preserved. Default is
  `FALSE`.

- include_mannose_branch:

  A logical scalar. If `TRUE`, also extract full subtrees rooted at Man
  (or generic Hex) directly attached to a core a3/a6 mannose. Default is
  `FALSE`.

## Value

A
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector containing the unique extracted branching motifs.

## Details

The function works by:

1.  Converting the input to a set of unique `glycan_structure` objects.

2.  Searching for the N-glycan branch pattern:
    `HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-`.

3.  For each match, identifying the root node of the branch (the
    leftmost HexNAc in the pattern).

4.  Extracting the full subtree rooted at that node.

5.  Preserving the correct anomeric configuration (e.g., "b1") by
    inspecting the linkage to the root node.

With `include_mannose_branch = TRUE`, the same procedure also extracts
Man/Hex-rooted branches. With `include_empty_branch = TRUE`, bare core
arms are returned as core paths so they remain valid `glycan_structure`
objects.

## Examples

``` r
glycans <- c(
  "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-",
  "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(a1-4)GlcNAc(b1-"
)
extract_branch_motif(glycans)
#> <glycan_structure[2]>
#> [1] Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-
#> [2] Gal(b1-4)GlcNAc(b1-
#> # Unique structures: 2
extract_branch_motif(
  glycans, include_empty_branch = TRUE, include_mannose_branch = TRUE
)
#> <glycan_structure[3]>
#> [1] Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-
#> [2] Gal(b1-4)GlcNAc(b1-
#> [3] Man(??-?)Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> # Unique structures: 3
```
