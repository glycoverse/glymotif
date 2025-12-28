# Extract Branch Motifs

An N-glycan branching motif if the substructures linked to either the
a3- or a6-core-mannose.

For example:

    Neu5Ac - Gal - GlcNAc - Man
    ~~~~~~~~~~~~~~~~~~~~~      \
      A branching motif         Man - GlcNAc - GlcNAc -
                               /
                            Man

This function returns all the unique branching motifs found in the input
glycans. It is very useful for analyzing N-glycan antennary patterns
combined with `glydet::quantify_motifs()`. When used with
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
or related functions, we recommend setting `alignment` to "exact".

## Usage

``` r
extract_branch_motif(glycans)
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
```
