# Count How Many Times Glycans have the Given Motif(s)

These functions are closely related to
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md).
However, instead of returning logical values, they return the number of
times the `glycans` have the `motif`(s).

- `count_motif()` counts a single motif in multiple glycans

- `count_motifs()` counts multiple motifs in multiple glycans

## Usage

``` r
count_motif(
  glycans,
  motif,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
)

count_motifs(
  glycans,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
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
    are accepted, including IUPAC-condensed, WURCS, GlycoCT, and others.

- motif:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    scalar.

  - A glycan structure string, supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - A known motif name (use
    [`db_motifs()`](https://glycoverse.github.io/glymotif/reference/db_motifs.md)
    to see all available motifs).

- alignment:

  A character string. Possible values are "substructure", "core",
  "terminal", and "whole". If not provided, the value will be decided
  based on the `motif` argument. If `motif` is a motif name, the
  alignment in the database will be used. Otherwise, "substructure" will
  be used.

- ignore_linkages:

  A logical value. If `TRUE`, linkages will be ignored in the
  comparison. Default is `FALSE`.

- strict_sub:

  A logical value. If `TRUE` (default), substituents will be matched in
  strict mode, which means if the glycan has a substituent in some
  residue, the motif must have the same substituent to be matched.

- match_degree:

  A logical vector indicating which motif nodes must match the glycan's
  in- and out-degree exactly. For
  [`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  `count_motif()`, and
  [`match_motif()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
  this must be a logical vector with length 1 or the number of motif
  nodes (length 1 is recycled). For
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  `count_motifs()`, and
  [`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
  this must be a list of logical vectors with length equal to `motifs`;
  each element follows the same length rules. When `match_degree` is
  provided, `alignment` and `alignments` are silently ignored.

- motifs:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector, supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - A character vector of motif names (use
    [`db_motifs()`](https://glycoverse.github.io/glymotif/reference/db_motifs.md)
    to see all available motifs).

- alignments:

  A character vector specifying alignment types for each motif. Can be a
  single value (applied to all motifs) or a vector of the same length as
  motifs.

## Value

- `count_motif()`: An integer vector indicating how many times each
  `glycan` has the `motif`.

- `count_motifs()`: An integer matrix where rows correspond to glycans
  and columns correspond to motifs. Row names contain glycan identifiers
  and column names contain motif identifiers.

## Details

This function actually perform v2f algorithm to get all possible matches
between `glycans` and `motif`. However, the result is not necessarily
the number of matches.

Think about the following example:

- glycan: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-4)GlcNAc(b1-`

- motif: `Gal(b1-?)[Gal(b1-?)]GlcNAc(b1-`

To draw the glycan out:

    Gal 1
       \ b1-? b1-4
        GlcNAc -- GlcNAc b1-
       / b1-?
    Gal 2

To draw the motif out:

    Gal 1
       \ b1-?
        GlcNAc b1-
       / b1-?
    Gal 2

To differentiate the galactoses, we number them as "Gal 1" and "Gal 2"
in both the glycan and the motif. The v2f subisomorphic algorithm will
return two matches:

- Gal 1 in the glycan matches Gal 1 in the motif, and Gal 2 matches Gal
  2.

- Gal 1 in the glycan matches Gal 2 in the motif, and Gal 2 matches Gal
  1.

However, from a biological perspective, the two matches are the same.
This function will take care of this, and return the "unique" number of
matches.

For other details about the handling of monosaccharide, linkages,
alignment, substituents, and implementation, see
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md).

## About Names

[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and `count_motif()` return a vector with no names. It is easy to trace
the names back to the original glycans.

[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and `count_motifs()` return a matrix with both row and column names. The
row names are the glycan names, and the column names are the motif
names. The names are decided according to the following rules:

1.  If `glycans` or `motifs` is a
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    object, the names are the IUPAC-condensed structure strings. (Sadly
    due to the constrains of the `vctrs` package
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    is built on, a
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector cannot have names.)

2.  If `glycans` or `motifs` is a character vector, either
    IUPAC-condensed structure strings or motif names, it will use the
    names of the character vector if exists, otherwise use the character
    vector itself as the names.

## See also

[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)

## Examples

``` r
library(glyparse)

count_motif("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-")
#> [1] 2
count_motif(
  "Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-",
  "Man(b1-?)[Man(b1-?)]GalNAc(b1-"
)
#> [1] 1
count_motif("Gal(b1-3)Gal(b1-", "Man(b1-")
#> [1] 0

# Vectorized usage with single motif
count_motif(c("Gal(b1-3)Gal(b1-3)GalNAc(b1-", "Gal(b1-3)GalNAc(b1-"), "Gal(b1-")
#> [1] 2 1

# Multiple motifs with count_motifs()
glycan1 <- parse_iupac_condensed("Gal(b1-3)Gal(b1-3)GalNAc(b1-")
glycan2 <- parse_iupac_condensed("Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-")
glycans <- c(glycan1, glycan2)

motifs <- c("Gal(b1-3)GalNAc(b1-", "Gal(b1-", "Man(b1-")
result <- count_motifs(glycans, motifs)
print(result)
#>                                            [,1] [,2] [,3]
#> Gal(b1-3)Gal(b1-3)GalNAc(b1-                  1    2    0
#> Man(b1-?)[Man(b1-?)]GalNAc(b1-4)GlcNAc(b1-    0    0    2

# Monosaccharide type matching examples
# Concrete glycan vs generic motif: matches (glycan converted to generic)
count_motif("Man(?1-", "Hex(?1-") # Returns 1
#> [1] 1

# Generic glycan vs concrete motif: doesn't match
count_motif("Hex(?1-", "Man(?1-") # Returns 0
#> [1] 0
```
