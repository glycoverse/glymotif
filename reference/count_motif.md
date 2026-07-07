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
  ...,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
)

count_motifs(
  glycans,
  motifs,
  ...,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  mode = c("strict", "lenient")
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

  - A GGM database motif name (use
    `db_motif_info() |> dplyr::filter(source_id == "GGM")` to inspect
    resolvable motif names).

- ...:

  These dots must be empty and are used only to force optional arguments
  to be supplied by name.

- alignment:

  A character string. Possible values are "substructure", "core",
  "terminal", and "whole". If not provided, the value will be decided
  based on the `motif` argument. If `motif` is a GGM motif name, the
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

- mode:

  Matching mode. `"strict"` preserves the default behavior where glycans
  cannot be more obscure than motifs. `"lenient"` treats glycan-side
  obscure fields as compatible with more specific motif fields while
  still rejecting concrete mismatches.

- motifs:

  One of:

  - A
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector.

  - A glycan structure string vector, supported by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - A character vector of GGM database motif names (use
    `db_motif_info() |> dplyr::filter(source_id == "GGM")` to inspect
    resolvable motif names).

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
and `count_motif()` perserve names from the input `glycans` vector.

[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and `count_motifs()` return a matrix with both row and column names. The
row names are the glycan names, and the column names are the motif
names.

Glycan names follow the same rule as
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and `count_motif()`.

Motif names have the following rules:

1.  If `motifs` have names, use the names.

2.  If `motifs` don't have names and are GGM database motif names (e.g.
    "N-glycan core"), use them.

3.  Otherwise, no colnames.

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
#> Warning: Matching lower-level `glycans` against higher-level `motifs` usually returns no
#> matches.
#> ℹ `glycans` have "basic" structure level, while `motifs` have "partial"
#>   structure level.
#> ℹ Use motifs at the same structure level as the glycans, or reduce motif
#>   structure levels before matching.
#> ℹ See `?get_structure_level` for details.
#> [1] 0
```
