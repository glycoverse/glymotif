# Add Motif Annotations

**\[deprecated\]**

`add_motifs_int()` and `add_motifs_lgl()` were deprecated in glymotif
0.17.0. For data frames, use
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
with
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
and
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
or
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md).
For
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
or
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
objects, use
[`glyexp::mutate_row()`](https://glycoverse.github.io/glyexp/reference/mutate_col.html)
with the same tibble expression. The deprecated methods for legacy
glyexp containers remain available for compatibility.

## Usage

``` r
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)

add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)

# S3 method for class 'glyexp_experiment'
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)

# S3 method for class 'glyexp_experiment'
add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)

# S3 method for class 'data.frame'
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)

# S3 method for class 'data.frame'
add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL,
  ...
)
```

## Arguments

- x:

  A legacy glyexp data container or a tibble with a structure column.

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
  [`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`match_motif()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
  this must be a logical vector with length 1 or the number of motif
  nodes (length 1 is recycled). For
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
  this must be a list of logical vectors with length equal to `motifs`;
  each element follows the same length rules. When `match_degree` is
  provided, `alignment` and `alignments` are silently ignored.

- ...:

  Additional arguments passed to the method.

## Value

The legacy glyexp data container or data frame supplied in `x`, with
motif annotations added to its variable information.

## About Names

The naming rule for the new columns follows these priorities:

1.  If `motifs` is a named vector (character or
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)),
    the names are used directly as column names.

2.  If `motifs` is unnamed and contains GGM database motif names (e.g.,
    "N-Glycan core"), the motif names are used as column names.

3.  If `motifs` is unnamed and contains
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    objects or IUPAC-condensed structure strings, the IUPAC-condensed
    strings are used as column names.

4.  If `motifs` is a motif spec from
    [`dynamic_motifs()`](https://glycoverse.github.io/glymotif/reference/dynamic_motifs.md)
    or
    [`branch_motifs()`](https://glycoverse.github.io/glymotif/reference/branch_motifs.md),
    the IUPAC-condensed strings of the extracted motifs are used as
    column names.

Note: This behavior differs from
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
and
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
which return matrices with NULL column names for unnamed IUPAC string or
structure motifs. The functions here always provide column names since
they are designed for adding motif annotations to data frames.

## See also

[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html),
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)

## Examples

``` r
library(glyexp)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(tibble)

glyco_se <- real_experiment2
motifs <- c(
  lacnac = "Gal(??-?)GlcNAc(??-",
  sia_lacnac = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-"
)

glyco_se |>
  mutate_row(as_tibble(count_motifs(glycan_structure, motifs)))
#> 
#> ── GlycomicSE ──────────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 144 samples, 67 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: glycan_composition <comp>, glycan_structure <struct>, lacnac <int>, sia_lacnac <int>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>

df <- tibble(
  glycan_structure = c(
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
  )
)
df |>
  mutate(as_tibble(count_motifs(glycan_structure, motifs)))
#> # A tibble: 2 × 3
#>   glycan_structure                lacnac sia_lacnac
#>   <chr>                            <int>      <int>
#> 1 Gal(b1-4)GlcNAc(b1-                  1          0
#> 2 Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-      1          1
```
