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
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects, use
[`glyexp::mutate_var()`](https://glycoverse.github.io/glyexp/reference/mutate_obs.html)
with the same tibble expression.

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

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object, or a tibble with a structure column.

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

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with motif annotations added to the variable information.

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
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)

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

exp <- real_experiment2
motifs <- c(
  lacnac = "Gal(??-?)GlcNAc(??-",
  sia_lacnac = "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-"
)

exp |>
  mutate_var(as_tibble(have_motifs(glycan_structure, motifs))) |>
  get_var_info()
#> # A tibble: 67 × 5
#>    variable                glycan_composition glycan_structure lacnac sia_lacnac
#>    <glue>                  <comp>             <struct>         <lgl>  <lgl>     
#>  1 Man(3)GlcNAc(3)         Man(3)GlcNAc(3)    GlcNAc(?1-?)Man… FALSE  FALSE     
#>  2 Man(3)GlcNAc(7)         Man(3)GlcNAc(7)    GlcNAc(?1-?)[Gl… FALSE  FALSE     
#>  3 Man(5)GlcNAc(2)         Man(5)GlcNAc(2)    Man(?1-?)[Man(?… FALSE  FALSE     
#>  4 Man(4)Gal(2)GlcNAc(4)N… Man(4)Gal(2)GlcNA… Neu5Ac(?2-?)Gal… TRUE   TRUE      
#>  5 Man(3)Gal(1)GlcNAc(3)   Man(3)Gal(1)GlcNA… Gal(?1-?)GlcNAc… TRUE   FALSE     
#>  6 Man(3)Gal(2)GlcNAc(4)F… Man(3)Gal(2)GlcNA… Gal(?1-?)GlcNAc… TRUE   FALSE     
#>  7 Man(3)GlcNAc(3)Fuc(1)   Man(3)GlcNAc(3)Fu… GlcNAc(?1-?)Man… FALSE  FALSE     
#>  8 Man(3)GlcNAc(4)         Man(3)GlcNAc(4)    GlcNAc(?1-?)Man… FALSE  FALSE     
#>  9 Man(3)Gal(2)GlcNAc(5)N… Man(3)Gal(2)GlcNA… Neu5Ac(?2-?)Gal… TRUE   TRUE      
#> 10 Man(3)Gal(1)GlcNAc(5)F… Man(3)Gal(1)GlcNA… Neu5Ac(?2-?)Gal… TRUE   TRUE      
#> # ℹ 57 more rows

df <- get_var_info(exp)
df |>
  mutate(as_tibble(count_motifs(glycan_structure, motifs)))
#> # A tibble: 67 × 5
#>    variable                glycan_composition glycan_structure lacnac sia_lacnac
#>    <glue>                  <comp>             <struct>          <int>      <int>
#>  1 Man(3)GlcNAc(3)         Man(3)GlcNAc(3)    GlcNAc(?1-?)Man…      0          0
#>  2 Man(3)GlcNAc(7)         Man(3)GlcNAc(7)    GlcNAc(?1-?)[Gl…      0          0
#>  3 Man(5)GlcNAc(2)         Man(5)GlcNAc(2)    Man(?1-?)[Man(?…      0          0
#>  4 Man(4)Gal(2)GlcNAc(4)N… Man(4)Gal(2)GlcNA… Neu5Ac(?2-?)Gal…      2          2
#>  5 Man(3)Gal(1)GlcNAc(3)   Man(3)Gal(1)GlcNA… Gal(?1-?)GlcNAc…      1          0
#>  6 Man(3)Gal(2)GlcNAc(4)F… Man(3)Gal(2)GlcNA… Gal(?1-?)GlcNAc…      2          0
#>  7 Man(3)GlcNAc(3)Fuc(1)   Man(3)GlcNAc(3)Fu… GlcNAc(?1-?)Man…      0          0
#>  8 Man(3)GlcNAc(4)         Man(3)GlcNAc(4)    GlcNAc(?1-?)Man…      0          0
#>  9 Man(3)Gal(2)GlcNAc(5)N… Man(3)Gal(2)GlcNA… Neu5Ac(?2-?)Gal…      2          1
#> 10 Man(3)Gal(1)GlcNAc(5)F… Man(3)Gal(1)GlcNA… Neu5Ac(?2-?)Gal…      1          1
#> # ℹ 57 more rows
```
