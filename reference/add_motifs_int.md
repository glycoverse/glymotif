# Add Motif Annotations

This function adds motif annotations to the variable information of a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
or a tibble with a structure column. `add_motifs_int()` adds integer
annotations (how many motifs are present). `add_motifs_lgl()` adds
boolean annotations (whether the motif is present).

## Usage

``` r
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)

add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)

# S3 method for class 'glyexp_experiment'
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)

# S3 method for class 'glyexp_experiment'
add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)

# S3 method for class 'data.frame'
add_motifs_int(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)

# S3 method for class 'data.frame'
add_motifs_lgl(
  x,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object, or a tibble with a structure column.

- motifs:

  A character vector of motif names, glycan structure strings, or a
  'glyrepr_structure' object. For glycan structure strings, all formats
  supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)
  are accepted, including IUPAC-condensed, WURCS, GlycoCT, and others.

- alignments:

  A character vector specifying alignment types for each motif. Can be a
  single value (applied to all motifs) or a vector of the same length as
  motifs.

- ignore_linkages:

  A logical value. If `TRUE`, linkages will be ignored in the
  comparison.

- strict_sub:

  A logical value. If `TRUE` (default), substituents will be matched in
  strict mode, which means if the glycan has a substituent in some
  residue, the motif must have the same substituent to be matched.

- ...:

  Additional arguments passed to the method.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with motif annotations added to the variable information.

## About Names

The naming rule for the new columns is similar to that of
[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md).
Briefly, you can use named character vector to name the motifs, and that
will be used as the new column names. The only catchup is that you
cannot pass a named
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
to `motifs`. This is a fundamental limitation of the `vctrs_rcrd` class,
which
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
is built on.

## Why do we need these functions

Adding one motif annotation to a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
is easy:

    exp |>
      mutate_var(has_hex = have_motif(glycan_structure, "Hex"))

However, adding multiple motifs is not as straightforward. You can still
use `mutate_var()` to add multiple motifs like this:

    exp |>
      mutate_var(
        n_hex = count_motif(glycan_structure, "Hex"),
        n_dhex = count_motif(glycan_structure, "dHex"),
        n_hexnac = count_motif(glycan_structure, "HexNAc"),
      )

This method has two problems:

1.  it has a lot of boilerplate code (a lot of typing)

2.  it is not very efficient, as each call to `count_motif` performs
    validation and conversion on `glycan_structure`, which is a
    time-consuming process.

Therefore, we think it would be better to have a function that adds
multiple motif annotations in a single call, in a more intuitive way.
That's why we provide these two functions.

Under the hood, they use a more straightforward approach for
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects:

1.  get the motif annotation matrix using
    [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
    or
    [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)

2.  convert the matrix to a tibble

3.  use
    [`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html)
    to add the tibble to the variable information

## See also

[`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
[`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
