# View motif matches on a glycan

Visualize where a motif matches a glycan structure.

## Usage

``` r
view_motif(
  glycan,
  motif,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
)
```

## Arguments

- glycan:

  One of:

  - A scalar
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html).

  - A glycan structure string. All formats supported by
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
    [`db_motif_info()`](https://glycoverse.github.io/glymotif/dev/reference/db_motif_info.md)
    to see all available motifs).

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
  [`have_motif()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
  [`count_motif()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
  and
  [`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md),
  this must be a logical vector with length 1 or the number of motif
  nodes (length 1 is recycled). For
  [`have_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md),
  [`count_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md),
  and
  [`match_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md),
  this must be a list of logical vectors with length equal to `motifs`;
  each element follows the same length rules. When `match_degree` is
  provided, `alignment` and `alignments` are silently ignored.

## Value

A `ggplot` object returned by
[`glydraw::draw_cartoon()`](https://glycoverse.github.io/glydraw/reference/draw_cartoon.html).
If no match is found, the glycan is drawn without highlighted residues
and a cli alert is emitted.

## Details

`view_motif()` matches one `motif` against one `glycan` with the same
matching rules used by
[`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md),
then draws the glycan with the matched residues highlighted.

## See also

[`match_motif()`](https://glycoverse.github.io/glymotif/dev/reference/match_motif.md),
[`glydraw::draw_cartoon()`](https://glycoverse.github.io/glydraw/reference/draw_cartoon.html)

## Examples

``` r
library(glyparse)
library(glyrepr)

glycan <- n_glycan_core()
motif <- parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-")

if (FALSE) { # \dontrun{
view_motif(glycan, motif)
} # }
```
