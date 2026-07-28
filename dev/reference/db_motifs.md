# Get All Motifs from the Database

This function returns a database motif specification. We use GlycoMotif
collections (https://glycomotif.glyomics.org/glycomotif/GlycoMotif) as
the source of the motifs. This function is useful to be integrated with
[`have_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md)
and
[`count_motifs()`](https://glycoverse.github.io/glymotif/dev/reference/count_motif.md).
For example, use `have_motifs(glycans, db_motifs())` to check against
the default GlyGen motif collection, or pass `source_id` to use another
collection.

## Usage

``` r
db_motifs(source_id = "GGM")
```

## Arguments

- source_id:

  A character vector of motif collection identifiers to use. Defaults to
  `"GGM"` for backward compatibility. Use
  `dplyr::distinct(db_motif_info(), source_id, source)` to get all
  available sources. You can use more than one motif collections like
  `c("GGM", "CCRC")`. To use all available motifs, use the `"GM"`
  collection directly.

## Value

A `db_motifs_spec` object.

## Details

Use
[`db_motif_info()`](https://glycoverse.github.io/glymotif/dev/reference/db_motif_info.md)
to inspect the motifs included in the database. You can use
`dplyr::distinct(db_motif_info(), source_id, source)` to get all
available sources.

## Data source and license

The bundled annotations are derived from the
[GlycoMotif](https://glycomotif.glyomics.org/glycomotif/GlycoMotif)
resource. GlyGen distributes its database sets under the [Creative
Commons Attribution 4.0 International
license](https://www.glygen.org/license/). See the package `COPYRIGHTS`
file for attribution and snapshot details.

## Examples

``` r
db_motifs()
#> <<db_motifs_spec>>
#> This object should be passed to the `motifs` argument of `have_motifs()`,
#> `count_motifs()`, or `match_motifs()`.
#> Configuration: uses GlycoMotif database source ID: "GGM"
```
