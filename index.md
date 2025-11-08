# glymotif

Provides comprehensive tools for glycan motif analysis and detection in
glycobioinformatics research. The package enables users to identify,
count, and match glycan motifs (recurring substructures) within complex
glycan structures using advanced subgraph isomorphism algorithms. It
includes a curated database of known motifs from the GlycoMotif GlyGen
Collection, supports both concrete and generic monosaccharide matching,
and offers flexible alignment options (core, terminal, or anywhere). Key
functionalities include motif presence detection, occurrence counting,
detailed node-to-node mapping, and batch analysis of multiple glycans
against multiple motifs. The package seamlessly integrates with the
glycoverse ecosystem, particularly ‘glyrepr’ and ‘glyparse’, making it
essential for structural glycomics analysis, biomarker discovery, and
understanding glycan-mediated biological processes.

## Installation

You can install the latest release of glymotif from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glycoverse/glymotif@*release")
```

Or install the development version:

``` r
remotes::install_github("glycoverse/glymotif")
```

## Documentation

- 🚀 Get started:
  [Here](https://glycoverse.github.io/glymotif/articles/glymotif.html)
- 🔧 Motif matching rules:
  [Here](https://glycoverse.github.io/glymotif/articles/motif-matching.html)
- 🔬 Working with
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html):
  [Here](https://glycoverse.github.io/glymotif/articles/with-exp.html)
- 📚 Reference:
  [Here](https://glycoverse.github.io/glymotif/reference/index.html)

## Role in `glycoverse`

`glymotif` provides possibilities for one important job in
glyco-bioinformatics: to detect motifs in glycans. The package is
designed to be used directly by users for structural analysis, as well
as to provide backend support for other packages in the `glycoverse`
ecosystem.

## Example

``` r
library(glymotif)
library(glyparse)
```

Say we have a glycan, …

``` r
(glycan <- parse_iupac_condensed("Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"))
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```

… and we want to check if it has the O-Glycan core 1 motif.

``` r
have_motif(glycan, "Gal(b1-3)GalNAc(a1-", alignment = "core")
#> [1] TRUE
```

Or use the motif name directly.

``` r
have_motif(glycan, "O-Glycan core 1")
#> [1] TRUE
```
