
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glymotif <a href="https://glycoverse.github.io/glymotif/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glymotif)](https://CRAN.R-project.org/package=glymotif)
[![R-CMD-check](https://github.com/glycoverse/glymotif/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glymotif/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glymotif/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glymotif)
<!-- badges: end -->

The goal of glymotif is to extract glycan motifs from glycan structures.
It works seemlessly with [glyrepr](https://github.com/glycoverse/glyrepr)
and [glyparse](https://github.com/glycoverse/glyparse).

## Installation

You can install the development version of glymotif from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("glycoverse/glymotif")
```

## Example

``` r
library(glymotif)
library(glyparse)
```

Say we have a glycan, …

``` r
(glycan <- parse_iupac_condensed("Gal(b1-3)GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"))
#> Glycan Graph (NE)
#> Gal: 2, GalNAc: 1, GlcNAc: 1
#> ------------------
#> GalNAc (a1-)
#> └─Gal (b1-3)
#>   └─GlcNAc (b1-3)
#>     └─Gal (b1-3)
```

… and we want to check if it has the O-Glycan core 1 motif.

``` r
has_motif(glycan, "Gal(b1-3)GalNAc(a1-", alignment = "core")
#> [1] TRUE
```

Or use the motif name directly.

``` r
has_motif(glycan, "O-Glycan core 1")
#> [1] TRUE
```
