# Working with glyexp

## Working with Motifs in Experiments

`glymotif` is most useful when combined with other tools in the
`glycoverse` ecosystem. If you’re already using `glyread` to import your
glycoproteomics results and `glyexp` to manage your experimental data,
you can add motif-level information directly to your experiment.
`glymotif` provides functions for adding motif information to
`experiment` objects.

**Important note:** This vignette assumes you’re familiar with the
`glyexp` package. If you haven’t used it yet, start with its
[introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html).

``` r

library(glymotif)
library(glyexp)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Real Glycoproteomics Data

Let’s work with a real-world dataset that shows how `glymotif` can be
used with experimental data. We will use `real_experiment` from the
`glyexp` package, a serum N-glycoproteomics study with 12 samples.
First, let’s use
[`glyclean::auto_clean()`](https://glycoverse.github.io/glyclean/reference/auto_clean.html)
to preprocess the data.

``` r

library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate

exp <- auto_clean(real_experiment)
#> 
#> ── Normalizing data ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ No QC samples found. Using all samples.
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 24 (0.56%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ No QC samples found. Using default imputation method based on sample size.
#> ℹ Sample size <= 30, using `impute_sample_min()`.
#> ✔ Imputation completed.
#> 
#> ── Aggregating data ──
#> 
#> ℹ Aggregating to "gfs" level
#> ✔ Aggregation completed.
#> 
#> ── Normalizing data again ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column  not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3979 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
```

Now, let’s inspect the variable and sample information.

``` r

get_var_info(exp)
#> # A tibble: 3,979 × 6
#>    variable       protein glycan_composition glycan_structure protein_site gene 
#>    <chr>          <chr>   <comp>             <struct>                <int> <chr>
#>  1 P08185-176-He… P08185  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          176 SERP…
#>  2 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  3 P04196-344-He… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          344 HRG  
#>  4 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  5 P10909-291-He… P10909  Hex(6)HexNAc(5)    Hex(??-?)HexNAc…          291 CLU  
#>  6 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  7 P04196-345-He… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          345 HRG  
#>  8 P04196-344-He… P04196  Hex(5)HexNAc(4)dH… dHex(??-?)Hex(?…          344 HRG  
#>  9 P04196-344-He… P04196  Hex(4)HexNAc(3)    Hex(??-?)HexNAc…          344 HRG  
#> 10 P04196-344-He… P04196  Hex(4)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#> # ℹ 3,969 more rows
```

``` r

get_sample_info(exp)
#> # A tibble: 12 × 2
#>    sample group
#>    <chr>  <fct>
#>  1 C1     C    
#>  2 C2     C    
#>  3 C3     C    
#>  4 H1     H    
#>  5 H2     H    
#>  6 H3     H    
#>  7 M1     M    
#>  8 M2     M    
#>  9 M3     M    
#> 10 Y1     Y    
#> 11 Y2     Y    
#> 12 Y3     Y
```

This is an N-glycoproteomics dataset with 500 PSMs (Peptide Spectrum
Matches) across 12 samples, which is a useful setting for motif
analysis.

**Tip:** In real-world data analysis, you will usually want to use
`glyclean` to perform data preprocessing before diving into any
analysis.

## Adding Motif Annotations to an Experiment

The variable information tibble contains details about each glycoform -
the proteins, sites, and glycan structures. We can add columns that
describe whether specific motifs are present in each glycan.

**The simple approach:** One motif at a time

``` r

exp |> 
  mutate_var(n_hex = have_motif(glycan_structure, "Hex(a1-")) |>
  get_var_info() |>
  select(variable, protein, glycan_structure, n_hex)
#> # A tibble: 3,979 × 4
#>    variable                             protein glycan_structure           n_hex
#>    <chr>                                <chr>   <struct>                   <lgl>
#>  1 P08185-176-Hex(5)HexNAc(4)NeuAc(2)   P08185  NeuAc(??-?)Hex(??-?)HexNA… FALSE
#>  2 P04196-344-Hex(5)HexNAc(4)NeuAc(1)-1 P04196  NeuAc(??-?)Hex(??-?)HexNA… FALSE
#>  3 P04196-344-Hex(5)HexNAc(4)           P04196  Hex(??-?)HexNAc(??-?)Hex(… FALSE
#>  4 P04196-344-Hex(5)HexNAc(4)NeuAc(1)-2 P04196  NeuAc(??-?)Hex(??-?)HexNA… FALSE
#>  5 P10909-291-Hex(6)HexNAc(5)-1         P10909  Hex(??-?)HexNAc(??-?)Hex(… FALSE
#>  6 P04196-344-Hex(5)HexNAc(4)NeuAc(2)   P04196  NeuAc(??-?)Hex(??-?)HexNA… FALSE
#>  7 P04196-345-Hex(5)HexNAc(4)           P04196  Hex(??-?)HexNAc(??-?)Hex(… FALSE
#>  8 P04196-344-Hex(5)HexNAc(4)dHex(2)    P04196  dHex(??-?)Hex(??-?)HexNAc… FALSE
#>  9 P04196-344-Hex(4)HexNAc(3)-1         P04196  Hex(??-?)HexNAc(??-?)Hex(… FALSE
#> 10 P04196-344-Hex(4)HexNAc(4)NeuAc(1)   P04196  NeuAc(??-?)Hex(??-?)HexNA… FALSE
#> # ℹ 3,969 more rows
```

**A repetitive approach:** Multiple motifs one call at a time

``` r

# Don't do this
exp |>
  mutate_var(
    n_hex = have_motif(glycan_structure, "Hex(a1-"),
    n_hexna = have_motif(glycan_structure, "HexNAc(a1-"),
    n_dhex = have_motif(glycan_structure, "dHex(a1-")
  )
```

While this approach works, it’s not very efficient. Those three separate
calls to
[`have_motif()`](https://glycoverse.github.io/glymotif/dev/reference/have_motif.md)
all perform time-consuming validations and conversions on the same set
of glycan structures. Plus, there’s a lot of repetitive typing. You have
to type `have_motif` and `glycan_structure` three times.

**The preferred approach:** Use the
[`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/dev/reference/add_motifs_int.md)
or
[`add_motifs_int()`](https://glycoverse.github.io/glymotif/dev/reference/add_motifs_int.md)
functions instead. They might look like simple syntactic sugar, but they
are optimized for this scenario:

``` r

exp2 <- exp |> 
  add_motifs_lgl(c(motif1 = "Hex(??-", motif2 = "HexNAc(??-", motif3 = "dHex(??-"))
```

The motif annotations are now integrated into your variable information
tibble.

``` r

exp2 |>
  get_var_info() |>
  select(variable, motif1, motif2, motif3)
#> # A tibble: 3,979 × 4
#>    variable                             motif1 motif2 motif3
#>    <chr>                                <lgl>  <lgl>  <lgl> 
#>  1 P08185-176-Hex(5)HexNAc(4)NeuAc(2)   TRUE   TRUE   FALSE 
#>  2 P04196-344-Hex(5)HexNAc(4)NeuAc(1)-1 TRUE   TRUE   FALSE 
#>  3 P04196-344-Hex(5)HexNAc(4)           TRUE   TRUE   FALSE 
#>  4 P04196-344-Hex(5)HexNAc(4)NeuAc(1)-2 TRUE   TRUE   FALSE 
#>  5 P10909-291-Hex(6)HexNAc(5)-1         TRUE   TRUE   FALSE 
#>  6 P04196-344-Hex(5)HexNAc(4)NeuAc(2)   TRUE   TRUE   FALSE 
#>  7 P04196-345-Hex(5)HexNAc(4)           TRUE   TRUE   FALSE 
#>  8 P04196-344-Hex(5)HexNAc(4)dHex(2)    TRUE   TRUE   TRUE  
#>  9 P04196-344-Hex(4)HexNAc(3)-1         TRUE   TRUE   FALSE 
#> 10 P04196-344-Hex(4)HexNAc(4)NeuAc(1)   TRUE   TRUE   FALSE 
#> # ℹ 3,969 more rows
```

You can use these new columns for downstream filtering and analysis. For
example:

``` r

# You can perform pathway enrichment on all glycoproteins containing some motif:
exp2 |>
  filter_var(motif1 == TRUE) |>
  gly_enrich_reactome()  # from the `glystats` package
```

## Motif Quantification in Experiments

Want to quantify the motifs in your experiment? Try the `glydet`
package. It provides the `quantify_motifs()` function to perform
relative and absolute motif quantification.
