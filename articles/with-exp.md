# Working with glyexp

## When Motifs Meet Experiments: A Perfect Partnership 🤝

The real power of `glymotif` shines brightest when it joins forces with
other tools in the `glycoverse` ecosystem. If you’re already using
`glyread` to import your glycoproteomics results and `glyexp` to manage
your experimental data, you’re in for a treat! `glymotif` provides some
incredibly useful functions to perform `experiment` manipulation about
motifs—think of it as adding a new lens to your analytical microscope.
🔬

**Important note:** This vignette assumes you’re familiar with the
`glyexp` package. If you haven’t met it yet, we highly recommend
checking out its
[introduction](https://glycoverse.github.io/glyexp/articles/glyexp.html)
first. Trust us—it’s worth the detour! 🚀

``` r
library(glymotif)
library(glyexp)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:glyexp':
#> 
#>     select_var
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Meet Our Star Player: Real Glycoproteomics Data 🌟

Time to roll up our sleeves and dive into the good stuff! 💪 Let’s work
with a real-world dataset that will showcase what `glymotif` can do when
it teams up with actual experimental data. We will use `real_experiment`
from the `glyexp` package, an serum N-glycoproteomics study with 12
samples. Firstly, let’s use
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
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [152ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [47ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values✔ Imputing missing values [25ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [1s]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [15ms]
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3880 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
```

Now, let’s peek under the hood and see what treasures we’re working
with! 👀

``` r
get_var_info(exp)
#> # A tibble: 3,880 × 6
#>    variable protein glycan_composition      glycan_structure  protein_site gene 
#>    <chr>    <chr>   <comp>                  <struct>                 <int> <chr>
#>  1 V1       P08185  Hex(5)HexNAc(4)NeuAc(2) NeuAc(??-?)Hex(?…          176 SERP…
#>  2 V2       P04196  Hex(5)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#>  3 V3       P04196  Hex(5)HexNAc(4)         Hex(??-?)HexNAc(…          344 HRG  
#>  4 V4       P04196  Hex(5)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#>  5 V5       P10909  Hex(6)HexNAc(5)         Hex(??-?)HexNAc(…          291 CLU  
#>  6 V6       P04196  Hex(5)HexNAc(4)NeuAc(2) NeuAc(??-?)Hex(?…          344 HRG  
#>  7 V7       P04196  Hex(5)HexNAc(4)         Hex(??-?)HexNAc(…          345 HRG  
#>  8 V8       P04196  Hex(5)HexNAc(4)dHex(2)  dHex(??-?)Hex(??…          344 HRG  
#>  9 V9       P04196  Hex(4)HexNAc(3)         Hex(??-?)HexNAc(…          344 HRG  
#> 10 V10      P04196  Hex(4)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#> # ℹ 3,870 more rows
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

What we have here is a beautiful N-glycoproteomics dataset featuring 500
PSMs (Peptide Spectrum Matches) across 12 samples — a perfect playground
for motif analysis! 🎮

**Pro tip:** 💡 In real-world data analysis, you’ll definitely want to
use `glyclean` to perform data preprocessing before diving into any
analysis. Think of it as washing your vegetables before
cooking—essential for the best results! 🥕

## Adding Motif Annotations to an Experiment 🏷️

Here’s where things get interesting! 🤔 We know that the variable
information tibble contains all the juicy details about each glycoform -
the proteins, sites, and glycan structures. But what if we want to
sprinkle some motif magic into the mix? What if we want to add more
columns that tell us about the motifs hiding in our glycans?

**The simple approach:** One motif at a time (very intuitive!) 🎯

``` r
exp |> 
  mutate_var(n_hex = have_motif(glycan_structure, "Hex(a1-")) |>
  get_var_info() |>
  select(variable, protein, glycan_structure, n_hex)
#> # A tibble: 3,880 × 4
#>    variable protein glycan_structure                                       n_hex
#>    <chr>    <chr>   <struct>                                               <lgl>
#>  1 V1       P08185  NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)… FALSE
#>  2 V2       P04196  NeuAc(??-?)Hex(??-?)HexNAc(??-?)[HexNAc(??-?)]Hex(??-… FALSE
#>  3 V3       P04196  Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)H… FALSE
#>  4 V4       P04196  NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)He… FALSE
#>  5 V5       P10909  Hex(??-?)HexNAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[H… FALSE
#>  6 V6       P04196  NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)… FALSE
#>  7 V7       P04196  Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)H… FALSE
#>  8 V8       P04196  dHex(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[dHex(??-?)He… FALSE
#>  9 V9       P04196  Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)Hex… FALSE
#> 10 V10      P04196  NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[HexNAc(??-?… FALSE
#> # ℹ 3,870 more rows
```

**The tempting approach:** Multiple motifs at once (you might be tempted
to do this…) 🤷‍♀️

``` r
# Don't do this
exp |>
  mutate_var(
    n_hex = have_motif(glycan_structure, "Hex(a1-"),
    n_hexna = have_motif(glycan_structure, "HexNAc(a1-"),
    n_dhex = have_motif(glycan_structure, "dHex(a1-")
  )
```

Hold on there, speed racer! 🛑 While this approach works, it’s not very
efficient, and your computer won’t thank you for it. Here’s why: those
three separate calls to
[`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md)
all perform time-consuming validations and conversions on the same set
of glycan structures. It’s like washing the same dishes three times
instead of doing them all at once! 🍽️ Plus, there’s a lot of repetitive
typing. You have to type `have_motif` and `glycan_structure` three times
— talk about finger fatigue! 😴

**The smart approach:** Use the
[`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
or
[`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
functions instead! ⚡ They might look like simple syntactic sugar, but
they’re actually optimized powerhouses designed specifically for this
exact scenario:

``` r
exp2 <- exp |> 
  add_motifs_lgl(c(motif1 = "Hex(??-", motif2 = "HexNAc(??-", motif3 = "dHex(??-"))
```

Voilà! 🎉 The motif annotations are now seamlessly integrated into your
variable information tibble.

``` r
exp2 |>
  get_var_info() |>
  select(variable, motif1, motif2, motif3)
#> # A tibble: 3,880 × 4
#>    variable motif1 motif2 motif3
#>    <chr>    <lgl>  <lgl>  <lgl> 
#>  1 V1       TRUE   TRUE   FALSE 
#>  2 V2       TRUE   TRUE   FALSE 
#>  3 V3       TRUE   TRUE   FALSE 
#>  4 V4       TRUE   TRUE   FALSE 
#>  5 V5       TRUE   TRUE   FALSE 
#>  6 V6       TRUE   TRUE   FALSE 
#>  7 V7       TRUE   TRUE   FALSE 
#>  8 V8       TRUE   TRUE   TRUE  
#>  9 V9       TRUE   TRUE   FALSE 
#> 10 V10      TRUE   TRUE   FALSE 
#> # ℹ 3,870 more rows
```

But wait, what can you actually *do* with these shiny new columns? 🤔
The possibilities are endless, but here’s a tantalizing example to get
your creative juices flowing:

``` r
# You can perform pathway enrichment on all glycoproteins containing some motif:
exp2 |>
  filter_var(motif1 == TRUE) |>
  gly_enrich_reactome()  # from the `glystats` package
```

## The Art of Motif Quantification in Experiments 📊

Want to quantify the motifs in your experiment? Try the `glydet`
package! It provides the `quantify_motifs()` function to perform
relative and absolute motif quantification.
