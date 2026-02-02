# Changelog

## glymotif (development version)

### Breaking changes

- Remove “exact” alignment type from `alignment` and `alignments`
  arguments. This alignment type was introduced in v0.12.0 but is no
  longer needed. Users who previously relied on “exact” alignment can
  use “substructure” as the closest available alternative, but note that
  it is less strict and may produce more matches, so results may differ.

### New features

- Add a `match_degree` argument to all motif matching functions.
  `match_degree` provides a more delicate way to control the alignment
  of each monosaccharide residue.

## glymotif 0.12.2

### Minor improvements and fixes

- Add validation to raise an error when `motifs` parameter contains
  duplicates in
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  [`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md),
  [`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  and
  [`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md).
  ([\#8](https://github.com/glycoverse/glymotif/issues/8))
- Update all motif functions to properly preserve and use names from
  input vectors, taking advantage of glyrepr 0.10.0’s support for names
  on glycan structure vectors.
  ([\#5](https://github.com/glycoverse/glymotif/issues/5))([\#6](https://github.com/glycoverse/glymotif/issues/6))

## glymotif 0.12.1

### Minor improvements and fixes

- Adapt to glyrepr 0.10.0.

## glymotif 0.12.0

We introduce the “dynamic motif” feature in this release. Instead of
specifying motifs manually or using motifs from the database, you can
now extract motifs existing in your set of glycans. Two functions,
`extract_motifs()` and `extract_motifs_lgl()`, are added for this
purpose.

### New features

- Add
  [`extract_motif()`](https://glycoverse.github.io/glymotif/reference/extract_motif.md)
  to find all unique substructures (motifs) in a set of glycans. This
  function is suitable for motif finding in small glycans like O-GalNAc
  glycans.
- Add
  [`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)
  to find all branching motifs in a set of N-glycans. This function is
  particular useful for N-glycan motif finding, where antennary patterns
  can be versatile.
- The `alignment` or `alignments` argument of all related functions now
  supports a new alignment type: “exact”. This type of alignment works
  the best with
  [`extract_branch_motif()`](https://glycoverse.github.io/glymotif/reference/extract_branch_motif.md)
  results.

## glymotif 0.11.2

### Minor improvements and fixes

- glymotif now depends on the CRAN version of glyparse.

## glymotif 0.11.1

### Minor improvements and fixes

- glymotif now depends on the CRAN version of glyrepr.

## glymotif 0.11.0

### New features

- [`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
  now support a new argument `strict_sub` to control the strictness of
  substituent matching.The argument defaults to `TRUE`, to be consistent
  with the previous behavior. When `strict_sub` is `FALSE`, the
  substituent is optional in the motif, so the glycan “Neu5Ac9Ac(a2-”
  can match the motif “Neu5Ac(a2-”.

## glymotif 0.10.0

### New features

- [`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  and
  [`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  now support data frames as input.

### Minor improvements and bug fixes

- Fix a bug introduced by igraph v2.2.0.
- Performance optimization for the core algorithm.

## glymotif 0.9.1

### Minor improvements and bug fixes

- Fix bugs introduced by the breaking changes in `glyexp` 0.10.0.

## glymotif 0.9.0

### Breaking changes

- Remove `quantify_motifs()`. This function is reimplemented in the
  `glydet` package, with more features and better performance.

## glymotif 0.8.1

### Minor improvements and bug fixes

- `quantify_motifs()` now returns a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with “traitomics” type for glycomics data, and
  “traitproteomics” type for glycoproteomics data, instead of
  “traitomics” for all input.

## glymotif 0.8.0

### New features

- `quantify_motifs()` has been totally rewritten. Serious bug about
  column aggregation with descriptional columns about glycan structures
  are fixed. The function now behaves like `glydet::derive_traits()`,
  adding back columns in `var_info` only when they have “many-to-one”
  relationship with glycosites (unique combinations of `protein` and
  `protein_site`).
- `quantify_motifs()`,
  [`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md),
  and
  [`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  now support a character vector ofglycan structure strings as the
  “glycan_structure” column in `var_info`.

### Minor improvements and bug fixes

- Add a section in the Get Started vignette about ambiguity handling.
- Update the “Working with glyexp” vignette to use
  [`glyexp::real_experiment`](https://glycoverse.github.io/glyexp/reference/real_experiment.html).
- Update the URL of GlycoMotif in the documentation of
  [`all_motifs()`](https://glycoverse.github.io/glymotif/reference/all_motifs.md).
- Remove old “N-glycans” vignette from README.

## glymotif 0.7.0

### Breaking changes

- Remove `add_comp_descriptions()`, `add_glycan_descriptions()`,
  `add_struct_descriptions()`, `describe_n_glycans()`,
  `has_bisecting()`, `is_n_glycan()`, `n_antennae()`, `n_arm_fuc()`,
  `n_core_fuc()`, `n_gal()`, `n_glycan_type()`, `n_terminal_gal()`.
  These functions are replaced by functions in the `glydet` package now.

### Minor improvements and bug fixes

- Update dependencies to depend on release versions of glycoverse
  packages.
- `quantify_motifs()` now returns a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with “traitomics” type instead of “motifomics”.

## glymotif 0.6.2

### Minor improvements and bug fixes

- Fix bugs introduced by the breaking changes in `glyrepr` v0.7.0.

## glymotif 0.6.1

### Minor improvements and bug fixes

- Fix some incorrect structure strings in tests and vignettes.

## glymotif 0.6.0

### Breaking changes

- `available_motifs()` is renamed to
  [`all_motifs()`](https://glycoverse.github.io/glymotif/reference/all_motifs.md).

### New features

- Add
  [`match_motif()`](https://glycoverse.github.io/glymotif/reference/match_motif.md)
  and
  [`match_motifs()`](https://glycoverse.github.io/glymotif/reference/match_motif.md)
  to find all occurrences of the given motif(s) in the glycans.
  Node-to-node mapping is returned for each match.

## glymotif 0.5.0

### Breaking changes

- Remove support for omitted reducing-end anomers in IUPAC-condensed
  strings. Strings like “Gal(b1-3)GlcNAc” are no longer valid. You must
  specify the anomer of the reducing-end monosaccharide,
  e.g. “Gal(b1-3)GlcNAc(b1-”.
- N-glycan functions now raise a warning when the input glycans are not
  N-glycans and return `NA` for those glycans, instead of throwing an
  error.

### New features

- Add support for multiple glycan structure string formats in all
  functions. This includes IUPAC-condensed, IUPAC-short, IUPAC-extended,
  GlycoCT, WURCS, pGlyco-style, and StrucGP-style.
- All N-glycan functions now support vectorization, including
  `is_n_glycan()`, `n_glycan_type()`, `has_bisecting()`, `n_antennae()`,
  `n_core_fuc()`, `n_arm_fuc()`, `n_gal()`, and `n_terminal_gal()`.

### Minor improvements and bug fixes

- Update vignettes to explicitly include reducing-end anomers in
  IUPAC-condensed strings.
- Fix a bug that pausimannose-type glycans are not correctly recognized
  in `describe_n_glycans()`.

## glymotif 0.4.4

### Minor improvements and bug fixes

- Update vignette “Working with glyexp” to reflect the changes in
  `glyread` v0.5.0.

## glymotif 0.4.3

### Minor improvements and bug fixes

- Huge (really huge) performance optimization: all motif matching
  functions in this package now speed up 1000x, thanks to the fix of a
  performance bug about monosaccharide type conversion.

## glymotif 0.4.2

### Minor improvements

- Better error messages for invalid input.

## glymotif 0.4.1

### Bug fixes

- Fix a major bug in `quantify_motifs()`.

## glymotif 0.4.0

### Major changes

- Add `quantify_motifs()` to quantify motifs in a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html).

### Minor improvements

- Split the “Get Started” vignette into three:
  - The new “Get Started” vignette, only covering basic motif matching
    functions.
  - A new “N-Glycans” vignette, covering N-glycan-specific functions.
  - A new “Working with
    [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)”
    vignette, covering the integration with
    [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html).

## glymotif 0.3.1

### Minor improvements

- Update the documentation to reflect the naming rules of the return
  values of
  [`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md).

## glymotif 0.3.0

### Major changes

- Add `add_glycan_descriptions()`, `add_struct_descriptions()`, and
  `add_comp_descriptions()`. These functions are moved from `glyexp` to
  `glymotif` for better separation of concerns.
- Add
  [`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  and
  [`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  for adding motif annotations to a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html).

### Minor improvements

- Remove the `parallel` argument from `add_glycan_descriptions()`,
  `add_struct_descriptions()`, and `add_comp_descriptions()`. We found
  it not very useful on a regular basis.
- Update the vignette to reflect the new functions above.

### Bug fixes

- Fix a bug in
  [`get_motif_structure()`](https://glycoverse.github.io/glymotif/reference/get_motif_structure.md),
  where the order of the results was not consistent with the order of
  the input motifs. This caused unexpected results in a wide range of
  functions including
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and the new
  [`add_motifs_lgl()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  and
  [`add_motifs_int()`](https://glycoverse.github.io/glymotif/reference/add_motifs_int.md)
  functions, when using motif names as input. Alone with
  [`get_motif_structure()`](https://glycoverse.github.io/glymotif/reference/get_motif_structure.md),
  similar bugs in
  [`get_motif_alignment()`](https://glycoverse.github.io/glymotif/reference/get_motif_structure.md)
  and
  [`get_motif_aglycon()`](https://glycoverse.github.io/glymotif/reference/get_motif_structure.md)
  are also fixed.

## glymotif 0.2.1

### Minor improvements

- [`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md)
  now support multiple substituents in the motif, to align with the
  updates in `glyrepr` v0.5.0.

## glymotif 0.2.0

### Bug fixes

- Fixed monosaccharide type matching logic in
  [`have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.md),
  [`count_motif()`](https://glycoverse.github.io/glymotif/reference/count_motif.md),
  and
  [`count_motifs()`](https://glycoverse.github.io/glymotif/reference/count_motif.md).
  - Generic glycans can no longer cause errors when compared with
    concrete motifs
  - Concrete glycans now properly match generic motifs (converted to
    generic first)
  - Generic glycans correctly return FALSE when compared with concrete
    motifs
  - This resolves incorrect results in mixed-type motif matrices

### Breaking changes

- Note: While technically backward compatible in API, the behavior
  change in monosaccharide type matching may affect code that depended
  on the previous (incorrect) error-throwing or wrong-result behavior.
