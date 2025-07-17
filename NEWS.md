# glymotif 0.4.2

## Minor improvements

* Better error messages for invalid input.

# glymotif 0.4.1

## Bug fixes

* Fix a major bug in `quantify_motifs()`.

# glymotif 0.4.0

## Major changes

* Add `quantify_motifs()` to quantify motifs in a `glyexp::experiment()`.

## Minor improvements

* Split the "Get Started" vignette into three:
  * The new "Get Started" vignette, only covering basic motif matching functions.
  * A new "N-Glycans" vignette, covering N-glycan-specific functions.
  * A new "Working with `glyexp::experiment()`" vignette, covering the integration with `glyexp::experiment()`.

# glymotif 0.3.1

## Minor improvements

* Update the documentation to reflect the naming rules of the return values 
  of `have_motif()`, `have_motifs()`, `count_motif()`, and `count_motifs()`.

# glymotif 0.3.0

## Major changes

* Add `add_glycan_descriptions()`, `add_struct_descriptions()`, and `add_comp_descriptions()`.
  These functions are moved from `glyexp` to `glymotif` for better separation of concerns.
* Add `add_motifs_int()` and `add_motifs_lgl()` for adding motif annotations to a `glyexp::experiment()`.

## Minor improvements

* Remove the `parallel` argument from `add_glycan_descriptions()`, `add_struct_descriptions()`, and `add_comp_descriptions()`. 
  We found it not very useful on a regular basis.
* Update the vignette to reflect the new functions above.

## Bug fixes

* Fix a bug in `get_motif_structure()`, where the order of the results was not consistent with
  the order of the input motifs.
  This caused unexpected results in a wide range of functions including `have_motifs()`, `count_motifs()`,
  and the new `add_motifs_lgl()` and `add_motifs_int()` functions,
  when using motif names as input.
  Alone with `get_motif_structure()`, similar bugs in `get_motif_alignment()` and `get_motif_aglycon()`
  are also fixed.

# glymotif 0.2.1

## Minor improvements

* `have_motif()`, `have_motifs()`, `count_motif()`, and `count_motifs()` now
  support multiple substituents in the motif,
  to align with the updates in `glyrepr` v0.5.0.


# glymotif 0.2.0

## Bug fixes

* Fixed monosaccharide type matching logic in `have_motif()`, `have_motifs()`, 
  `count_motif()`, and `count_motifs()`.
  - Generic glycans can no longer cause errors when compared with concrete motifs
  - Concrete glycans now properly match generic motifs (converted to generic first)
  - Generic glycans correctly return FALSE when compared with concrete motifs
  - This resolves incorrect results in mixed-type motif matrices

## Breaking changes

* Note: While technically backward compatible in API, the behavior change in 
  monosaccharide type matching may affect code that depended on the previous
  (incorrect) error-throwing or wrong-result behavior.
