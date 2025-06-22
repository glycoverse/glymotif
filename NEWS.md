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