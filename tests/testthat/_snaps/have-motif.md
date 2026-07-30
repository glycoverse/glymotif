# all-generic glycans and all-concrete motifs warn

    Matching generic `glycans` against concrete `motifs` always returns no matches in strict mode.
    i Convert motifs with `glyrepr::convert_to_generic()`, or use `mode = "lenient"` when generic identities should match their concrete counterparts.

---

    Matching generic `glycans` against concrete `motifs` always returns no matches in strict mode.
    i Convert motifs with `glyrepr::convert_to_generic()`, or use `mode = "lenient"` when generic identities should match their concrete counterparts.

---

    Matching generic `glycans` against concrete `motifs` always returns no matches in strict mode.
    i Convert motifs with `glyrepr::convert_to_generic()`, or use `mode = "lenient"` when generic identities should match their concrete counterparts.

# all-topological glycans and higher-level motifs warn

    Matching lower-level `glycans` against higher-level `motifs` usually returns no matches.
    i `glycans` have "topological" structure level, while `motifs` have "intact" structure level.
    i Use motifs at the same structure level as the glycans, or remove motif linkage constraints with `glyrepr::remove_linkages()`.
    i See `?get_structure_level` for details.

---

    Matching lower-level `glycans` against higher-level `motifs` usually returns no matches.
    i `glycans` have "topological" structure level, while `motifs` have "intact/partial" structure level.
    i Use motifs at the same structure level as the glycans, or remove motif linkage constraints with `glyrepr::remove_linkages()`.
    i See `?get_structure_level` for details.

# have_motif: generic glycan does not match concrete motif

    Matching generic `glycans` against concrete `motifs` always returns no matches in strict mode.
    i Convert motifs with `glyrepr::convert_to_generic()`, or use `mode = "lenient"` when generic identities should match their concrete counterparts.

# have_motif: complex structures with type conversion

    Matching generic `glycans` against concrete `motifs` always returns no matches in strict mode.
    i Convert motifs with `glyrepr::convert_to_generic()`, or use `mode = "lenient"` when generic identities should match their concrete counterparts.

# floating motifs are rejected clearly

    Code
      have_motif("Gal(a1-3)Man(a1-", "{Gal(a1-3)}Man(a1-3)GlcNAc(b1-")
    Condition
      Error in `have_motif()`:
      ! `motifs` cannot contain unresolved floating parts.
      i Localize floating motif parts before matching.
