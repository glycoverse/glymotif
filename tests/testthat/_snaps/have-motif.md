# wrong glycan types

    Code
      have_motif(glycan, motif)
    Condition
      Error in `valid_glycans_arg()`:
      ! `glycans` must be a 'glyrepr_structure' object or an IUPAC-condensed structure character.

# wrong motif types

    Code
      have_motif(glycan, motif)
    Condition
      Error in `valid_motif_arg()`:
      ! `motif` must be either a 'glyrepr_structure' object with length 1, an IUPAC-condensed structure character scalar, or a known motif name.

# unkown motif name used as input

    Code
      have_motif(glycan, motif)
    Condition
      Error in `value[[3L]]()`:
      ! `motif` must be either a 'glyrepr_structure' object with length 1, an IUPAC-condensed structure character scalar, or a known motif name.

# bad glycan IUPAC

    Code
      have_motif(glycan, motif)
    Condition
      Error in `value[[3L]]()`:
      ! `glycan` could not be parsed as a valid IUPAC-condensed structure.

# warning when user-provided alignment is different from database

    Code
      have_motif(glycan, motif, alignment = "terminal")
    Condition
      Warning:
      The provided alignment type "terminal" is different from the motif's alignment type "core" in database.
    Output
      [1] FALSE

# generic glycan and concrete motif

    Code
      have_motif(glycan, motif)
    Condition
      Error in `ensure_glycans_mono_type()`:
      ! `generic` glycans cannot be compared with `concrete` motifs.

# custom alignment is different from database

    Code
      have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc", "O-Glycan core 1", alignment = "substructure")
    Condition
      Warning:
      The provided alignment type "substructure" is different from the motif's alignment type "core" in database.
    Output
      [1] TRUE

