# warning when user-provided alignment is different from database

    Code
      has_motif(glycan, motif, alignment = "terminal")
    Condition
      Warning:
      The provided alignment type "terminal" is different from the motif's alignment type "core" in database.
    Output
      [1] FALSE

# simple glycan and concrete motif

    Code
      has_motif(glycan, motif)
    Condition
      Error in `ensure_glycan_mono_type()`:
      ! The monosaccharide type of `glycan` cannot be obscurer than `motif`.
      x "simple" is obscurer than "concrete".

# simple glycan and generic motif

    Code
      has_motif(glycan, motif)
    Condition
      Error in `ensure_glycan_mono_type()`:
      ! The monosaccharide type of `glycan` cannot be obscurer than `motif`.
      x "simple" is obscurer than "generic".

# generic glycan and concrete motif

    Code
      has_motif(glycan, motif)
    Condition
      Error in `ensure_glycan_mono_type()`:
      ! The monosaccharide type of `glycan` cannot be obscurer than `motif`.
      x "generic" is obscurer than "concrete".

# custom alignment is different from database

    Code
      has_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc", "O-Glycan core 1", alignment = "substructure")
    Condition
      Warning:
      The provided alignment type "substructure" is different from the motif's alignment type "core" in database.
    Output
      [1] TRUE

