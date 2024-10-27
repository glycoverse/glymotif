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

