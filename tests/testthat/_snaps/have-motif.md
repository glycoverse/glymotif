# wrong glycan types

    Code
      have_motif(glycan, motif)
    Condition
      Error in `prepare_have_motif_args()`:
      ! `glycans` must be a 'glyrepr_structure' object or an IUPAC-condensed structure character.
      x The input is of class <igraph>.

# wrong motif types

    Code
      have_motif(glycan, motif)
    Condition
      Error in `prepare_have_motif_args()`:
      ! The `motif` argument must be a scalar vector.
      x The input is of length 0.

# unkown motif name used as input

    Code
      have_motif(glycan, motif)
    Condition
      Error in `prepare_have_motif_args()`:
      ! `motifs` must be a 'glyrepr_structure' object,a character vector of IUPAC-condensed structure strings,or a character vector of known motif names.
      x Some motifs are neither valid IUPAC-condensed structures nor known motif names.
      i Use `available_motifs()` to see all valid motif names.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: {.val {x}}
      i IUPAC-condensed string cannot contain whitespace

# bad glycan IUPAC

    Code
      have_motif(glycan, motif)
    Condition
      Error in `prepare_have_motif_args()`:
      ! `glycans` must be a 'glyrepr_structure' object or an IUPAC-condensed structure character.
      x Some glycans could not be parsed as valid IUPAC-condensed structures.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: {.val {x}}
      i IUPAC-condensed string cannot contain whitespace

# warning when user-provided alignment is different from database

    Code
      have_motif(glycan, motif, alignment = "terminal")
    Condition
      Warning:
      The provided alignment type "terminal" is different from the motif's alignment type "core" in database for motif "O-Glycan core 1".
    Output
      [1] FALSE

# custom alignment is different from database

    Code
      have_motif("Gal(b1-3)GalNAc(a1-3)GlcNAc", "O-Glycan core 1", alignment = "substructure")
    Condition
      Warning:
      The provided alignment type "substructure" is different from the motif's alignment type "core" in database for motif "O-Glycan core 1".
    Output
      [1] TRUE

