# has_motifs: motif names with custom alignments

    Code
      has_motifs(glycan, motifs, alignments = "substructure")
    Condition
      Warning:
      Use user-provided alignments, not the ones in the database.
    Output
      O-Glycan core 1 O-Glycan core 2 
                 TRUE            TRUE 

# has_motifs: custom alignments with wrong length

    Code
      has_motifs(glycan, motifs, alignments = alignments)
    Condition
      Error in `valid_alignments_arg()`:
      ! `alignments` must be either a single character string or a character vector of the same length as `motifs`.
      i `motif` length: 3, `alignments` length: 2

# has_motifs: some bad IUPAC

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `ensure_motifs_are_graphs()`:
      ! Motifs at indices 2 and 3 are neither known motif names or able to be parsed as IUPAC-condensed structure strings.

# has_motifs: some bad motif names

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `get_motifs_type()`:
      ! Motifs "bad1" and "bad2" are not known motif names.

# has_motifs: wrong types

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `valid_motifs_arg()`:
      ! `motifs` must be either a list of 'glycan_graph' objects, a character vector of IUPAC-condensed structure strings, or a character vector of known motif names.

# has_motifs: missing motifs uses default motifs

    Code
      result[1:10]
    Output
      Blood group H (type 2) - Lewis y                        i antigen 
                                 FALSE                            FALSE 
                              LacdiNAc                              GT2 
                                 FALSE                            FALSE 
      Blood group B (type 1) - Lewis b                            LcGg4 
                                 FALSE                            FALSE 
                Sialosyl paragloboside                   Sialyl Lewis x 
                                 FALSE                            FALSE 
                    A antigen (type 3)                       Type 1 LN2 
                                 FALSE                            FALSE 

# has_motifs: motifs have different mono types

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `has_motifs()`:
      ! All motifs must have the same monosaccharide type.

# has_motifs: glycan has obscurer mono type than motifs

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `ensure_glycan_mono_type()`:
      ! The monosaccharide type of `glycan` cannot be obscurer than `motif`.
      x "generic" is obscurer than "concrete".

# have_motif: motif names and custom alignment

    Code
      have_motif(glycans, motif, alignment = "whole")
    Condition
      Warning:
      The provided alignment type "whole" is different from the motif's alignment type "core" in database.
    Output
         G1    G2 
      FALSE FALSE 

# have_motif: some bad glycans

    Code
      have_motif(glycans, motif)
    Condition
      Error in `ensure_glycans_are_graphs()`:
      ! Glycans at indices 2 and 3 are not able to be parsed as IUPAC-condensed structure strings.

# have_motif: glycans have different mono types

    Code
      have_motif(glycans, motif)
    Condition
      Error in `have_motif()`:
      ! All glycans must have the same monosaccharide type.

# have_motif: glycans have obscurer mono types than motif

    Code
      have_motif(glycans, motif)
    Condition
      Error in `ensure_glycans_mono_type()`:
      ! The monosaccharide type of `glycans` cannot be obscurer than `motif`.
      x "generic" is obscurer than "concrete".

# have_motifs: alignments provided for known motif names

    Code
      have_motifs(glycans, motifs, alignments = "whole")
    Condition
      Warning:
      Use user-provided alignments, not the ones in the database.
    Output
                                                    O-Glycan core 1 O-Glycan core 2
      Gal(b1-3)GalNAc(a1-3)GlcNAc(b1-                         FALSE           FALSE
      GalNAc(a1-3)Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-           FALSE           FALSE

# have_motifs: alignments have wrong length

    Code
      have_motifs(glycans, motifs, alignments = alignments)
    Condition
      Error in `valid_alignments_arg()`:
      ! `alignments` must be either a single character string or a character vector of the same length as `motifs`.
      i `motif` length: 3, `alignments` length: 2

# have_motifs: bad motif names

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `get_motifs_type()`:
      ! Motifs "bad1" and "bad2" are not known motif names.

# have_motifs: bad IUPAC in motifs

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `ensure_motifs_are_graphs()`:
      ! Motifs at indices 2 and 3 are neither known motif names or able to be parsed as IUPAC-condensed structure strings.

# have_motifs: bad IUPAC in glycans

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `ensure_glycans_are_graphs()`:
      ! Glycans at indices 2 and 3 are not able to be parsed as IUPAC-condensed structure strings.

# have_motifs: glycans have different mono types

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `have_motifs()`:
      ! All glycans must have the same monosaccharide type.

# have_motifs: motifs have different mono types

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `have_motifs()`:
      ! All motifs must have the same monosaccharide type.

# have_motifs: glycans have obscurer mono types than motifs

    Code
      have_motifs(glycans, motifs)
    Condition
      Error in `ensure_glycans_mono_type()`:
      ! The monosaccharide type of `glycans` cannot be obscurer than `motif`.
      x "generic" is obscurer than "concrete".

