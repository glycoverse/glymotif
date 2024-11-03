# `has_motifs()` works for motif names with custom alignments

    Code
      has_motifs(glycan, motifs, alignments = "substructure")
    Condition
      Warning:
      Use user-provided alignments, not the ones in the database.
    Output
      O-Glycan core 1 O-Glycan core 2 
                 TRUE            TRUE 

# some bad IUPAC

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `has_motifs()`:
      ! Motifs at indices 2 and 3 are neither known motif names or able to be parsed as IUPAC-condensed structure strings.

# some bad motif names

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `get_motifs_type()`:
      ! Motifs "bad1" and "bad2" are not known motif names.

# wrong types

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `get_motifs_type()`:
      ! `motifs` must be either 'glycan_graph' objects, a character vector of IUPAC-condensed structure strings, or a character vector of known motif names.

# missing motifs uses default motifs

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

# `have_motif()` with motif name and custom alignment

    Code
      have_motif(glycans, motif, alignment = "whole")
    Condition
      Warning:
      The provided alignment type "whole" is different from the motif's alignment type "core" in database.
    Output
         G1    G2 
      FALSE FALSE 

# some bad glycans

    Code
      have_motif(glycans, motif)
    Condition
      Error in `have_motif()`:
      ! Glycans at indices 2 and 3 are not able to be parsed as IUPAC-condensed structure strings.

