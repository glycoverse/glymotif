# some bad motifs

    Code
      has_motifs(glycan, motifs)
    Condition
      Error in `has_motifs()`:
      ! Motifs at indices 2 and 3 are neither known motif names or able to be parsed as IUPAC-condensed structure strings.

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

# some bad glycans

    Code
      have_motif(glycans, motif)
    Condition
      Error in `have_motif()`:
      ! Glycans at indices 2 and 3 are not able to be parsed as IUPAC-condensed structure strings.

