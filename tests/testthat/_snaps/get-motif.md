# db_motifs() creates a db_motifs_spec object

    Code
      db_motifs()
    Message
      <<db_motifs_spec>>
      This object should be passed to the `motifs` argument of `have_motifs()`,
      `count_motifs()`, `match_motifs()`, `add_motifs_lgl()`, or `add_motifs_int()`.
      Configuration: uses all GlyGen GlycoMotif database motifs

# getting motif structure works

    Code
      result
    Output
      <glycan_structure[1]>
      [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-
      # Unique structures: 1

