# describing glycans works

    Code
      describe_n_glycans(glycans)
    Output
      # A tibble: 3 x 8
        glycan             glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal
        <chr>              <chr>       <lgl>          <int>      <int>     <int> <int>
      1 Hex(??-?)[Hex(??-~ highmannose FALSE             NA          0         0     0
      2 Hex(??-?)HexNAc(?~ complex     FALSE              2          0         0     1
      3 dHex(??-?)HexNAc(~ complex     FALSE              2          1         1     0
      # i 1 more variable: n_terminal_gal <int>

# strictly describing glycans works

    Code
      describe_n_glycans(glycans, strict = TRUE)
    Output
      # A tibble: 3 x 8
        glycan             glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal
        <chr>              <chr>       <lgl>          <int>      <int>     <int> <int>
      1 Man(a1-3)[Man(a1-~ highmannose FALSE             NA          0         0     0
      2 Gal(b1-4)GlcNAc(b~ complex     FALSE              2          0         0     1
      3 Fuc(a1-3)GlcNAc(b~ complex     FALSE              2          1         1     0
      # i 1 more variable: n_terminal_gal <int>

# describe not N-glycans

    Code
      describe_n_glycans(glycans)
    Condition
      Error in `describe_n_glycans()`:
      ! Glycans at indices 1 are not N-glycans.

# describe N-glycans with no name

    Code
      result
    Output
      # A tibble: 3 x 8
        glycan             glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal
        <chr>              <chr>       <lgl>          <int>      <int>     <int> <int>
      1 Hex(??-?)[Hex(??-~ highmannose FALSE             NA          0         0     0
      2 Hex(??-?)HexNAc(?~ complex     FALSE              2          0         0     1
      3 dHex(??-?)HexNAc(~ complex     FALSE              2          1         1     0
      # i 1 more variable: n_terminal_gal <int>

