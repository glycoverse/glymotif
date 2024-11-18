# describing glycans works

    Code
      describe_n_glycans(glycans)
    Output
      # A tibble: 3 x 8
        glycan glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal
        <chr>  <chr>       <lgl>          <int>      <int>     <int> <int>
      1 H5N2   highmannose FALSE             NA          0         0     0
      2 H4N4   complex     FALSE              2          0         0     1
      3 H3N4F2 complex     FALSE              2          1         1     0
      # i 1 more variable: n_terminal_gal <int>

# strictly describing glycans works

    Code
      describe_n_glycans(glycans, strict = TRUE)
    Output
      # A tibble: 3 x 8
        glycan glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal
        <chr>  <chr>       <lgl>          <int>      <int>     <int> <int>
      1 H5N2   highmannose FALSE             NA          0         0     0
      2 H4N4   complex     FALSE              2          0         0     1
      3 H3N4F2 complex     FALSE              2          1         1     0
      # i 1 more variable: n_terminal_gal <int>

# describe not N-glycans

    Code
      describe_n_glycans(glycans)
    Condition
      Error in `describe_n_glycans()`:
      ! Glycans at indices 1 are not N-glycans.

