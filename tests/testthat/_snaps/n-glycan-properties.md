# describing glycans works

    Code
      describe_n_glycans(glycans)
    Output
      # A tibble: 3 x 7
        glycan glycan_type bisecting antennae core_fuc arm_fuc terminal_gal
        <chr>  <chr>       <lgl>        <int>    <int>   <int>        <int>
      1 H5N2   highmannose FALSE           NA        0       0            0
      2 H4N4   complex     FALSE            2        0       0            1
      3 H3N4F2 complex     FALSE            2        1       1            0

# strictly describing glycans works

    Code
      describe_n_glycans(glycans, strict = TRUE)
    Output
      # A tibble: 3 x 7
        glycan glycan_type bisecting antennae core_fuc arm_fuc terminal_gal
        <chr>  <chr>       <lgl>        <int>    <int>   <int>        <int>
      1 H5N2   highmannose FALSE           NA        0       0            0
      2 H4N4   complex     FALSE            2        0       0            1
      3 H3N4F2 complex     FALSE            2        1       1            0

# describe not N-glycans

    Code
      describe_n_glycans(glycans)
    Condition
      Error in `describe_n_glycans()`:
      ! Glycans at indices 1 are not N-glycans.

