# describing glycans works

    Code
      describe_n_glycans(glycans)
    Output
      # A tibble: 3 x 7
        glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal n_terminal_gal
        <chr>       <lgl>          <int>      <int>     <int> <int>          <int>
      1 highmannose FALSE             NA          0         0     0              0
      2 complex     FALSE              2          0         0     1              1
      3 complex     FALSE              2          1         1     0              0

# strictly describing glycans works

    Code
      describe_n_glycans(glycans, strict = TRUE)
    Output
      # A tibble: 3 x 7
        glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal n_terminal_gal
        <chr>       <lgl>          <int>      <int>     <int> <int>          <int>
      1 highmannose FALSE             NA          0         0     0              0
      2 complex     FALSE              2          0         0     1              1
      3 complex     FALSE              2          1         1     0              0

# describe N-glycans with no name

    Code
      result
    Output
      # A tibble: 3 x 7
        glycan_type bisecting n_antennae n_core_fuc n_arm_fuc n_gal n_terminal_gal
        <chr>       <lgl>          <int>      <int>     <int> <int>          <int>
      1 highmannose FALSE             NA          0         0     0              0
      2 complex     FALSE              2          0         0     1              1
      3 complex     FALSE              2          1         1     0              0

