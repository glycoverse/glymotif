# extract_branch_motif warns for glycans other than N-glycans

    ! Some of `glycans` do not have the N-glycan core motif.
    ! Those glycans are skipped.

# floating branch attributes retain the reference extraction path

    Code
      extract_branch_motif(glycan)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `floating_graph_info()`:
      ! A floating glycan structure must contain exactly one main tree.

