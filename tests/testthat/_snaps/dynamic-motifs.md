# special branch options require logical scalars

    Code
      branch_motifs(include_empty_branch = NA)
    Condition
      Error in `branch_motifs()`:
      ! Assertion on 'include_empty_branch' failed: May not be NA.

---

    Code
      branch_motifs(include_mannose_branch = 1)
    Condition
      Error in `branch_motifs()`:
      ! Assertion on 'include_mannose_branch' failed: Must be of type 'logical flag', not 'double'.

---

    Code
      extract_branch_motif(glycan, include_empty_branch = NA)
    Condition
      Error in `extract_branch_motif()`:
      ! Assertion on 'include_empty_branch' failed: May not be NA.

---

    Code
      extract_branch_motif(glycan, include_mannose_branch = 1)
    Condition
      Error in `extract_branch_motif()`:
      ! Assertion on 'include_mannose_branch' failed: Must be of type 'logical flag', not 'double'.
