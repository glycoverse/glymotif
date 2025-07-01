# add_comp_descriptions adds correct columns

    Code
      exp_with_desc <- add_comp_descriptions(exp)

# add_comp_descriptions handles A and G correctly

    Code
      exp_with_desc <- add_comp_descriptions(exp)

# add_struct_descriptions adds description columns

    Code
      exp_with_desc <- add_struct_descriptions(exp)

# add_glycan_descriptions adds both composition and structure descriptions

    Code
      exp_with_desc <- add_glycan_descriptions(exp)
    Message
      v Structure descriptions added.
      v Composition descriptions added.

# add_comp_descriptions handles repeated calls correctly

    Code
      exp <- add_comp_descriptions(exp)

# add_struct_descriptions handles repeated calls correctly

    Code
      exp <- add_struct_descriptions(exp)

# add_glycan_descriptions handles repeated calls correctly

    Code
      exp <- add_glycan_descriptions(exp)
    Message
      v Structure descriptions added.
      v Composition descriptions added.

---

    Code
      exp2 <- add_glycan_descriptions(exp)
    Message
      i Structure descriptions already added. Skipping.
      i Composition descriptions already added. Skipping.

