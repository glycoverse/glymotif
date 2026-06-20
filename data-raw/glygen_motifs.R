library(jsonlite)
library(tidyverse)
library(glyparse)
library(glyrepr)

df <- read_json("data-raw/glygen_motifs.json", simplifyDataFrame = TRUE)[[
  "results"
]][["bindings"]] |>
  as_tibble() |>
  mutate(across(everything(), ~ .x$value)) |>
  mutate(
    alignment = case_match(
      alignment,
      "Substructure" ~ "substructure",
      "Whole-Glycan" ~ "whole",
      "Nonreducing-End" ~ "terminal",
      "Core" ~ "core"
    )
  ) |>
  mutate(glycan_structure = parse_wurcs(wurcs, on_failure = "na"))

good_df <- df |>
  filter(!is.na(glycan_structure))

good_df[good_df$accession == "001020", "name"] <- "O-Fucose Core 1"
good_df[good_df$accession == "001021", "name"] <- "O-Fucose Core 2"
good_df[
  good_df$accession == "001022",
  "name"
] <- "N-glycan core, GlcNAc truncated"
good_df[good_df$accession == "001041", "name"] <- "N-glycan core, Man truncated"
good_df[good_df$accession == "000106", "name"] <- "GD1alpha"

glygen_motifs <- good_df |>
  mutate(glycan_structure = fill_anomer_pos(glycan_structure)) |>
  select(accession, name, alignment, glycan_structure)

usethis::use_data(glygen_motifs, internal = TRUE, overwrite = TRUE)
