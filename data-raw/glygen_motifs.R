library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(glyparse)
library(glyrepr)

source_metadata <- tribble(
  ~source_id , ~source                 ,
  "GGM"      , "GlyGen Motifs"         ,
  "CCRC"     , "CCRC Motifs"           ,
  "GTC"      , "GlyTouCan Motifs"      ,
  "GE"       , "GlycoEpitope Epitopes" ,
  "GD"       , "Glydin"                ,
  "GDB"      , "Glydin - BiOligo"      ,
  "GDV"      , "Glydin - Cermav"       ,
  "GDC"      , "Glydin - Cummings"     ,
  "GDH"      , "Glydin - Hayes"        ,
  "GDSB"     , "Glydin - SugarBind"    ,
  "UCM"      , "UniCarbKB Motifs"      ,
  "GM"       , "All Motifs"
)

read_motif_collection <- function(path) {
  source_id <- tools::file_path_sans_ext(basename(path))
  source <- source_metadata$source[match(source_id, source_metadata$source_id)]

  if (is.na(source)) {
    cli::cli_abort("Unknown motif source file: {.file {path}}.")
  }

  motif_data <- read_json(path, simplifyDataFrame = TRUE)[[
    "results"
  ]][["bindings"]] |>
    as_tibble() |>
    mutate(across(everything(), ~ .x$value))

  if (!"name" %in% names(motif_data)) {
    motif_data$name <- NA_character_
  }

  motif_data |>
    mutate(
      source_id = source_id,
      source = source,
      alignment = case_when(
        .data$alignment == "Substructure" ~ "substructure",
        .data$alignment == "Whole-Glycan" ~ "whole",
        .data$alignment == "Nonreducing-End" ~ "terminal",
        .data$alignment == "Core" ~ "core",
        .default = NA_character_
      )
    )
}

df <- list.files("data-raw", pattern = "[.]json$", full.names = TRUE) |>
  keep(~ basename(.x) != "glygen_motifs.json") |>
  map_dfr(read_motif_collection) |>
  mutate(glycan_structure = parse_wurcs(.data$wurcs, on_failure = "na"))

good_df <- df |>
  filter(!is.na(.data$glycan_structure))

good_df[good_df$source_id == "GGM" & good_df$accession == "001020", "name"] <-
  "O-Fucose Core 1"
good_df[good_df$source_id == "GGM" & good_df$accession == "001021", "name"] <-
  "O-Fucose Core 2"
good_df[
  good_df$source_id == "GGM" & good_df$accession == "001022",
  "name"
] <- "N-glycan core, GlcNAc truncated"
good_df[good_df$source_id == "GGM" & good_df$accession == "001041", "name"] <-
  "N-glycan core, Man truncated"
good_df[good_df$source_id == "GGM" & good_df$accession == "000106", "name"] <-
  "GD1alpha"

glygen_motifs <- good_df |>
  mutate(glycan_structure = fill_anomer_pos(.data$glycan_structure)) |>
  select(all_of(c(
    "source_id",
    "source",
    "accession",
    "name",
    "alignment",
    "glycan_structure"
  )))

collection_counts <- glygen_motifs |>
  count(.data$source_id, .data$source, name = "motifs")

print(collection_counts, n = Inf)

usethis::use_data(glygen_motifs, internal = TRUE, overwrite = TRUE)
