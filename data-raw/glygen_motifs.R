library(jsonlite)
library(tidyverse)

GFC_PATH <- "/Users/fubin/GlycanFormatConverter-cli/target/GlycanFormatConverter-cli.jar"

wurcs_to_iupac <- function(wurcs, gfc_path) {
  command <- paste(
    "java -jar",
    gfc_path,
    "-i WURCS -e IUPAC-Condensed",
    "-seq",
    paste0("'", wurcs, "'")
  )
  system(command, intern = TRUE)
}


df <- read_json("data-raw/glygen_motifs.json", simplifyDataFrame = TRUE)[["results"]][["bindings"]] %>%
  as_tibble() %>%
  mutate(across(everything(), ~ .x$value)) %>%
  mutate(iupac = map_chr(wurcs, wurcs_to_iupac, gfc_path = GFC_PATH)) %>%
  mutate(alignment = case_match(
    alignment,
    "Substructure" ~ "substructure",
    "Whole-Glycan" ~ "whole",
    "Nonreducing-End" ~ "terminal",
    "Core" ~ "core"
  ))

bad_iupac <- map_lgl(
  df$iupac,
  ~ tryCatch({
      glyparse::parse_iupac_condensed(.x)
      FALSE
    }
    , error = function(e) TRUE
  )
)
which(bad_iupac)
df[bad_iupac, c("accession", "iupac")]

good_df <- df
good_df[bad_iupac, ]$iupac <- paste0(good_df[bad_iupac, ]$iupac, "(?1-")
good_df[good_df$accession == "001024", "iupac"] <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-?)]GlcNAc(?1-"
good_df[good_df$accession == "001020", "name"] <- "O-Fucose Core 1"
good_df[good_df$accession == "001021", "name"] <- "O-Fucose Core 2"
good_df[good_df$accession == "001022", "name"] <- "N-glycan core, GlcNAc truncated"
good_df[good_df$accession == "001041", "name"] <- "N-glycan core, Man truncated"
good_df[good_df$accession == "000106", "name"] <- "GD1alpha"
good_df <- good_df[!good_df$accession %in% c("001026", "001027"), ]

glygen_motifs <- select(good_df, -wurcs)

usethis::use_data(glygen_motifs, internal = TRUE, overwrite = TRUE)
