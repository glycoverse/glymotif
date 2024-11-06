# Paucimannose 1
ng_motifs <- c(
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     └─Man (a1-3)
  pman1 = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  pman2 = "Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-3)
  #     └─Man (a1-6)
  #       ├─Man (a1-6)
  #       └─Man (a1-3)
  hman = "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  hybrid = "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",

  # same as pman1
  core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
)
ng_motifs <- purrr::map(ng_motifs, glyparse::parse_iupac_condensed)
ng_s_motifs <- purrr::map(ng_motifs, glyrepr::convert_glycan_mono_type, to = "simple")


#' @importFrom magrittr %>%
n_glycan_type <- function(glycan) {
  valid_glycan_arg(glycan)
  glycan <- glycan %>%
    ensure_glycan_is_graph() %>%
    glyrepr::convert_glycan_mono_type(to = "simple", strict = FALSE)
  if (has_motif_(glycan, ng_s_motifs$pman1, alignment = "whole", ignore_linkages = TRUE) ||
      has_motif_(glycan, ng_s_motifs$pman2, alignment = "whole", ignore_linkages = TRUE)) {
    "paucimannose"
  } else if (has_motif_(glycan, ng_s_motifs$hybrid, alignment = "core", ignore_linkages = TRUE)) {
    "hybrid"
  } else if (has_motif_(glycan, ng_s_motifs$hman, alignment = "core", ignore_linkages = TRUE)) {
    "highmannose"
  } else if (has_motif_(glycan, ng_s_motifs$core, alignment = "core", ignore_linkages = TRUE)){
    "complex"
  } else {
    rlang::abort("Not an N-glycan")
  }
}
