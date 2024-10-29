make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  glycan <- glyrepr::convert_glycan_mono_type(glycan, to = mono_type, strict = FALSE)
  if (!linkage) {
    glycan <- glyrepr::remove_linkages(glycan)
  }
  glycan
}
