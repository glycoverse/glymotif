make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  glycan <- glyrepr::convert_mono_type(glycan, to = mono_type)
  if (!linkage) {
    glycan <- glyrepr::remove_linkages(glycan)
  }
  glycan
}
