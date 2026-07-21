make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  if (mono_type == "generic") {
    glycan <- glyrepr::convert_to_generic(glycan)
  }
  if (!linkage) {
    glycan <- glyrepr::remove_linkages(glycan)
  }
  glycan
}
