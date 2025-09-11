.onLoad <- function(libname, pkgname) {
  rlang::run_on_load()

  # 设置memoise缓存以优化性能
  # 按照memoise包推荐的最佳实践，在包加载时设置缓存
  match_single_sub <<- memoise::memoise(match_single_sub)
  match_anomer <<- memoise::memoise(match_anomer)
  parse_linkage <<- memoise::memoise(parse_linkage)
  parse_anomer <<- memoise::memoise(parse_anomer)
  parse_pos2 <<- memoise::memoise(parse_pos2)
  get_motif_structure <<- memoise::memoise(get_motif_structure)
  get_motif_alignment <<- memoise::memoise(get_motif_alignment)
  get_motif_aglycon <<- memoise::memoise(get_motif_aglycon)
  terminal_nodes <<- memoise::memoise(terminal_nodes)
  core_node <<- memoise::memoise(core_node)
}
