
is_gds <- function(x) {
  inherits(x, "SeqVarGDSClass")
}

cosine <- function(x, y) {
  stopifnot(is.vector(x),
            is.vector(y),
            length(x) == length(y))
  return(c(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))))
}



combn_tbl <- function(x, m, name = 'combs') {
  tibble(!!name :=  combn(x, m) %>%
          apply(2, list) %>%
          purrr::flatten())
}

#' @importFrom dplyr mutate group_by summarise select lead
#' @importFrom tidyr expand_grid replace_na
perm_tbl <- function(x, m, name = 'perms') {

  map(seq_len(m), ~ x) %>%
    setNames(str_c('i', seq_len(m))) %>%
    do.call(expand_grid, .) %>%
    mutate(row = seq_len(n())) %>%
    pivot_longer(starts_with('i')) %>%
    group_by(row) %>%
    filter(!any(replace_na(value == lead(value), FALSE))) %>%
    summarise(!!name := list(value)) %>%
    select(-row)
}

#' @importFrom dplyr n_distinct distinct transmute
#' @importFrom tidyr chop
#' @importFrom purrr reduce
enum_phase_ratios <- function(n_phase, min_copy_num = n_phase, max_copy_num = n_phase) {

  map_df(seq.int(min_copy_num, max_copy_num), function(ploidy) {
    rep(list(seq_len(n_phase)), ploidy) %>%
      setNames(seq_len(ploidy)) %>%
      do.call(expand_grid, .) %>%
      mutate(id = seq_len(n())) %>%
      pivot_longer(-id, names_to = 'copy', values_to = 'hap') %>%
      group_by(id, hap) %>%
      summarise(hap_count = n()) %>%
      group_by(id) %>%
      filter(n_distinct(hap) == n_phase,
             sum(hap_count) == ploidy) %>%
      mutate(ratio = `if`(n() == 1, 1L, as.integer(hap_count / reduce(hap_count, pracma::gcd)))) %>%
      ungroup() %>%
      select(id, hap, ratio) %>%
      pivot_wider(names_from = hap,
                  names_prefix = 'p',
                  values_from = ratio) %>%
      select(-id) %>%
      distinct() %>%
      mutate(copy_num = ploidy,
             id = seq_len(n())) %>%
      pivot_longer(starts_with('p'),
                   names_to = 'phase',
                   values_to = 'ratio') %>%
      group_by(copy_num, id) %>%
      summarise(ratio = list(ratio)) %>%
      ungroup() %>%
      select(-id)
    }) %>%
    mutate(ratio_hash = map_chr(ratio, digest::digest)) %>%
    group_by(ratio_hash) %>%
    summarise(ratio = ratio[1],
              copy_num = list(copy_num)) %>%
    ungroup() %>%
    select(-ratio_hash)
}


