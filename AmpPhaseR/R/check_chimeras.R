
# return a tibble with columns: seq, count, freq, is_candidate, is_chimera, parents
#' @importFrom purrr map_df map2_chr keep discard map_lgl
#' @importFrom magrittr "%<>%" set_names
check_chimeras <- function(seq_tbl, max_phases, max_breakpoints) {

  stopifnot(is.data.frame(seq_tbl),
            all(c('seq', 'count', 'freq') %in% colnames(seq_tbl)),
            is.integer(seq_tbl$count),
            is.integer(max_phases) & max_phases >= 1,
            is.character(seq_tbl$seq),
            is.numeric(seq_tbl$freq),
            all(nchar(seq_tbl$seq) == nchar(seq_tbl$seq[1])),
            !any(is.na(seq_tbl)))

  # rules:
  # 1) sum of parents > sum of chimeras
  # 2) All phases >= min parent must be assigned

  if (nrow(seq_tbl) == 1 || max_phases == 1) {
    return(tibble(n_phases = integer(),
                  seq_parent = list(),
                  seq_child = list()))
  }

  res <-
    tibble(n_phases = seq.int(2, min(max_phases, nrow(seq_tbl)))) %>%
    mutate(data = map(n_phases, function(n) {
      combn_tbl(seq_tbl$seq, n, name = 'parent')
    })) %>%
    unnest(data) %>%
    mutate(non_parent = map(parent, ~setdiff(seq_tbl$seq, .))) %>%
    mutate(child_data = map(non_parent, function(s) {
      map_df(seq.int(0, length(s)), function(n){
        combn_tbl(s, n, name = 'child') %>%
          mutate(other = map(child, ~setdiff(s, .)))
      })
    })) %>%
    select(-non_parent) %>%
    unnest(child_data) %>%
    mutate(id = seq_len(n())) %>%
    pivot_longer(c(parent, child, other),
                 values_to = 'seq') %>%
    unchop(seq) %>%
    left_join(select(seq_tbl, seq, freq), 'seq') %>%
    chop(c(seq, freq)) %>%
    pivot_wider(names_from = name,
                values_from = c(seq, freq)) %>%
    (function(x) {
      if (!'seq_child' %in% colnames(x)) {
        x <- mutate(x, seq_child = list(character()), freq_child = list(double()))
      }
      if (!'seq_other' %in% colnames(x)) {
        x <- mutate(x, seq_other = list(character()), freq_other = list(double()))
      }
      return(x)
    }) %>%
    filter(map_dbl(freq_parent, sum) > map_dbl(freq_child, sum),
           map_dbl(freq_parent, min) > suppressWarnings(map_dbl(freq_other, max))) %>%
    filter(map2_lgl(seq_parent, seq_child, function(parents,children) {
      all(map_lgl(children, ~ are_parents(., parents, max_breakpoints)))
    })) %>%
    mutate(cov_p = map_dbl(freq_parent, sum),
           cov_c = map_dbl(freq_child, sum)) %>%
    arrange(n_phases, desc(cov_p + cov_c), cov_p) %>%
    group_by(n_phases) %>%
    slice(1) %>%
    ungroup() %>%
    select(n_phases, seq_parent, seq_child)

  return(res)
}

are_parents <- function(subject, pool, max_breakpoints) {

  stopifnot(is_scalar_character(subject),
            is_character(pool),
            all(nchar(subject) == nchar(pool)),
            is_scalar_integerish(max_breakpoints))

  n_par <- length(pool)

  amat <-
    c(subject, pool) %>%
    str_split('', simplify = T)

  is_unique <- map_lgl(seq_len(ncol(amat)), function(i) !any(amat[1, i] == amat[-1, i]))

  if (any(is_unique)) {
    # not a chimera
    return(FALSE)
  }

  # paths
  result <-
    tibble(pos = seq_len(ncol(amat)),
           par = map(pos, ~ which(amat[1, .] == amat[1+ seq_len(n_par), .]))) %>%
    filter(lengths(par) < n_par) %>%
    mutate(par_hash = map_chr(par, digest::digest),
           group_start = replace_na(par_hash != lag(par_hash), TRUE),
           group = cumsum(group_start)) %>%
    group_by(group) %>%
    summarise(pos = list(pos), par = par[1]) %>%
    with(set_names(par, str_c('G', group))) %>%
    do.call(expand_grid, .) %>%
    mutate(id = seq_len(n())) %>%
    pivot_longer(starts_with('G'),
                 names_prefix = 'G',
                 names_transform = list(group = as.integer),
                 names_to = 'group',
                 values_to = 'parent') %>%
    arrange(id, group) %>%
    group_by(id) %>%
    summarise(parents = list(parent %>% sort %>% unique), breaks = sum(parent != lead(parent), na.rm = T))

  return(any(result$breaks) < max_breakpoints)
}


# return list of possible parents
#' @importFrom rlang is_scalar_character is_character
#' @importFrom dplyr lag lead
chimera_parents <- function(subject, pool, max_breakpoints) {

  stopifnot(is_scalar_character(subject),
            is_character(pool),
            all(nchar(subject) == nchar(pool)),
            is_scalar_integerish(max_breakpoints))

  n_par <- length(pool)

  amat <-
    c(subject, pool) %>%
    str_split('', simplify = T)

  is_unique <- map_lgl(seq_len(ncol(amat)), function(i) !any(amat[1, i] == amat[-1, i]))

  if (any(is_unique)) {
    # not a chimera
    return(list())
  }

  # paths
  result <-
    tibble(pos = seq_len(ncol(amat)),
           par = map(pos, ~ which(amat[1, .] == amat[1+ seq_len(n_par), .]))) %>%
    filter(lengths(par) < n_par) %>%
    mutate(par_hash = map_chr(par, digest::digest),
           group_start = replace_na(par_hash != lag(par_hash), TRUE),
           group = cumsum(group_start)) %>%
    group_by(group) %>%
    summarise(pos = list(pos), par = par[1]) %>%
    with(set_names(par, str_c('G', group))) %>%
    do.call(expand_grid, .) %>%
    mutate(id = seq_len(n())) %>%
    pivot_longer(starts_with('G'),
                 names_prefix = 'G',
                 names_transform = list(group = as.integer),
                 names_to = 'group',
                 values_to = 'parent') %>%
    arrange(id, group) %>%
    group_by(id) %>%
    summarise(parents = list(parent %>% sort %>% unique), breaks = sum(parent != lead(parent), na.rm = T)) %>%
    filter(breaks <= max_breakpoints) %>%
    pull(parents)

  return(result)
}

est_parent_freq <- function(seq_tbl) {
  stopifnot(is.data.frame(seq_tbl),
            all(c('seq', 'count', 'is_candidate') %in% colnames(seq_tbl)),
            is.integer(seq_tbl$count),
            is.character(seq_tbl$seq),
            is.logical(seq_tbl$is_candidate),
            all(nchar(seq_tbl$seq) == nchar(seq_tbl$seq[1])),
            !any(is.na(seq_tbl)))

  a_mat <-
    seq_tbl$seq %>%
    (function(s) str_split(s, '', simplify = T) %>% set_rownames(s))

  parents <- seq_tbl %>% filter(is_candidate) %>% pull(seq)
  pi <- which(rownames(a_mat) %in% parents)

  p_combs <- map_lgl(seq_len(nrow(a_mat)), function(i) {
    all(map_lgl(seq_len(ncol(a_mat)), function(j){
      a_mat[i, j] %in% a_mat[pi, j]
    }))
  })
  a_mat <- a_mat[p_combs, , drop = FALSE]

  freq_est <-
    map_df(parents, function(seq) {
      uniq_pos <-
        map_lgl(seq_len(ncol(a_mat)),
                function(j) !any(a_mat[seq, j] == a_mat[setdiff(parents, seq), j])) %>%
        which()

      if (length(uniq_pos) > 0) {
        sweep(a_mat[, uniq_pos, drop=FALSE], 2, a_mat[seq, uniq_pos, drop=FALSE], '==' ) %>%
          set_colnames(str_c('P', seq_along(uniq_pos))) %>%
          as_tibble() %>%
          mutate(count = seq_tbl$count[p_combs]) %>%
          pivot_longer(starts_with('P'),
                       names_to = 'index',
                       names_prefix = 'P',
                       names_transform = list(index = as.integer),
                       values_to  = 'is_match') %>%
          group_by(index, is_match) %>%
          summarise(count = sum(count)) %>%
          summarise(freq = sum(count[-is_match]) / sum(count)) %>%
          summarise(freq = mean(freq, na.rm = TRUE)) %>%
          mutate(seq = seq)
      } else {
        tibble(seq = seq, freq = NA_real_)
      }
    }) %>%
    (function(x) {
      if (!any(is.na(x$freq))) {
        mutate(x, freq = freq / sum(freq))
      } else if (sum(is.na(x$freq) == 1) && sum(x$freq, na.rm = TRUE) < 1) {
        mutate(x, freq = replace(freq, which(is.na(freq)), 1 - sum(freq, na.rm = TRUE)))
      } else {
        seq_tbl %>%
          filter(is_candidate) %>%
          mutate(freq = count / sum(count)) %>%
          select(freq, seq)
      }
    })


  return(freq_est)

}
