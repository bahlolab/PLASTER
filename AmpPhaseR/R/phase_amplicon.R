
#' @export
#' @importFrom rlang is_scalar_integerish is_scalar_double is_scalar_character
#' @importFrom dplyr filter pull desc n bind_rows slice bind_cols full_join first last
#' @importFrom purrr map_chr map_dbl
#' @importFrom ggplot2 ggplot aes geom_tile geom_rect element_blank margin scale_fill_manual theme ggtitle guides guide_legend coord_flip geom_bar geom_col
phase_amplicon <- function(gds,
                           min_copy_num = 1,
                           max_copy_num = 2,
                           max_freq_abs_delta = 0.50,
                           max_freq_rel_delta = 0.50,
                           max_breakpoints = 3,
                           max_allele_freq = 0.90,
                           min_phase_freq = 0.05,
                           name = NULL) {

  stopifnot(is_scalar_integerish(min_copy_num) && min_copy_num >=1,
            is_scalar_integerish(max_copy_num) && max_copy_num >= min_copy_num,
            is_scalar_integerish(max_breakpoints) && max_breakpoints >= 1,
            is_scalar_double(max_freq_abs_delta) && max_freq_abs_delta > 0,
            is_scalar_double(max_freq_rel_delta) && max_freq_rel_delta > 0,
            is_scalar_double(max_allele_freq) && max_allele_freq > 0,
            is_scalar_double(min_phase_freq) && min_phase_freq > 0,
            is.null(name) || is_scalar_character(name))

  alleles <- allele_matrix(gds)

  if (any(is.na(alleles))) {
    stop('NA variant calls not yet implemented.')
    # NAs are unacceptable
    # Need strategy to remove
    # suggest calculate delta if variant removed or sample removed
    # remove highest point recursively
  }

  var_smry <- summarise_vars(alleles)
  vars_keep <-
    var_smry %>%
    group_by(variant) %>%
    filter(all(freq < max_allele_freq)) %>%
    ungroup() %>%
    pull(variant) %>%
    unique() %>%
    str_remove('V') %>%
    as.integer() %>%
    sort()

  if (length(vars_keep) > 0) {
    alleles_sub <- alleles[, vars_keep, drop = FALSE]
  } else {
    alleles_sub <- matrix(rep('R', nrow(alleles)), ncol = 1)
  }

  derep_tbl <-
    tibble(read_index = seq_len(nrow(alleles_sub)),
           seq = map_chr(read_index, ~str_c(alleles_sub[., ], collapse = ''))) %>%
    group_by(seq) %>%
    summarise(count = n(), indices = list(read_index)) %>%
    arrange(desc(count)) %>%
    mutate(is_singleton = count == 1,
           freq = if_else(!is_singleton,
                          count / sum(count[!is_singleton]),
                          NA_real_))

  if (sum(derep_tbl$freq > min_phase_freq, na.rm = T)==0) {
    return(list(success = FALSE,
                note = 'all Singletons or too few reads'))
  }

  cand_phases <-
    derep_tbl %>%
    select(-indices) %>%
    filter(freq > min_phase_freq) %>%
    check_chimeras(max_breakpoints = max_breakpoints,
                    max_phases = max_copy_num)

  phased <-
    `if`(nrow(cand_phases),
         cand_phases %>%
           select(n_phase = n_phases, seq = seq_parent) %>%
           pmap_df(function(n_phase, seq) {
             derep_tbl %>%
               filter(seq %in% !!seq) %>%
               select(seq, count) %>%
               (function(x) {
                 # estimate the original frequencies of each phase without chimeras
                 freq_est <-
                   filter(derep_tbl, !is_singleton) %>%
                   select(seq, count) %>%
                   anti_join(x, 'seq') %>%
                   mutate(is_candidate = FALSE) %>%
                   bind_rows(mutate(x, is_candidate = TRUE)) %>%
                   est_parent_freq()

                 ratio <-
                   enum_phase_ratios(n_phase,
                                     min_copy_num = max(n_phase, min_copy_num),
                                     max_copy_num = max_copy_num) %>%
                   mutate(abs_freq_delta = map_dbl(ratio,
                                                   function(r) sum(abs( freq_est$freq - (r / sum(r))))),
                          abs_rel_delta = map_dbl(ratio,
                                                  function(r) max(abs(freq_est$freq - (r / sum(r))) / (r / sum(r))))) %>%
                   slice(which.min(abs_freq_delta))

                 left_join(x, freq_est, by = 'seq') %>%
                   mutate(n_phase = n_phase) %>%
                   chop(-n_phase) %>%
                   bind_cols(ratio)
               })
           }) %>%
           filter(abs_freq_delta < max_freq_abs_delta,
                  abs_rel_delta < max_freq_rel_delta) %>%
           slice(which.max(n_phase)),
         NULL)

  if (is.null(phased) || nrow(phased) == 0) {
    phased <-
      slice(derep_tbl, 1) %>%
      select(seq, count) %>%
      mutate(n_phase = 1,
             freq = 1,
             ratio = 1,
             copy_num = min_copy_num,
             abs_freq_delta = 0,
             abs_rel_delta = 0) %>%
      chop(c(seq, count, freq, ratio, copy_num))
  }

  if (phased$n_phase >= 2) {
    chimeras <-
      derep_tbl %>%
      filter(!seq  %in% phased$seq[[1]],
             !is_singleton) %>%
      filter(map_lgl(seq, function(x) {
        are_parents(x, phased$seq[[1]], max_breakpoints = max_breakpoints)
      })) %>%
      select(seq)
  } else {
    chimeras <- tibble(seq = character())
  }


  read_phase <-
    phased %>%
    select(seq) %>%
    unnest(seq) %>%
    mutate(phase = str_c('phase_', seq_along(seq))) %>%
    # TODO: actually ensure that parents are only phase_1 and phase_2 ?
    bind_rows(chimeras %>%  mutate(phase = 'chimera')) %>%
    full_join(select(derep_tbl, seq, indices, is_singleton), 'seq') %>%
    mutate(is_singleton = replace_na(is_singleton, FALSE),
           phase = replace(phase, is_singleton, 'singleton'),
           phase = replace_na(phase, 'other')) %>%
    select(-seq) %>%
    unnest(indices) %>%
    mutate(read_name = rownames(alleles)[indices]) %>%
    select(read_name, phase) %>%
    arrange(read_name)

  ### plot results
  read_order <-
    as.data.frame(alleles) %>%
    dplyr::mutate_all(as.factor) %>%
    cluster::daisy(metric = 'gower') %>%
    hclust() %>%
    {rownames(alleles)[.$order]}

  phase_summary <-
    phased %>%
    mutate(copy_num = map_chr(copy_num, ~str_c(., collapse = ';'))) %>%
    unnest(c(seq, count, freq, ratio)) %>%
    mutate(phase = str_c('phase_', seq_along(seq))) %>%
    select(phase, count, freq, ratio, copy_num, abs_freq_delta, abs_rel_delta)


  phase_plot <-
    as.data.frame(alleles) %>%
    rownames_to_column('read_name') %>%
    left_join(read_phase, 'read_name') %>%
    mutate(phase = factor(phase,
                          levels = c(setdiff(phase, c('other', 'chimera', 'singleton')) %>% sort(),
                                     c('chimera', 'other', 'singleton')))) %>%
    mutate(read_id = match(read_name, read_order)) %>%
    arrange(desc(phase), read_id) %>%
    mutate(read_id = seq_len(n())) %>%
    select(-read_name) %>%
    pivot_longer(starts_with('V'),
                 names_to = 'variant',
                 names_prefix = 'V',
                 names_transform = list(variant = as.integer),
                 values_to = 'allele',
                 values_transform = list(allele = as.character)) %>%
    mutate(allele = if_else(allele == 'R', 'REF', allele),
           allele = if_else(allele == 'D', 'DEL', allele),
           allele = factor(allele, c('REF', 'DEL', 'A', 'C', 'G', 'T'))) %>%
    (function(x) {
      bounds <-
        group_by(x, phase) %>%
        summarise(xmin = min(variant)-0.5, xmax = max(variant)+0.5,
                  ymin = min(read_id)-0.5, ymax = max(read_id) + 0.5)
      p1 <-
        ggplot(x, aes(variant, read_id, fill = allele)) +
        geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  data = bounds, col = NA, fill = 'gray70', inherit.aes = FALSE) +
        geom_tile() +
        scale_fill_manual(values = c(REF = 'gray70',
                                     DEL = 'gray25',
                                     `T` = 'brown3',
                                     A = 'chartreuse4',
                                     G = 'darkgoldenrod3',
                                     C = 'dodgerblue3')) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank(),
              plot.margin = margin(),
              legend.position = 'top') +
        geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  data = bounds, col = 'white', fill = NA, inherit.aes = FALSE)

      l1 <- cowplot::get_legend(p1)
      p1 <- p1 + guides(fill = FALSE)

      phase_pal <-
        c(setdiff(bounds$phase, c('other', 'chimera', 'singleton')) %>% sort(),
          c('chimera', 'other', 'singleton')) %>%
        set_names(scales::hue_pal()(length(.)), .)

      p2 <-
        bounds %>%
        ggplot(aes(fill = phase)) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank(),
              plot.margin = margin(),
              legend.position = 'top') +
        guides(fill = guide_legend(nrow = 2)) +
        scale_fill_manual(values = phase_pal) +
        geom_rect(aes(xmin=0.5, xmax=1.5, ymin=ymin, ymax=ymax),
                  data = bounds, col = 'white')

      l2 <- cowplot::get_legend(p2)
      p2 <- p2 + guides(fill = FALSE)

      pleg <- cowplot::plot_grid(l1, l2, nrow=1)
      p12 <- cowplot::plot_grid(p1, p2, rel_widths = c(9, 1), align = 'h', axis = 'tb')

      p3 <-
        read_phase %>%
        mutate(phase = factor(phase, rev(levels(x$phase)))) %>%
        ggplot(aes(phase, fill = phase)) +
        geom_bar() +
        coord_flip() +
        scale_fill_manual(values = phase_pal) +
        guides(fill = FALSE)

      p4a <-
        phase_summary %>%
        mutate(prop_raw = count / sum(count),
               prop_corrected = freq) %>%
        select(phase, starts_with('prop')) %>%
        pivot_longer(starts_with('prop'),
                     names_prefix = 'prop_',
                     names_to = 'prop') %>%
        ggplot(aes(prop, value, fill = phase)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = phase_pal) +
        xlab('proportion') +
        ylab('')

      p4b <-
        phase_summary %>%
        ggplot(aes(x = 'assigned', ratio, fill = phase)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = phase_pal) +
        theme(axis.title.y = element_blank()) +
        guides(fill = FALSE)

      p4ab <- cowplot::plot_grid(p4a, p4b, ncol = 1, rel_heights=c(5,3), align = 'v', axis = 'lr')

      p34 <- cowplot::plot_grid(p3, p4ab, ncol = 2, rel_widths = c(1,1))

      cowplot::plot_grid(pleg, p12, p34, ncol = 1, rel_heights = c(1,4,2))
    })


  return(list(success = TRUE,
              read_phase = read_phase,
              phase_summary = phase_summary,
              phase_plot = phase_plot))
}
