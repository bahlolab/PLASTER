#!/usr/bin/env Rscript

stopifnot(require(readr),
          require(stringr),
          require(dplyr),
          require(magrittr),
          require(purrr),
          require(tidyr),
          require(stringr),
          require(docopt))

"
Usage:
  merge_lima_smry.R <ccs_smry> <sr_smry> <output> [options]
" -> doc

opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

stopifnot(file.exists(opts$ccs_smry), 
          file.exists(opts$sr_smry),
          !is.null(opts$output))

tibble(read_type = c('ccs', 'subread'),
       filename = c(opts$ccs_smry, opts$sr_smry)) %>% 
  mutate(data = map(filename, function(fn) {
    read_lines(fn) %>%
      str_split_fixed(':',2) %>%
      set_colnames(c('tag', 'value')) %>%
      as_tibble() %>%
      filter(value != '', tag != '') %>%
      mutate(tag = str_to_lower(tag) %>% str_remove_all('\\(.\\)') %>% str_replace_all('\\s+', '_') %>% str_remove('_$'),
             value = str_remove_all(value,' |%|\\([0-9.]+%?\\)')) %>%
      spread(tag, value)
  })) %>% 
  select(-filename) %>% 
  unnest(data) %>% 
  write_tsv(opts$output)

