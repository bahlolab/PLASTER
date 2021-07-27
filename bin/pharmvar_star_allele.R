#!/usr/bin/env Rscript

stopifnot(require(httr),
          require(SeqArray),
          require(jsonlite),
          require(tidyverse),
          require(assertthat),
          require(magrittr),
          require(docopt))

doc <- "
Usage:
  pharmvar_star_allele.R <sample_vcf> <pharmvar_vcf> [options]

Options:
  sample_vcf                  vcf(.gz) file with sample variants.
  pharmvar_vcf                vcf(.gz) file with PharmVar star allele variants.
  --amplicon=<f>              Json file with data about amplicon.
  --out-pref=<f>              Output file prefix.
"
opts <- docopt(doc)

## parse options, check inputs
amplicon <- fromJSON(opts$amplicon)
stopifnot(file.exists(opts$pharmvar_vcf),
          file.exists(opts$sample_vcf),
          is_scalar_character(amplicon$pharmvar_gene),
          is_scalar_character(opts$out_pref),
          is_scalar_integer(amplicon$start),
          is_scalar_integer(amplicon$end))

## function definitions
seq_open_vcf <- function(vcf_fn) {
  stopifnot(is_scalar_character(vcf_fn), 
            file.exists(vcf_fn))
  
  gds_fn <- str_replace(vcf_fn, '.vcf.gz', '.gds')
  if (! file.exists(gds_fn)) {
    seqVCF2GDS(vcf.fn = vcf_fn, out.fn = gds_fn, storage.option = 'ZIP_RA', ignore.chr.prefix = '')
  }
  seqOpen(gds_fn, allow.duplicate = T)
}

get_var_info <- function(gds) {
  tibble(
    variant_id = seqGetData(gds, 'variant.id'),
    chrom = seqGetData(gds, 'chromosome'),
    pos = seqGetData(gds, 'position'),
    allele = seqGetData(gds, 'allele')) %>% 
    mutate(allele_data = map(allele, ~ {
      alleles <- c(str_split(., ',', simplify = TRUE))
      tibble(ref = alleles[1], alt = alleles[-1], allele_index = seq_along(alt))
    })) %>% 
    select(-allele) %>% 
    unnest(allele_data) %>% 
    mutate(ref = str_to_upper(ref),
           alt = str_to_upper(alt))
}

get_vep_info <- function(gds, amplicon, ann_tag = 'CSQ') {
  
  vep_ann_names <- 
    header(gds)$INFO %>%
    as.data.frame() %>% rownames_to_column() %>% as_tibble() %>%
    filter(str_detect(rowname, ann_tag),
           str_detect(Description, 'Ensembl VEP'))%>%
    pull(Description) %>%
    str_remove('.+Format: ') %>% 
    str_split('\\|', simplify = T) %>% 
    c()
  
  vep_ann <- seqGetData(gds, str_c('annotation/info/', ann_tag))
  
  tibble(
    variant_id = seqGetData(gds, 'variant.id'),
    chromosome = seqGetData(gds, 'chromosome'),
    position = seqGetData(gds, 'position')) %>% {
      .[rep(1:nrow(.), times = vep_ann$length),] } %>%
    bind_cols(
      str_split_fixed(vep_ann$data, '\\|', length(vep_ann_names)) %>%
        set_colnames(vep_ann_names) %>% 
        as_tibble()) %>% 
    set_colnames(., str_to_lower(names(.)))  %>% 
    type_convert(., col_types = cols()) %>% 
    select_if(., ~ !all(is.na(.))) %>% 
    mutate(impact = ordered(impact, c('MODIFIER', 'LOW', 'MODERATE',  'HIGH'))) %>% 
    (function(x) {
      names(amplicon) %>% 
        keep(~ str_starts(., 'vep_')) %>% 
        str_remove('vep_') %>% 
        intersect(colnames(x)) %>% 
        (function(fields) {
          if (length(fields)) {
            fields %>% 
              map(function(field) {
                x[[field]] == amplicon[str_c('vep_', field)]
              }) %>% 
              accumulate('&') %>%
              unlist() %>% 
              { filter(x, .) } 
          } else {
            x
          } }) }) %>% 
    arrange(variant_id, desc(impact)) %>% 
    group_by(variant_id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(annot = str_c(str_extract(amino_acids, '^.+(?=/)'),
                         protein_position,
                         str_extract(amino_acids, '(?<=/).+$')),
           annot = case_when(str_detect(consequence, 'frameshift') ~ str_replace(annot, 'X$', 'fs'),
                             str_detect(consequence, 'splice') ~ 'splicing',
                             TRUE ~ annot))
}

get_pharmvar_func <- function(gene, alleles) {
  map_chr(alleles, function(allele) {
    response <- 
      str_c("https://www.pharmvar.org/api-service/alleles/",
            gene, '*', allele, "/function") %>% 
      GET()
    if (response$status_code != 200) {
      warning("failed to retrieve function for ", gene, '*', allele)
      return('error')
    }
    content(response, as = 'text')
  })
}

ji <- function(set1, set2) { length(intersect(set1, set2)) / length(union(set1, set2)) }

# load pharmvar allele info
pharmvar_gds <- seq_open_vcf(opts$pharmvar_vcf)

pv_var_info <- 
  SeqArray::info(pharmvar_gds) %>% 
  as_tibble() %>% 
  select(starts_with(str_c(amplicon$pharmvar_gene, '_'))) %>% 
  mutate(variant_id = seq_len(n())) %>% 
  pivot_longer(-variant_id, names_to = 'allele', names_prefix = str_c(amplicon$pharmvar_gene, '_')) %>% 
  filter(value) %>% 
  select(-value) %>% 
  chop(allele) %>% 
  inner_join(get_var_info(pharmvar_gds) %>% 
               select(-allele_index) %>% 
               mutate(vid = str_c(str_remove(chrom, 'chr'), '-', pos, '-', ref, '-', alt)),
             by = 'variant_id') %>% 
  mutate(in_range = pos >= amplicon$start & pos <= amplicon$end) %T>% 
  with({
    if (sum(!in_range)) {
      message('Warning: ', sum(!in_range), '/', length(in_range), ' variants outside amplicon region')
    }
  })

pv_alleles <-
  pv_var_info %>% 
  filter(in_range) %>% 
  select(allele, vid) %>% 
  unnest(allele) %>%
  chop(vid) %>% 
  mutate(core_id = str_extract(allele, '[^.]+'),
         sub_id = str_extract(allele, '(?<=\\.)[^.]+') %>% replace_na('0'),
         is_core = !str_detect(allele, '\\.')) %>% 
  add_row(allele = '1.001',
          core_id = '1',
          sub_id = '001',
          is_core = FALSE, 
          vid = list(character()) %>% vctrs::as_list_of()) %>% 
  arrange(as.numeric(core_id), as.numeric(sub_id)) %>% 
  (function(x) {
    ambig <-
      filter(x, !is_core) %>% 
      select(allele, vid) %>%
      chop(allele) %>% 
      filter(lengths(allele) > 1) %>% 
      mutate(rep = map_chr(allele, first)) %T>% 
      with(map_chr(allele, ~str_c(., collapse = ',')) %>% 
             str_c('(', ., ')') %>% 
             str_c(collapse = '; ') %>% 
             message('Warning: ambiguous alleles due to amplicon region are ', .)) %>% 
      mutate(ambig = map2(allele, rep, ~setdiff(.x, .y))) %>% 
      select(allele = rep, ambig)
    x %>% 
      filter(!allele %in% unlist(ambig$ambig)) %>% 
      left_join(ambig, 'allele')
  }) %>% 
  nest(data=-core_id) %>% 
  mutate(`function` = get_pharmvar_func(amplicon$pharmvar_gene, core_id) %>% 
           ordered(c('no function', 'decreased function', 'normal function'))) %>% 
  unnest(data)

# process sample allele data
sample_gds <- seq_open_vcf(opts$sample_vcf)

var_info <-
  get_var_info(sample_gds) %T>% 
  with(assert_that(all(allele_index == 1L))) %>% 
  select(-allele_index) %>% 
  mutate(vid = str_c(str_remove(chrom, 'chr'), '-', pos, '-', ref, '-', alt)) %>% 
  left_join(get_vep_info(sample_gds, amplicon) %>% 
              select(variant_id, symbol, gene, consequence, impact, feature, annot),
            by = "variant_id")

impacting_vars <-
  var_info %>% 
  filter(impact >= ordered('MODERATE', levels(impact))) %>% 
  pull(vid) %>% unique()

sample_gts <-
  seqGetData(sample_gds, 'genotype') %>% 
  {.[1, ,  ] } %>% 
  set_colnames(seqGetData(sample_gds, 'variant.id')) %>% 
  set_rownames(seqGetData(sample_gds, 'sample.id')) %>% 
  as.data.frame() %>% 
  rownames_to_column('sample_phase') %>% 
  separate(sample_phase, c('sample', 'phase'), '_ph') %>%
  mutate(phase = as.integer(phase)) %>% 
  pivot_longer(c(-sample, -phase),
               names_to = 'variant_id',
               names_transform = list(variant_id = as.integer),
               values_to = 'genotype')  %>% 
  left_join(select(var_info, variant_id, vid),
            by = "variant_id")

sample_phase_match <-
  sample_gts %>% 
  filter(genotype != 0) %>% 
  select(sample, phase, vid) %>% 
  chop(vid) %>% 
  complete(sample_gts %>% select(sample, phase) %>% distinct()) %>% 
  arrange_all() %>% 
  nest(sample_data = -vid) %>% 
  mutate(match = map(vid, function(vs) {
        # 1.001? i.e no variants
    if (length(vs) == 0) {
      return(
        tibble(sub_id = '001', core_id = '1',
               sim = 1, max_add_impact = ordered(NA, levels(var_info$impact)))
      )
    }
    # match core allele
    cid <-
      pv_alleles %>%
      filter(is_core) %>%
      filter(map_lgl(vid, function(x) all(x %in% vs))) %>%
      arrange(`function`,
              desc(lengths(vid)),
              as.numeric(core_id)) %>% 
      slice(1) %>%
      pull(core_id) %>%
      { `if`(length(.) == 0, '1', .)}

    # find closest matching sub allele
    sub_match <-
      pv_alleles %>%
      filter(!is_core,
             core_id == cid) %>% 
      mutate(sim = map_dbl(vid, ~ ji(., vs))) %>%
      arrange(desc(sim), as.numeric(sub_id))

    all_sub_vars <- sub_match$vid %>% unlist() %>% sort() %>% unique()

    max_add_impact <-
      var_info %>%
      filter(vid %in% setdiff(vs, all_sub_vars)) %>%
      pull(impact) %>%
      { suppressWarnings(max(.)) }

    return(
      sub_match %>%
        slice(1) %>%
        select(sub_id, core_id, sim) %>%
        mutate(max_add_impact = max_add_impact)
    )
  })) %>% 
  unnest(match) %>% 
  arrange(core_id, sub_id, sim) %>% 
  mutate(is_sub_match = sim == 1,
         is_func_equiv = { max_add_impact < ordered("MODERATE", levels(max_add_impact)) } %>% 
           replace_na(TRUE),
         vid_impact = map(vid, ~ intersect(., impacting_vars))) %>% 
  (function(x) {
    bind_rows(
      filter(x, is_func_equiv) %>% 
        mutate(core_id = str_c(core_id)) %>% 
        group_by(core_id, sub_id, is_sub_match)  %>% 
        mutate(sub_id = if_else(is_sub_match,
                                sub_id, 
                                str_c('novel-', str_pad(seq_along(sub_id), 2, 'l', '0')))) %>%
        ungroup() %>% 
        mutate(allele = str_c('*', core_id, '.', sub_id)),
      filter(x, !is_func_equiv) %>% 
        chop(-vid_impact) %>% 
        arrange(desc(map_int(sample_data, ~ sum(map_int(., nrow))))) %>% 
        mutate(core_id = str_c('N', seq_along(core_id))) %>% 
        unchop(-vid_impact) %>% 
        group_by(core_id) %>% 
        mutate(sub_id = str_c('novel-', str_pad(seq_along(sub_id), 2, 'l', '0'))) %>%
        ungroup() %>% 
        mutate(allele = str_c(core_id, '.', sub_id)))
  }) %>% 
  select(allele, core_id, sub_id, vid, vid_impact, sample_data) %>% 
  left_join(pv_alleles %>% select(core_id, `function`) %>% distinct(),
            by = 'core_id')

# write allele definitions
sample_phase_match %>% 
  mutate(variants_no_impact = 
           map2(vid, vid_impact, setdiff) %>% 
           map_chr(~ str_c(., collapse = ';')),
         vid_impact = map_chr(vid_impact, function(vi) {
           filter(var_info, vid %in% vi) %>%
             { `if`(nrow(.),
                    with(., str_c(vid, '(', annot, ')') %>% str_c(collapse = ';')),
                    '')
             }
         })) %>% 
  select(allele, `function`, variants_no_impact, variants_impact = vid_impact) %>% 
  arrange(allele) %>% 
  write_csv(str_c(opts$out_pref, '.allele_definition.csv'))

# write sample phase allele assignments
sample_phase_match %>% 
  mutate(core_allele = if_else(str_starts(core_id, 'N'), core_id, str_c('*', core_id))) %>% 
  unnest(sample_data) %>% 
  select(sample, phase, core_allele,  allele, `function`) %>% 
  arrange(sample, phase) %>% 
  write_csv(str_c(opts$out_pref, '.sample_phase_alleles.csv'))
