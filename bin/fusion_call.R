#!/usr/bin/env Rscript

stopifnot(require(Biostrings),
          require(DECIPHER),
          require(Rsamtools),
          require(tidyverse),
          require(magrittr),
          require(digest),
          require(docopt))

doc <- "
Usage:
  fusion_call.R [options]

Options:
  --primary=<f>                 BAM file aligned to primary amplicon.
  --secondary=<f>               BAM file aligned to secondary amplicon.
  --amplicons=<f>               Json file with data about amplicons.
  --seq=<f>                     Fasta file with amplicon sequences.
  --out=<f>                     Prefix of output files [default: output].
  --max-reads=<f>               Maximum number of reads to process, randomly sampled to this number [default: Infinity].
  --seed=<f>                    Random seed for sampling down to maximum number of reads [default: 1].
  --min-score-delta=<f>         Minimum improvement in alignment score to call a fusion [default: 10].
  --min-fusion-reads=<f>        Minimum number of clustered fusion reads for a call [default: 5].
  --min-total-prop=<f>          Minimum proportion of total reads assigned to fusion for call [default: 0.01].
  --min-chim-prop=<f>           Minimum proportion of chimeric reads assigned to fusion for call [default: 0.25].
  --fusion-window=<f>           Window around breakpoint region for clustering [default: 10]
"
# opts <- '--primary LB-test-run-1.SM-NA17246.AM-CYP2D6.bam --secondary LB-test-run-1.SM-NA17246.AM-CYP2D7.bam \
#     --amplicons \'{"CYP2D6":{"chrom":"chr22","start":42125398,"end":42131503,"strand":"-"},"CYP2D7":{"chrom":"chr22","start":42137550,"end":42145176,"strand":"-"}}\' \
#     --seq chr22.fa \
#     --max-reads 400 \
#     --min-fusion-reads 5 \
#     --min-total-prop 0.01 \
#     --min-chim-prop 0.25 \
#     --min-score-delta 10 \
#     --fusion-window 10 \
#     --out SM-NA17246.fusion_check' %>% 
#   str_split('\\s+', simplify = T) %>% 
#   c() %>% 
#   docopt(doc, .)
opts <- docopt(doc)

# parse options
out_prefix <- opts$out
max_reads <- as.numeric(opts$max_reads)
seed <- as.integer(opts$seed)
min_score_delta <- as.integer(opts$min_score_delta)
min_fusion_reads <- as.integer(opts$min_fusion_reads)
min_total_prop <- as.numeric(opts$min_total_prop)
min_chim_prop <- as.numeric(opts$min_chim_prop)
fusion_window <- as.integer(opts$fusion_window)
amplicons <- jsonlite::parse_json(opts$amplicons)
bams <- c(opts$primary, opts$secondary)

stopifnot(length(bams) > 0,
          all(file.exists(bams)),
          file.exists(opts$seq))

seq <- readDNAStringSet(opts$seq)

ref_seq <-
  map2(amplicons, c('primary', 'secondary'), function(amp, nm) {
    subseq(seq, amp$start, amp$end) %>% 
      setNames(nm) %>% 
      { `if`(amp$strand == '-', reverseComplement(.), .) }
  }) %>% 
  unname() %>% 
  do.call(c, .)

bam_tbl <-
  tibble(ref = c('A', 'B'),
         amplicon = names(amplicons),
         is_rev = map_lgl(amplicon, ~ amplicons[[.]]$strand == '-'),
         file = list(opts$primary, opts$secondary)) %>% 
  unnest(file) %>% 
  mutate(reads = map2(file, is_rev, function(fn, rev) {
    scanBam(fn, param = ScanBamParam(what = c('seq', 'qname'))) %>% 
      first() %>% 
      with(setNames(seq, qname)) %>% 
      { `if`(rev, reverseComplement(.), .) }
  })) %>% 
  mutate(qname = map(reads, names),
         n_reads = lengths(qname))

reads <- do.call('c', bam_tbl$reads)

if (length(reads) > max_reads) {
  set.seed(seed)
  reads <- sample(reads, max_reads)
}
n_reads <- length(reads)

amp_names <- names(amplicons)
types <- c(A = amp_names[1],
           B = amp_names[2],
           AB = str_c(amp_names, collapse = '-'),
           BA = str_c(rev(amp_names), collapse = '-'))

# find fusion breakpoint window for all reads
chim_bp_tbl <-
  map_dfr(names(reads), function(i) {
    
    aln_fus_tbl <-
      c(ref_seq, setNames(reads[i], 'query')) %>% 
      AlignSeqs(verbose = FALSE) %>% 
      as.matrix() %>% 
      t() %>% 
      as_tibble() %>% 
      mutate(pos_A = cumsum(primary != '-'),
             pos_B = cumsum(secondary != '-')) %>% 
      filter(query != '-') %>% 
      mutate(pos_query = seq_len(n()),
             eq_A = query == primary,
             eq_B = query == secondary,
             A_f_cs = cumsum(eq_A),
             A_r_cs = rev(cumsum(rev(eq_A)) - rev(eq_A)),
             B_f_cs = cumsum(eq_B),
             B_r_cs = rev(cumsum(rev(eq_B)) - rev(eq_B)),
             A_B = A_f_cs + B_r_cs,
             B_A = B_f_cs + A_r_cs) %>% 
      select(eq_A, eq_B, starts_with('pos'), A_B, B_A)
    
    score_A <- sum(aln_fus_tbl$eq_A)
    score_B <- sum(aln_fus_tbl$eq_B)
    
    fus_res <-
      aln_fus_tbl %>% 
      mutate(max_score = pmax(B_A, A_B)) %>% 
      filter(max_score >= max(score_A, score_B) + min_score_delta,
             max_score == max(max_score)) %>%
      mutate(type = if_else(A_B > B_A, 'AB', 'BA')) %>% 
      select(pos_A, pos_query, type, max_score)
    
    if (nrow(fus_res)) {
      fus_res %>% 
        select(-type, -max_score) %>% 
        pivot_longer(everything(),
                     names_to = 'target',
                     names_prefix = 'pos_',
                     values_to = 'pos') %>% 
        arrange_all() %>% 
        group_by(target) %>% 
        summarise(start = first(pos),
                  end = last(pos),
                  .groups = 'drop') %>% 
        pivot_wider(names_from = target, values_from = c(start, end)) %>% 
        mutate(type = fus_res$type[1], score = fus_res$max_score[1], qname = i) %>% 
        select(qname, type, score, everything())
    } else if (score_A > score_B) {
      tibble(qname = i, type = 'A', score = score_A,
             start_A = NA_integer_, end_A = NA_integer_,
             start_query = NA_integer_, end_query = NA_integer_)
    } else {
      tibble(qname = i, type = 'B', score = score_B,
             start_A = NA_integer_, end_A = NA_integer_,
             start_query = NA_integer_, end_query = NA_integer_)
    }
  }) %>% 
  rename(start = start_A, end = end_A) %>%
  mutate(mid = round((start + end) / 2)) %>%
  left_join(bam_tbl %>% select(ref, amplicon, qname) %>% unnest(qname), by = "qname")

potential_fusions <-
  chim_bp_tbl %>% 
  select(qname, type, start, end, mid) %>% 
  filter(type %in% c('AB', 'BA')) %>% 
  nest(data = c(qname, start, end, mid)) %>% 
  mutate(data = map(data, function(data) {
      data %>% 
      select(pos = mid, qname) %>% 
      chop(qname) %>% 
      mutate(n = lengths(qname)) %>% 
      (function(x) {
        x %>% 
          mutate(pos = map(pos, ~ seq.int(. - fusion_window, . + fusion_window))) %>%
          unnest(pos) %>%
          group_by(pos) %>%
          summarise(n = sum(n),
                    qname = list(unlist(qname) %>% unique() %>% sort()),
                    .groups = 'drop') %>% 
          mutate(qhash = map_chr(qname, digest))
      }) %>% 
      filter(n == max(n)) %>% 
      mutate(block = cumsum(
        replace_na((pos != lag(pos) + 1) | qhash != lag(qhash), TRUE))) %>% 
      group_by(block) %>% 
      summarise(n = first(n),
                start = first(pos),
                end = last(pos),
                qname = qname[1],
                .groups = 'drop') %>% 
      unnest(qname) %>% 
      left_join(select(data, qname, mid), by = 'qname') %>% 
      mutate(dist_mid = abs(mid - (start+end)/2)) %>% 
      group_by(qname) %>% 
      slice(which.max(dist_mid)) %>% 
      group_by(block, start, end) %>% 
      summarise(n = n(),
                p_tot = n / n_reads,
                p_fus = n / nrow(data),
                qname = list(qname),
                .groups = 'drop') %>% 
      filter(n == max(n))
  })) %>% 
  {`if`(nrow(.),
        unnest(., data),
        tibble(type = character(),
               block = integer(),
               start = integer(),
               end = integer(),
               n = integer(),
               p_tot = double(),
               p_fus = double(),
               qname = list()))
  } %>% 
  mutate(pass = n > min_fusion_reads & p_tot > min_total_prop & p_fus > min_chim_prop,
         mid = round((start + end)/2) %>% as.integer()) %T>%
  { select(., -qname) %>% 
      mutate(type = types[type]) %>% 
      write_csv(str_c(out_prefix, '.candidates.csv')) }
  

# write outputs
chim_bp_tbl %>% 
  filter(!is.na(start)) %>% 
  select(qname, type, score, start, end, start_query, end_query) %>% 
  mutate(type = types[type]) %>% 
  write_csv(str_c(out_prefix, '.breakpoints.csv.gz'))
  

chim_bp_tbl %>%
  select(qname, ref, type, amplicon) %>%
  left_join(filter(potential_fusions, pass) %>% 
              select(type, breakpoint = mid, qname) %>% 
              unnest(qname),
            by = c("qname", "type")) %>% 
  chop(qname) %>% 
  mutate(type = case_when(
    type == ref             ~ str_c('clean-', amplicon),
    !is.na(breakpoint)      ~ str_c('fusion-', types[type], '-', breakpoint),
    type %in% c('AB', 'BA') ~ str_c('chimera-', types[type]),
    )) %>% 
  filter(!is.na(type)) %>%
  unchop(qname) %>% 
  select(type, ref, qname) %>% 
  nest(data = -type) %>% 
  mutate(bam_file = map2_chr(type, data, function(type, data) {
    dest_file <- str_c(out_prefix, type, 'bam', sep = '.')
    filter <- FilterRules(list(qn = function(read) {read$qname %in% data$qname}))
    if (n_distinct(data$ref) > 1) {
      tmp_bams <-
        filter(bam_tbl, ref %in% unique(data$ref)) %>% 
        pull(file) %>% 
        map_chr(function(file) {
          if (!file.exists(str_c(file, '.bai'))) {
            indexBam(file)
          }
          filterBam(file, tempfile(pattern = 'tmp_', tmpdir = '.', fileext = '.bam'), filter = filter)
        })
      mergeBam(tmp_bams, dest_file, overwrite = TRUE)
      file.remove(tmp_bams)
    } else {
      file <- filter(bam_tbl, ref == data$ref[1]) %>% pull(file)
      if (!file.exists(str_c(file, '.bai'))) {
        indexBam(file)
      }
      filterBam(file, dest_file, filter = filter)
    }
    return(dest_file)
  })) %>% 
  mutate(n_reads = map_int(data, nrow)) %>% 
  select(-data) %>% 
  write_csv(str_c(out_prefix, '.fus_smry.csv'))

