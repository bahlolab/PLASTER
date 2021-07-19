#!/usr/bin/env Rscript

stopifnot(require(tidyverse),
          require(docopt),
          require(AmpPhaseR))

doc <- "
Usage:
  AmpPhaseR_wrapper.R <vcf> <out_prefix> [options]

Options:
  vcf                         vcf(.gz) file with sample variants.
  out_prefix                  Prefix of output files.
  --copy-num=<f>              Copy number of sample [default: 2].
  --max-breakpoints=<f>       Maximum number of breakpoints for chimera detection [default: 2].
  --max-allele-freq=<f>       Maximum major allele frequency at a given variant position [default: 0.95].
  --min-phase-freq=<f>        Minimum frequency of a candidate phase [default: 0.05].
  --max-freq-abs-delta=<f>    Maximum absolute difference between observed and theoretical phase frequencies [default: 0.75].
  --max-freq-rel-delta=<f>    Maximum absolute difference between observed and theoretical phase frequencies relative to theoretical [default: 0.75].
"

opts <- docopt(doc)

copy_num <- as.integer(opts$copy_num)
max_freq_abs_delta <- as.numeric(opts$max_freq_abs_delta)
max_freq_rel_delta <- as.numeric(opts$max_freq_rel_delta)
max_breakpoints <- as.integer(opts$max_breakpoints)
max_allele_freq <- as.numeric(opts$max_allele_freq)
min_phase_freq <- as.numeric(opts$min_phase_freq)

stopifnot(file.exists(opts$vcf))

gds_fn <- stringr::str_replace(opts$vcf, '\\.vcf(\\.gz)?', '.gds')
if (! file.exists(gds_fn)) {
  SeqArray::seqVCF2GDS(vcf.fn = opts$vcf, out.fn = gds_fn, storage.option = 'ZIP_RA', ignore.chr.prefix = '')
}
gds <- SeqArray::seqOpen(gds_fn, allow.duplicate = TRUE)
name <- basename(opts$out_prefix)

result <- phase_amplicon(gds,
                         min_copy_num = copy_num,
                         max_copy_num = copy_num,
                         max_freq_abs_delta = max_freq_abs_delta,
                         max_freq_rel_delta = max_freq_rel_delta,
                         max_breakpoints = max_breakpoints,
                         max_allele_freq = max_allele_freq,
                         min_phase_freq = min_phase_freq,
                         name = name)
if (result$success) {
  result$read_phase %>% 
    count(phase) %>% 
    mutate(freq = n / sum(n)) %>% 
    arrange(desc(n)) %>% 
    write_tsv(str_c(opts$out_prefix, '_read_phase_smry.tsv.gz'))

  result$read_phase %>% 
    select(qname = read_name, phase) %>% 
    filter(str_detect(phase, 'phase_')) %>% 
    mutate(phase = str_remove(phase, 'phase_')) %>% 
    write_tsv(str_c(opts$out_prefix, '_phased.tsv.gz'))
  result$phase_summary %>%
    select(phase, count, freq, ratio, copy_num) %>%
    write_tsv(str_c(opts$out_prefix, '_phase_summary.tsv'))

  ggsave(plot = result$phase_plot, str_c(opts$out_prefix, '_phase_plot.png'), width =8, height = 11)
} else {
  # create empty files so nextflow process doesn't fail
  tibble(phase = character(), n = integer(), freq = double()) %>% 
    write_tsv(str_c(opts$out_prefix, '_read_phase_smry.tsv.gz'))

  tibble(qname = character(), phase = integer()) %>% 
    write_tsv(str_c(opts$out_prefix, '_phased.tsv.gz'))

  tibble(phase = character(), count = integer(), freq = double(),
         ratio = integer(), copy_num = integer()) %>%
    write_tsv(str_c(opts$out_prefix, '_phase_summary.tsv'))

  file.create(str_c(opts$out_prefix, '_phase_plot.png'))
}
