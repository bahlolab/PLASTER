
#' @importFrom magrittr set_rownames set_colnames "%>%"
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom tibble rownames_to_column column_to_rownames as_tibble tibble
#' @importFrom dplyr select mutate left_join case_when if_else starts_with count group_by ungroup summarise arrange
#' @importFrom purrr map map2
#' @importFrom stringr str_c str_split str_remove
allele_matrix <- function(gds) {

  # check gds has required data
  stopifnot(is_gds(gds),
            "genotype" %in% gdsfmt::ls.gdsn(gds, recursive = TRUE),
            "annotation/format/DP" %in% gdsfmt::ls.gdsn(gds, recursive = TRUE))

  gts <- SeqArray::seqGetData(gds, 'genotype')
  # check is haploid calls
  stopifnot(dim(gts[1]) == 1)
  gts <-
    abind::adrop(gts, 1) %>%
    set_rownames(SeqArray::seqGetData(gds, 'sample.id')) %>%
    set_colnames(str_c('V', SeqArray::seqGetData(gds, 'variant.id')))

  varinf <-
    SeqVarTools::variantInfo(gds) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    mutate(allele_data = map2(ref, alt, function(ref, alt) {
      `if`(nchar(alt) == 0,
           tibble(allele = ref, gt = 0L),
           tibble(allele = c(ref, str_split(alt, ',', simplify = T)),
                  gt = seq_along(allele) - 1L))
    }))

  depth <-
    SeqArray::seqGetData(gds, 'annotation/format/DP')$data %>%
    set_rownames(SeqArray::seqGetData(gds, 'sample.id')) %>%
    set_colnames(str_c('V', SeqArray::seqGetData(gds, 'variant.id'))) %>%
    as.data.frame() %>%
    rownames_to_column('read') %>%
    pivot_longer(starts_with('V'),
                 names_to = 'variant',
                 names_transform = list(variant = as.integer),
                 names_prefix = 'V',
                 values_to = 'dp')

  read_allele_tbl <-
    as.data.frame(gts) %>%
    rownames_to_column('read') %>%
    pivot_longer(starts_with('V'),
                 names_to = 'variant',
                 names_transform = list(variant = as.integer),
                 names_prefix = 'V',
                 values_to = 'gt') %>%
    left_join(varinf %>%
                select(variant = variant_id, allele_data) %>%
                unnest(allele_data),
              by = c("variant", "gt")) %>%
    left_join(depth,
              by = c("read", "variant")) %>%
    mutate(allele = case_when(dp == 0L ~ 'D',
                              gt == 0L ~'R',
                              TRUE ~ allele))

  read_allele_mat <-
    read_allele_tbl %>%
    select(read, variant, allele) %>%
    pivot_wider(names_from = variant,
                names_prefix = 'V',
                values_from = allele) %>%
    as.data.frame() %>%
    column_to_rownames('read') %>%
    as.matrix()

  return(read_allele_mat)
}

summarise_vars <- function(read_allele_mat) {

  stopifnot(is.matrix(read_allele_mat),
            !is.null(rownames(read_allele_mat)),
            !is.null(colnames(read_allele_mat)),
            all(dim(read_allele_mat) > 0))

  smry <-
    as.data.frame(read_allele_mat) %>%
    rownames_to_column('read') %>%
    pivot_longer(starts_with('V'),
                 names_to = 'variant',
                 values_to = 'allele') %>%
    count(variant, allele) %>%
    group_by(variant) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    arrange(as.integer(str_remove(variant, 'V')), allele)

  return(smry)
}
