#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

## INSTALL CRAN PACKAGES
cran_pkgs <-
  unlist(
    lapply(args[grepl('^CRAN:', args)],
           function(x) scan(gsub('^CRAN:', '', x), what = character())))

if (length(cran_pkgs)) {
  install.packages(cran_pkgs,
                   repos='https://cloud.r-project.org',
                   clean = TRUE)
  stopifnot(all(cran_pkgs %in% rownames(installed.packages())))
}

## INSTALL BIOCONDUCTOR PACAKGES
bioc_packages <-
  unlist(
    lapply(args[grepl('^BIOC:', args)],
           function(x) scan(gsub('^BIOC:', '', x), what = character())))

if (length(bioc_packages)) {
  if (!'BiocManager' %in% rownames(installed.packages())){
    install.packages('BiocManager',
                     repos='https://cloud.r-project.org',
                     clean = TRUE)
  }
  BiocManager::install(bioc_packages)
  stopifnot(all(bioc_packages %in% rownames(installed.packages())))
}

## INSTALL GITHUB PACKAGES
github_packages <-
  unlist(
    lapply(args[grepl('^GITHUB:', args)],
           function(x) scan(gsub('^GITHUB:', '', x), what = character())))

if (length(github_packages)) {
  if (!'devtools' %in% rownames(installed.packages())){
    install.packages('devtools',
                     repos='https://cloud.r-project.org',
                     clean = TRUE)
  }
  devtools::install_github(github_packages, force = TRUE, upgrade = 'never')
  stopifnot(all(gsub('@.+$', '', gsub('^[^/]+/', '', github_packages)) %in%
                  rownames(installed.packages())))
}
