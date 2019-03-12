#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if(!require(rmarkdown)){
    install.packages("rmarkdown")
    library(rmarkdown)
}
if(!require(kableExtra)){
  install.packages("kableExtra")
  library(kableExtra)
}
if(!require(revealjs)){
    location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE)
    install.packages("revealjs", lib=location, repos='https://cran.csiro.au')
    library(revealjs, lib.loc=location)
}
rmarkdown::render(args[1], output_format = args[2])
