#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--chr",
        type = "character",
        help = "What chromosome to run. Should match BAM headers"
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT(
    chr = opt$chr
)
