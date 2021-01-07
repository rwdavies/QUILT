#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--bamfile",
        type = "character",
        help = "Path to bamfile to analyze"
    ), 
    make_option(
        "--region",
        type = "character",
        help = "HLA region to be analyzed, for example A for HLA-A"
    ), 
    make_option(
        "--outputdir",
        type = "character",
        help = "Optional, what output directory to use. Use only if not specifying finaloutputfile [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--prepared_hla_reference_dir",
        type = "character",
        help = "Output directory containing HLA reference material necessary for QUILT HLA [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--quilt_hla_haplotype_panelfile",
        type = "character",
        help = "Prepared HLA haplotype reference panel file [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--finaloutputfile",
        type = "character",
        help = "Final output file path [default NA] ",
        default = NA
    ), 
    make_option(
        "--nGibbsSamples",
        type = "integer",
        help = "How many QUILT Gibbs samples to perform [default 15] ",
        default = 15
    ), 
    make_option(
        "--n_seek_iterations",
        type = "integer",
        help = "How many seek iterations to use in QUILT Gibbs sampling [default 3] ",
        default = 3
    ), 
    make_option(
        "--quilt_seed",
        type = "integer",
        help = "When running QUILT Gibbs sampling, what seed to use, optionally [default NA] ",
        default = NA
    ), 
    make_option(
        "--chr",
        type = "character",
        help = "What chromosome, probably chr6 or maybe 6 [default 'chr6'] ",
        default = 'chr6'
    ), 
    make_option(
        "--quilt_buffer",
        type = "integer",
        help = "For QUILT Gibbs sampling, what buffer to include around center of gene [default 500000] ",
        default = 500000
    ), 
    make_option(
        "--quilt_bqFilter",
        type = "integer",
        help = "For QUILT Gibbs sampling, do not consider sequence information if the base quality is below this threshold [default 10] ",
        default = 10
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT_HLA(
    bamfile = opt$bamfile,
    region = opt$region,
    outputdir = opt$outputdir,
    prepared_hla_reference_dir = opt$prepared_hla_reference_dir,
    quilt_hla_haplotype_panelfile = opt$quilt_hla_haplotype_panelfile,
    finaloutputfile = opt$finaloutputfile,
    nGibbsSamples = opt$nGibbsSamples,
    n_seek_iterations = opt$n_seek_iterations,
    quilt_seed = opt$quilt_seed,
    chr = opt$chr,
    quilt_buffer = opt$quilt_buffer,
    quilt_bqFilter = opt$quilt_bqFilter
)
