#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--bamlist",
        type = "character",
        help = "Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai"
    ), 
    make_option(
        "--region",
        type = "character",
        help = "HLA region to be analyzed, for example A for HLA-A"
    ), 
    make_option(
        "--dict_file",
        type = "character",
        help = "Path to dictionary file for reference build"
    ), 
    make_option(
        "--hla_gene_region_file",
        type = "character",
        help = "For reference packages built after QUILT 1.0.2, this is not used. For older reference packages, this is needed, and is a path to file with gene boundaries. 4 columns, named Name Chr Start End, with respectively gene name (e.g. HLA-A), chromsome (e.g. chr6), and 1 based start and end positions of gene [default NULL] ",
        default = NULL
    ), 
    make_option(
        "--outputdir",
        type = "character",
        help = "What output directory to use. Otherwise defaults to current directory [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--summary_output_file_prefix",
        type = "character",
        help = "Prefix for output text summary files [default 'quilt.hla.output'] ",
        default = 'quilt.hla.output'
    ), 
    make_option(
        "--nCores",
        type = "integer",
        help = "How many cores to use [default 1] ",
        default = 1
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
        "--final_output_RData_file",
        type = "character",
        help = "Final output RData file path, if desired [default NA] ",
        default = NA
    ), 
    make_option(
        "--write_summary_text_files",
        type = "logical",
        help = "Whether to write out final summary text files or not [default TRUE] ",
        default = TRUE
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
    ), 
    make_option(
        "--summary_best_alleles_threshold",
        type = "double",
        help = "When reporting results, give results until posterior probability exceeds this value [default 0.99] ",
        default = 0.99
    ), 
    make_option(
        "--downsampleToCov",
        type = "double",
        help = "For imputing states specifically using QUILT, what coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage. This is not used in the direct read mapping [default 30] ",
        default = 30
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT_HLA(
    bamlist = opt$bamlist,
    region = opt$region,
    dict_file = opt$dict_file,
    hla_gene_region_file = opt$hla_gene_region_file,
    outputdir = opt$outputdir,
    summary_output_file_prefix = opt$summary_output_file_prefix,
    nCores = opt$nCores,
    prepared_hla_reference_dir = opt$prepared_hla_reference_dir,
    quilt_hla_haplotype_panelfile = opt$quilt_hla_haplotype_panelfile,
    final_output_RData_file = opt$final_output_RData_file,
    write_summary_text_files = opt$write_summary_text_files,
    nGibbsSamples = opt$nGibbsSamples,
    n_seek_iterations = opt$n_seek_iterations,
    quilt_seed = opt$quilt_seed,
    chr = opt$chr,
    quilt_buffer = opt$quilt_buffer,
    quilt_bqFilter = opt$quilt_bqFilter,
    summary_best_alleles_threshold = opt$summary_best_alleles_threshold,
    downsampleToCov = opt$downsampleToCov
)
