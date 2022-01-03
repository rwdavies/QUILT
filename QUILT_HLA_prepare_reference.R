#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--outputdir",
        type = "character",
        help = "What output directory to use"
    ), 
    make_option(
        "--nGen",
        type = "double",
        help = "Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure, where K is number of haplotypes."
    ), 
    make_option(
        "--hla_types_panel",
        type = "character",
        help = "Path to file with 1000 Genomes formatted HLA types (see example for format details)"
    ), 
    make_option(
        "--ipd_igmt_alignments_zip_file",
        type = "character",
        help = "Path to zip file with alignments from IPD-IGMT (see README and example for more details)"
    ), 
    make_option(
        "--ref_fasta",
        type = "character",
        help = "Path to reference genome fasta"
    ), 
    make_option(
        "--refseq_table_file",
        type = "character",
        help = "Path to file with UCSC refseq gene information (see README and example for more details)"
    ), 
    make_option(
        "--full_regionStart",
        type = "double",
        help = "When building HLA full reference panel, start of maximal region spanning all HLA genes. The 1-based position x is kept if regionStart <= x <= regionEnd"
    ), 
    make_option(
        "--full_regionEnd",
        type = "double",
        help = "As above, but end of maximal region spanning all HLA genes"
    ), 
    make_option(
        "--buffer",
        type = "integer",
        help = "Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd"
    ), 
    make_option(
        "--region_exclude_file",
        type = "character",
        help = "File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference)"
    ), 
    make_option(
        "--genetic_map_file",
        type = "character",
        help = "Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_haplotype_file",
        type = "character",
        help = "Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)"
    ), 
    make_option(
        "--reference_legend_file",
        type = "character",
        help = "Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)"
    ), 
    make_option(
        "--reference_sample_file",
        type = "character",
        help = "Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)"
    ), 
    make_option(
        "--reference_exclude_samplelist_file",
        type = "character",
        help = "File with one column of samples to exclude from final reference panels. These samples can be removed at one of two points, depending on reference_exclude_samples_for_initial_phasing. If reference_exclude_samples_for_initial_phasing is FALSE, then these samples, if present, are used for the initial phasing and allele assignment. If TRUE, then they are removed immediately [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_exclude_samples_for_initial_phasing",
        type = "logical",
        help = "See help for reference_exclude_samplelist_file [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--all_hla_regions",
        type = "character",
        help = "Character vector of all HLA genes for which IPD-IMGT files are available for [default c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y')] ",
        default = "c('A','B','C','DMA','DMB','DOA','DOB','DPA1','DPA2','DPB1','DPB2','DQA1','DQA2','DQB1','DRA','DRB1','DRB3','DRB4','DRB5','E','F','G','HFE','H','J','K','L','MICA','MICB','N','P','S','TAP1','TAP2','T','U','V','W','Y')"
    ), 
    make_option(
        "--hla_regions_to_prepare",
        type = "character",
        help = "Character vector of HLA genes to prepare for imputation [default c('A','B','C','DQA1','DQB1','DRB1')] ",
        default = "c('A','B','C','DQA1','DQB1','DRB1')"
    ), 
    make_option(
        "--chr",
        type = "character",
        help = "What chromosome to run (probably chr6 or similar) [default 'chr6'] ",
        default = 'chr6'
    ), 
    make_option(
        "--minRate",
        type = "double",
        help = "Minimum recomb rate cM/Mb [default 0.1] ",
        default = 0.1
    ), 
    make_option(
        "--nCores",
        type = "integer",
        help = "How many cores to use [default 1] ",
        default = 1
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT_HLA_prepare_reference(
    outputdir = opt$outputdir,
    nGen = opt$nGen,
    hla_types_panel = opt$hla_types_panel,
    ipd_igmt_alignments_zip_file = opt$ipd_igmt_alignments_zip_file,
    ref_fasta = opt$ref_fasta,
    refseq_table_file = opt$refseq_table_file,
    full_regionStart = opt$full_regionStart,
    full_regionEnd = opt$full_regionEnd,
    buffer = opt$buffer,
    region_exclude_file = opt$region_exclude_file,
    genetic_map_file = opt$genetic_map_file,
    reference_haplotype_file = opt$reference_haplotype_file,
    reference_legend_file = opt$reference_legend_file,
    reference_sample_file = opt$reference_sample_file,
    reference_exclude_samplelist_file = opt$reference_exclude_samplelist_file,
    reference_exclude_samples_for_initial_phasing = opt$reference_exclude_samples_for_initial_phasing,
    all_hla_regions = eval(parse(text=opt$all_hla_regions)),
    hla_regions_to_prepare = eval(parse(text=opt$hla_regions_to_prepare)),
    chr = opt$chr,
    minRate = opt$minRate,
    nCores = opt$nCores
)
