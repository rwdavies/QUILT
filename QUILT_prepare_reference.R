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
        "--chr",
        type = "character",
        help = "What chromosome to run. Should match BAM headers"
    ), 
    make_option(
        "--nGen",
        type = "double",
        help = "Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure"
    ), 
    make_option(
        "--regionStart",
        type = "integer",
        help = "When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd [default NA] ",
        default = NA
    ), 
    make_option(
        "--regionEnd",
        type = "integer",
        help = "When running imputation, where to stop [default NA] ",
        default = NA
    ), 
    make_option(
        "--buffer",
        type = "integer",
        help = "Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd [default NA] ",
        default = NA
    ), 
    make_option(
        "--output_file",
        type = "character",
        help = "Path to output RData file containing prepared haplotypes (has default value that works with QUILT) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_vcf_file",
        type = "character",
        help = "Path to reference haplotype file in VCF format (values must be 0 or 1) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_haplotype_file",
        type = "character",
        help = "Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_legend_file",
        type = "character",
        help = "Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_sample_file",
        type = "character",
        help = "Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference_populations",
        type = "character",
        help = "Vector with character populations to include from reference_sample_file e.g. CHB, CHS [default NA] ",
        default = NA
    ), 
    make_option(
        "--reference_phred",
        type = "integer",
        help = "Phred scaled likelihood or an error of reference haplotype. Higher means more confidence in reference haplotype genotypes, lower means less confidence [default 30] ",
        default = 30
    ), 
    make_option(
        "--reference_exclude_samplelist_file",
        type = "character",
        help = "File with one column of samples to exclude from reference samples e.g. in validation, the samples you are imputing [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--region_exclude_file",
        type = "character",
        help = "File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--genetic_map_file",
        type = "character",
        help = "Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate) [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--nMaxDH",
        type = "integer",
        help = "Integer Maximum number of distinct haplotypes to store in reduced form. Recommended to keep as 2 ** N - 1 where N is an integer greater than 0 i.e. 255, 511, etc [default NA] ",
        default = NA
    ), 
    make_option(
        "--tempdir",
        type = "character",
        help = "What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/ [default NA] ",
        default = NA
    ), 
    make_option(
        "--make_fake_vcf_with_sites_list",
        type = "logical",
        help = "Whether to output a list of sites as a minimal VCF, for example to use with GATK 3 to genotype given sites [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--output_sites_filename",
        type = "character",
        help = "If make_fake_vcf_with_sites_list is TRUE, optional desired filename where to output sites VCF [default NA] ",
        default = NA
    ), 
    make_option(
        "--expRate",
        type = "double",
        help = "Expected recombination rate in cM/Mb [default 1] ",
        default = 1
    ), 
    make_option(
        "--maxRate",
        type = "double",
        help = "Maximum recomb rate cM/Mb [default 100] ",
        default = 100
    ), 
    make_option(
        "--minRate",
        type = "double",
        help = "Minimum recomb rate cM/Mb [default 0.1] ",
        default = 0.1
    ), 
    make_option(
        "--use_mspbwt",
        type = "logical",
        help = "Build mspbwt indices to be used in imputation [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--mspbwt_nindices",
        type = "integer",
        help = "How many mspbwt indices to build [default 4L] ",
        default = 4L
    ), 
    make_option(
        "--override_use_eMatDH_special_symbols",
        type = "integer",
        help = "Not for general use. If NA will choose version appropriately depending on whether a PBWT flavour is used. [default NA] ",
        default = NA
    ), 
    make_option(
        "--use_hapMatcherR",
        type = "logical",
        help = "Used for nMaxDH less than or equal to 255. Use R raw format to hold hapMatcherR. Lowers RAM use [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--impute_rare_common",
        type = "logical",
        help = "Whether to use common SNPs first for imputation, followed by a round of rare imputation [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--rare_af_threshold",
        type = "double",
        help = "Allele frequency threshold under which SNPs are considered rare, otherwise they are considered common. [default 0.001] ",
        default = 0.001
    ), 
    make_option(
        "--use_list_of_columns_of_A",
        type = "logical",
        help = "If when using mspbwt, use columns of A rather than the whole thing, to speed up this version [default TRUE] ",
        default = TRUE
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT_prepare_reference(
    outputdir = opt$outputdir,
    chr = opt$chr,
    nGen = opt$nGen,
    regionStart = opt$regionStart,
    regionEnd = opt$regionEnd,
    buffer = opt$buffer,
    output_file = opt$output_file,
    reference_vcf_file = opt$reference_vcf_file,
    reference_haplotype_file = opt$reference_haplotype_file,
    reference_legend_file = opt$reference_legend_file,
    reference_sample_file = opt$reference_sample_file,
    reference_populations = eval(parse(text=opt$reference_populations)),
    reference_phred = opt$reference_phred,
    reference_exclude_samplelist_file = opt$reference_exclude_samplelist_file,
    region_exclude_file = opt$region_exclude_file,
    genetic_map_file = opt$genetic_map_file,
    nMaxDH = opt$nMaxDH,
    tempdir = opt$tempdir,
    make_fake_vcf_with_sites_list = opt$make_fake_vcf_with_sites_list,
    output_sites_filename = opt$output_sites_filename,
    expRate = opt$expRate,
    maxRate = opt$maxRate,
    minRate = opt$minRate,
    use_mspbwt = opt$use_mspbwt,
    mspbwt_nindices = opt$mspbwt_nindices,
    override_use_eMatDH_special_symbols = opt$override_use_eMatDH_special_symbols,
    use_hapMatcherR = opt$use_hapMatcherR,
    impute_rare_common = opt$impute_rare_common,
    rare_af_threshold = opt$rare_af_threshold,
    use_list_of_columns_of_A = opt$use_list_of_columns_of_A
)
