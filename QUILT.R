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
        "--bamlist",
        type = "character",
        help = "Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--cramlist",
        type = "character",
        help = "Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into QUILT [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--sampleNames_file",
        type = "character",
        help = "Optional, if not specified, sampleNames are taken from the SM tag in the header of the BAM / CRAM file. This argument is the path to file with sampleNames for samples. It is used directly to name samples in the order they appear in the bamlist / cramlist [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference",
        type = "character",
        help = "Path to reference fasta used for making cram files. Only required if cramlist is defined [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--nCores",
        type = "integer",
        help = "How many cores to use [default 1] ",
        default = 1
    ), 
    make_option(
        "--nGibbsSamples",
        type = "integer",
        help = "How many Gibbs samples to use [default 7] ",
        default = 7
    ), 
    make_option(
        "--n_seek_its",
        type = "integer",
        help = "How many iterations between first using current haplotypes to update read labels, and using current read labels to get new reference haplotypes, to perform [default 3] ",
        default = 3
    ), 
    make_option(
        "--Ksubset",
        type = "integer",
        help = "How many haplotypes to use in the faster Gibbs sampling [default 400] ",
        default = 400
    ), 
    make_option(
        "--Knew",
        type = "integer",
        help = "How many haplotypes to replace per-iteration after doing the full reference panel imputation [default 100] ",
        default = 100
    ), 
    make_option(
        "--K_top_matches",
        type = "integer",
        help = "How many top haplotypes to store in each grid site when looking for good matches in the full haplotype reference panel. Large values potentially bring in more haplotype diversity, but risk losing haplotypes that are good matches over shorter distances [default 5] ",
        default = 5
    ), 
    make_option(
        "--heuristic_match_thin",
        type = "double",
        help = "What fraction of grid sites to use when looking for good matches in the full haplotype reference panel. Smaller values run faster but potentially miss haplotypes [default 0.1] ",
        default = 0.1
    ), 
    make_option(
        "--output_filename",
        type = "character",
        help = "Override the default bgzip-VCF / bgen output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple QUILT runs are processing on the same region then they may over-write each others inputs and outputs. [default NULL] ",
        default = NULL
    ), 
    make_option(
        "--RData_objects_to_save",
        type = "character",
        help = "Can be used to name interim and misc results from imputation to save an an RData file. Default NULL means do not save such output [default NULL] ",
        default = NULL
    ), 
    make_option(
        "--output_RData_filename",
        type = "character",
        help = "Override the default location for miscellaneous outputs saved in RData format [default NULL] ",
        default = NULL
    ), 
    make_option(
        "--tempdir",
        type = "character",
        help = "What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/ [default NA] ",
        default = NA
    ), 
    make_option(
        "--bqFilter",
        type = "double",
        help = "Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used [default as.integer(17)] ",
        default = as.integer(17)
    ), 
    make_option(
        "--panel_size",
        type = "integer",
        help = "Integer number of reference haplotypes to use, set to NA to use all of them [default NA] ",
        default = NA
    ), 
    make_option(
        "--posfile",
        type = "character",
        help = "Optional, only needed when using genfile. File with positions of where to impute, lining up one-to-one with genfile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab> [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--genfile",
        type = "character",
        help = "Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--phasefile",
        type = "character",
        help = "Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header. For NIPT imputation, there are 3 columns, representing maternal transmitted, maternal untransmitted, and paternal transmitted [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--maxDifferenceBetweenReads",
        type = "double",
        help = "How much of a difference to allow the reads to make in the forward backward probability calculation. For example, if P(read | state 1)=1 and P(read | state 2)=1e-6, re-scale so that their ratio is this value. This helps prevent any individual read as having too much of an influence on state changes, helping prevent against influence by false positive SNPs [default 1e10] ",
        default = 1e10
    ), 
    make_option(
        "--make_plots",
        type = "logical",
        help = "Whether to make some plots of per-sample imputation. Especially nice when truth data. This is pretty slow though so useful more for debugging and understanding and visualizing performance [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--verbose",
        type = "logical",
        help = "whether to be more verbose when running [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--shuffle_bin_radius",
        type = "integer",
        help = "Parameter that controls how to detect ancestral haplotypes that are shuffled during EM for possible re-setting. If set (not NULL), then recombination rate is calculated around pairs of SNPs in window of twice this value, and those that exceed what should be the maximum (defined by nGen and maxRate) are checked for whether they are shuffled [default 5000] ",
        default = 5000
    ), 
    make_option(
        "--iSizeUpperLimit",
        type = "double",
        help = "Do not use reads with an insert size of more than this value [default 1e6] ",
        default = 1e6
    ), 
    make_option(
        "--record_read_label_usage",
        type = "logical",
        help = "Whether to store what read labels were used during the Gibbs samplings (i.e. whether reads were assigned to arbitrary labelled haplotype 1 or 2) [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--record_interim_dosages",
        type = "logical",
        help = "Whether to record interim dosages or not [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--use_bx_tag",
        type = "logical",
        help = "Whether to try and use BX tag in same to indicate that reads come from the same underlying molecule [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--bxTagUpperLimit",
        type = "integer",
        help = "When using BX tag, at what distance between reads to consider reads with the same BX tag to come from different molecules [default 50000] ",
        default = 50000
    ), 
    make_option(
        "--addOptimalHapsToVCF",
        type = "logical",
        help = "Whether to add optimal haplotypes to vcf when phasing information is present, where optimal is imputation done when read label origin is known [default FALSE] ",
        default = FALSE
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT(
    outputdir = opt$outputdir,
    chr = opt$chr,
    regionStart = opt$regionStart,
    regionEnd = opt$regionEnd,
    buffer = opt$buffer,
    bamlist = opt$bamlist,
    cramlist = opt$cramlist,
    sampleNames_file = opt$sampleNames_file,
    reference = opt$reference,
    nCores = opt$nCores,
    nGibbsSamples = opt$nGibbsSamples,
    n_seek_its = opt$n_seek_its,
    Ksubset = opt$Ksubset,
    Knew = opt$Knew,
    K_top_matches = opt$K_top_matches,
    heuristic_match_thin = opt$heuristic_match_thin,
    output_filename = opt$output_filename,
    RData_objects_to_save = eval(parse(text=opt$RData_objects_to_save)),
    output_RData_filename = opt$output_RData_filename,
    tempdir = opt$tempdir,
    bqFilter = opt$bqFilter,
    panel_size = opt$panel_size,
    posfile = opt$posfile,
    genfile = opt$genfile,
    phasefile = opt$phasefile,
    maxDifferenceBetweenReads = opt$maxDifferenceBetweenReads,
    make_plots = opt$make_plots,
    verbose = opt$verbose,
    shuffle_bin_radius = opt$shuffle_bin_radius,
    iSizeUpperLimit = opt$iSizeUpperLimit,
    record_read_label_usage = opt$record_read_label_usage,
    record_interim_dosages = opt$record_interim_dosages,
    use_bx_tag = opt$use_bx_tag,
    bxTagUpperLimit = opt$bxTagUpperLimit,
    addOptimalHapsToVCF = opt$addOptimalHapsToVCF
)
