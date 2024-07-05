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
        "--method",
        type = "character",
        help = "What method to run (diploid or nipt) [default diploid] ",
        default = "diploid"
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
        "--fflist",
        type = "character",
        help = "Path to file with fetal fraction values, one row per entry, in the same order as the bamlist [default \"\"] ",
        default = ""
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
        "--n_burn_in_seek_its",
        type = "integer",
        help = "How many iterations of the seek_its should be burn in. As an example, if n_seek_its is 3 and n_burn_in_seek_its is 2, then only the dosage from the final round is included. If n_seek_its is 4 and n_burn_in_seek_its is 2, then dosages from the last two rounds are used. Default value NA sets n_burn_in_seek_its to n_seek_its minus 1 [default NA] ",
        default = NA
    ), 
    make_option(
        "--Ksubset",
        type = "integer",
        help = "How many haplotypes to use in the faster Gibbs sampling [default 600] ",
        default = 600
    ), 
    make_option(
        "--Knew",
        type = "integer",
        help = "How many haplotypes to replace per-iteration after doing the full reference panel imputation [default 600] ",
        default = 600
    ), 
    make_option(
        "--K_top_matches",
        type = "integer",
        help = "How many top haplotypes to store in each grid site when looking for good matches in the full haplotype reference panel. Large values potentially bring in more haplotype diversity, but risk losing haplotypes that are good matches over shorter distances [default 5] ",
        default = 5
    ), 
    make_option(
        "--output_gt_phased_genotypes",
        type = "logical",
        help = "When TRUE, output GT entry contains phased genotypes (haplotypes). When FALSE, it is from the genotype posteriors, and masked when the maximum genotype posterior entry is less than 0.9 [default TRUE] ",
        default = TRUE
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
        "--prepared_reference_filename",
        type = "character",
        help = "Optional path to prepared RData file with reference objects. Can be used instead of outputdir to coordinate use of QUILT_prepare_reference and QUILT [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--save_prepared_reference",
        type = "logical",
        help = "If preparing reference as part of running QUILT, whether to save the prepared reference output file. Note that if the reference was already made using QUILT_prepare_reference, this is ignored [default FALSE] ",
        default = FALSE
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
        help = "Optional, only needed when using genfile or phasefile. File with positions of where to impute, lining up one-to-one with genfile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab> [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--genfile",
        type = "character",
        help = "Path to gen file with high coverage results. Empty for no genfile. If both genfile and phasefile are given, only phasefile is used, as genfile (unphased genotypes) is derivative to phasefile (phased genotypes). File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--phasefile",
        type = "character",
        help = "Path to phase file with truth phasing results. Empty for no phasefile. Supercedes genfile if both options given. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header. [default \"\"] ",
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
        "--make_plots_block_gibbs",
        type = "logical",
        help = "Whether to make some plots of per-sample imputation looking at how the block Gibbs is performing. This can be extremely slow so use for debugging or visualizing performance on one-off situations not for general runs [default FALSE] ",
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
        help = "Whether to store what read labels were used during the Gibbs samplings (i.e. whether reads were assigned to arbitrary labelled haplotype 1 or 2) [default FALSE] ",
        default = FALSE
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
    ), 
    make_option(
        "--estimate_bq_using_truth_read_labels",
        type = "logical",
        help = "When using phasefile with known truth haplotypes, infer truth read labels, and use them to infer the real base quality against the bam recorded base qualities [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--override_default_params_for_small_ref_panel",
        type = "logical",
        help = "When set to TRUE, then when using a smaller reference panel size (fewer haplotypes than Ksubset), parameter choices are reset appropriately. When set to FALSE, original values are used, which might crash QUILT [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--gamma_physically_closest_to",
        type = "integer",
        help = "For HLA imputation, the physical position closest to the centre of the gene [default NA] ",
        default = NA
    ), 
    make_option(
        "--seed",
        type = "integer",
        help = "The seed that controls random number generation. When NA, not used# [default NA] ",
        default = NA
    ), 
    make_option(
        "--hla_run",
        type = "logical",
        help = "Whether to use QUILT to generate posterior state probabilities as part of QUILT-HLA [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--downsampleToCov",
        type = "double",
        help = "What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage [default 30] ",
        default = 30
    ), 
    make_option(
        "--minGLValue",
        type = "double",
        help = "For non-Gibbs full imputation, minimum allowed value in haplotype gl, after normalization. In effect, becomes 1/minGLValue becomes maximum difference allowed between genotype likelihoods [default 1e-10] ",
        default = 1e-10
    ), 
    make_option(
        "--minimum_number_of_sample_reads",
        type = "integer",
        help = "Minimum number of sample reads a sample must have for imputation to proceed. Samples that have fewer reads than this will not be imputed in a given region and all output will be set to missing [default 2] ",
        default = 2
    ), 
    make_option(
        "--nGen",
        type = "double",
        help = "Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure [default NA] ",
        default = NA
    ), 
    make_option(
        "--reference_vcf_file",
        type = "character",
        help = "Path to reference VCF file with haplotypes, matching the reference haplotype and legend file [default \"\"] ",
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
        "--print_extra_timing_information",
        type = "logical",
        help = "Print extra timing information, i.e. how long sub-processes take, to better understand why things take as long as they do [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--small_ref_panel_block_gibbs_iterations",
        type = "character",
        help = "What iterations to perform block Gibbs sampling for the Gibbs sampler [default c(3, 6, 9)] ",
        default = "c(3, 6, 9)"
    ), 
    make_option(
        "--small_ref_panel_gibbs_iterations",
        type = "integer",
        help = "How many iterations to run the Gibbs sampler for each time it is run (i.e. how many full passes to run the Gibbs sampler over all the reads) [default 20] ",
        default = 20
    ), 
    make_option(
        "--plot_per_sample_likelihoods",
        type = "logical",
        help = "Plot per sample likelihoods i.e. the likelihood as the method progresses through the Gibbs sampling iterations [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--use_small_eHapsCurrent_tc",
        type = "logical",
        help = "For testing purposes only [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--mspbwtL",
        type = "integer",
        help = "How many neighouring haplotypes to scan up and down at each grid. [default 3] ",
        default = 3
    ), 
    make_option(
        "--mspbwtM",
        type = "integer",
        help = "Minimun long grids matches [default 1] ",
        default = 1
    ), 
    make_option(
        "--use_mspbwt",
        type = "logical",
        help = "Use msPBWT to select new haplotypes [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--mspbwt_nindices",
        type = "integer",
        help = "How many mspbwt indices to build [default 4L] ",
        default = 4L
    ), 
    make_option(
        "--use_splitreadgl",
        type = "logical",
        help = "Use split real GL in hap selection and imputation [default FALSE] ",
        default = FALSE
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
        "--shard_check_every_pair",
        type = "logical",
        help = "When using shard gibbs sampler, whether to check every pair of SNPs, or not [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--use_eigen",
        type = "logical",
        help = "Use eigen library for per haploid full li and stephens pass of full haplotype reference panel [default TRUE] ",
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
        help = "Allele frequency yhreshold under which SNPs are considered rare, otherwise they are considered common [default 0.001] ",
        default = 0.001
    ), 
    make_option(
        "--make_heuristic_plot",
        type = "logical",
        help = "Whether to make a plot for understanding heuristic performance [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--heuristic_approach",
        type = "character",
        help = "Which heuristic to use [default 'A'] ",
        default = 'A'
    ), 
    make_option(
        "--use_list_of_columns_of_A",
        type = "logical",
        help = "If when using mspbwt, use columns of A rather than the whole thing, to speed up this version [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--calculate_gamma_on_the_fly",
        type = "logical",
        help = "If when calculating genProbs, calculate gamma on the fly rather than saving [default TRUE] ",
        default = TRUE
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(QUILT))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
QUILT(
    outputdir = opt$outputdir,
    chr = opt$chr,
    method = opt$method,
    regionStart = opt$regionStart,
    regionEnd = opt$regionEnd,
    buffer = opt$buffer,
    fflist = opt$fflist,
    bamlist = opt$bamlist,
    cramlist = opt$cramlist,
    sampleNames_file = opt$sampleNames_file,
    reference = opt$reference,
    nCores = opt$nCores,
    nGibbsSamples = opt$nGibbsSamples,
    n_seek_its = opt$n_seek_its,
    n_burn_in_seek_its = opt$n_burn_in_seek_its,
    Ksubset = opt$Ksubset,
    Knew = opt$Knew,
    K_top_matches = opt$K_top_matches,
    output_gt_phased_genotypes = opt$output_gt_phased_genotypes,
    heuristic_match_thin = opt$heuristic_match_thin,
    output_filename = opt$output_filename,
    RData_objects_to_save = eval(parse(text=opt$RData_objects_to_save)),
    output_RData_filename = opt$output_RData_filename,
    prepared_reference_filename = opt$prepared_reference_filename,
    save_prepared_reference = opt$save_prepared_reference,
    tempdir = opt$tempdir,
    bqFilter = opt$bqFilter,
    panel_size = opt$panel_size,
    posfile = opt$posfile,
    genfile = opt$genfile,
    phasefile = opt$phasefile,
    maxDifferenceBetweenReads = opt$maxDifferenceBetweenReads,
    make_plots = opt$make_plots,
    make_plots_block_gibbs = opt$make_plots_block_gibbs,
    verbose = opt$verbose,
    shuffle_bin_radius = opt$shuffle_bin_radius,
    iSizeUpperLimit = opt$iSizeUpperLimit,
    record_read_label_usage = opt$record_read_label_usage,
    record_interim_dosages = opt$record_interim_dosages,
    use_bx_tag = opt$use_bx_tag,
    bxTagUpperLimit = opt$bxTagUpperLimit,
    addOptimalHapsToVCF = opt$addOptimalHapsToVCF,
    estimate_bq_using_truth_read_labels = opt$estimate_bq_using_truth_read_labels,
    override_default_params_for_small_ref_panel = opt$override_default_params_for_small_ref_panel,
    gamma_physically_closest_to = opt$gamma_physically_closest_to,
    seed = opt$seed,
    hla_run = opt$hla_run,
    downsampleToCov = opt$downsampleToCov,
    minGLValue = opt$minGLValue,
    minimum_number_of_sample_reads = opt$minimum_number_of_sample_reads,
    nGen = opt$nGen,
    reference_vcf_file = opt$reference_vcf_file,
    reference_haplotype_file = opt$reference_haplotype_file,
    reference_legend_file = opt$reference_legend_file,
    reference_sample_file = opt$reference_sample_file,
    reference_populations = opt$reference_populations,
    reference_phred = opt$reference_phred,
    reference_exclude_samplelist_file = opt$reference_exclude_samplelist_file,
    region_exclude_file = opt$region_exclude_file,
    genetic_map_file = opt$genetic_map_file,
    nMaxDH = opt$nMaxDH,
    make_fake_vcf_with_sites_list = opt$make_fake_vcf_with_sites_list,
    output_sites_filename = opt$output_sites_filename,
    expRate = opt$expRate,
    maxRate = opt$maxRate,
    minRate = opt$minRate,
    print_extra_timing_information = opt$print_extra_timing_information,
    small_ref_panel_block_gibbs_iterations = eval(parse(text=opt$small_ref_panel_block_gibbs_iterations)),
    small_ref_panel_gibbs_iterations = opt$small_ref_panel_gibbs_iterations,
    plot_per_sample_likelihoods = opt$plot_per_sample_likelihoods,
    use_small_eHapsCurrent_tc = opt$use_small_eHapsCurrent_tc,
    mspbwtL = opt$mspbwtL,
    mspbwtM = opt$mspbwtM,
    use_mspbwt = opt$use_mspbwt,
    mspbwt_nindices = opt$mspbwt_nindices,
    use_splitreadgl = opt$use_splitreadgl,
    override_use_eMatDH_special_symbols = opt$override_use_eMatDH_special_symbols,
    use_hapMatcherR = opt$use_hapMatcherR,
    shard_check_every_pair = opt$shard_check_every_pair,
    use_eigen = opt$use_eigen,
    impute_rare_common = opt$impute_rare_common,
    rare_af_threshold = opt$rare_af_threshold,
    make_heuristic_plot = opt$make_heuristic_plot,
    heuristic_approach = opt$heuristic_approach,
    use_list_of_columns_of_A = opt$use_list_of_columns_of_A,
    calculate_gamma_on_the_fly = opt$calculate_gamma_on_the_fly
)
