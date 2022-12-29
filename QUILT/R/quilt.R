#' @title QUILT
#' @param outputdir What output directory to use
#' @param chr What chromosome to run. Should match BAM headers
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param bamlist Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai
#' @param cramlist Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into QUILT
#' @param sampleNames_file Optional, if not specified, sampleNames are taken from the SM tag in the header of the BAM / CRAM file. This argument is the path to file with sampleNames for samples. It is used directly to name samples in the order they appear in the bamlist / cramlist
#' @param reference Path to reference fasta used for making cram files. Only required if cramlist is defined
#' @param nCores How many cores to use
#' @param nGibbsSamples How many Gibbs samples to use
#' @param n_seek_its How many iterations between first using current haplotypes to update read labels, and using current read labels to get new reference haplotypes, to perform
#' @param n_burn_in_seek_its How many iterations of the seek_its should be burn in. As an example, if n_seek_its is 3 and n_burn_in_seek_its is 2, then only the dosage from the final round is included. If n_seek_its is 4 and n_burn_in_seek_its is 2, then dosages from the last two rounds are used. Default value NA sets n_burn_in_seek_its to n_seek_its minus 1
#' @param Ksubset How many haplotypes to use in the faster Gibbs sampling
#' @param Knew How many haplotypes to replace per-iteration after doing the full reference panel imputation
#' @param K_top_matches How many top haplotypes to store in each grid site when looking for good matches in the full haplotype reference panel. Large values potentially bring in more haplotype diversity, but risk losing haplotypes that are good matches over shorter distances
#' @param output_gt_phased_genotypes When TRUE, output GT entry contains phased genotypes (haplotypes). When FALSE, it is from the genotype posteriors, and masked when the maximum genotype posterior entry is less than 0.9
#' @param heuristic_match_thin What fraction of grid sites to use when looking for good matches in the full haplotype reference panel. Smaller values run faster but potentially miss haplotypes
#' @param output_filename Override the default bgzip-VCF / bgen output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple QUILT runs are processing on the same region then they may over-write each others inputs and outputs.
#' @param RData_objects_to_save Can be used to name interim and misc results from imputation to save an an RData file. Default NULL means do not save such output
#' @param output_RData_filename Override the default location for miscellaneous outputs saved in RData format
#' @param prepared_reference_filename Optional path to prepared RData file with reference objects. Can be used instead of outputdir to coordinate use of QUILT_prepare_reference and QUILT
#' @param save_prepared_reference If preparing reference as part of running QUILT, whether to save the prepared reference output file. Note that if the reference was already made using QUILT_prepare_reference, this is ignored
#' @param tempdir What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/
#' @param bqFilter Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used
#' @param panel_size Integer number of reference haplotypes to use, set to NA to use all of them
#'
#' @param posfile Optional, only needed when using genfile or phasefile. File with positions of where to impute, lining up one-to-one with genfile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>
#' @param genfile Path to gen file with high coverage results. Empty for no genfile. If both genfile and phasefile are given, only phasefile is used, as genfile (unphased genotypes) is derivative to phasefile (phased genotypes). File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header
#' @param phasefile Path to phase file with truth phasing results. Empty for no phasefile. Supercedes genfile if both options given. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header.
#' @param maxDifferenceBetweenReads How much of a difference to allow the reads to make in the forward backward probability calculation. For example, if P(read | state 1)=1 and P(read | state 2)=1e-6, re-scale so that their ratio is this value. This helps prevent any individual read as having too much of an influence on state changes, helping prevent against influence by false positive SNPs
#' @param make_plots Whether to make some plots of per-sample imputation. Especially nice when truth data. This is pretty slow though so useful more for debugging and understanding and visualizing performance
#' @param make_plots_block_gibbs Whether to make some plots of gibbs block sampling. This is pretty slow though so useful more for debugging and understanding and visualizing performance
#' @param verbose whether to be more verbose when running
#' @param shuffle_bin_radius Parameter that controls how to detect ancestral haplotypes that are shuffled during EM for possible re-setting. If set (not NULL), then recombination rate is calculated around pairs of SNPs in window of twice this value, and those that exceed what should be the maximum (defined by nGen and maxRate) are checked for whether they are shuffled
#' @param iSizeUpperLimit Do not use reads with an insert size of more than this value
#' @param record_read_label_usage Whether to store what read labels were used during the Gibbs samplings (i.e. whether reads were assigned to arbitrary labelled haplotype 1 or 2)
#' @param record_interim_dosages Whether to record interim dosages or not
#' @param use_bx_tag Whether to try and use BX tag in same to indicate that reads come from the same underlying molecule
#' @param bxTagUpperLimit When using BX tag, at what distance between reads to consider reads with the same BX tag to come from different molecules
#' @param addOptimalHapsToVCF Whether to add optimal haplotypes to vcf when phasing information is present, where optimal is imputation done when read label origin is known
#' @param estimate_bq_using_truth_read_labels When using phasefile with known truth haplotypes, infer truth read labels, and use them to infer the real base quality against the bam recorded base qualities
#' @param override_default_params_for_small_ref_panel When set to TRUE, then when using a smaller reference panel size (fewer haplotypes than Ksubset), parameter choices are reset appropriately. When set to FALSE, original values are used, which might crash QUILT
#' @param gamma_physically_closest_to For HLA imputation, the physical position closest to the centre of the gene
#' @param seed The seed that controls random number generation. When NA, not used#
#' @param hla_run Whether to use QUILT to generate posterior state probabilities as part of QUILT-HLA
#' @param downsampleToCov What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage
#' @param minGLValue For non-Gibbs full imputation, minimum allowed value in haplotype gl, after normalization. In effect, becomes 1/minGLValue becomes maximum difference allowed between genotype likelihoods
#' @param minimum_number_of_sample_reads Minimum number of sample reads a sample must have for imputation to proceed. Samples that have fewer reads than this will not be imputed in a given region and all output will be set to missing
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure
#' @param reference_vcf_file Path to reference VCF file with haplotypes, matching the reference haplotype and legend file
#' @param reference_haplotype_file Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)
#' @param reference_legend_file Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)
#' @param reference_sample_file Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)
#' @param reference_populations Vector with character populations to include from reference_sample_file e.g. CHB, CHS
#' @param reference_phred Phred scaled likelihood or an error of reference haplotype. Higher means more confidence in reference haplotype genotypes, lower means less confidence
#' @param reference_exclude_samplelist_file File with one column of samples to exclude from reference samples e.g. in validation, the samples you are imputing
#' @param region_exclude_file File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference)
#' @param genetic_map_file Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate)
#' @param nMaxDH Integer Maximum number of distinct haplotypes to store in reduced form. Recommended to keep as 2 ** N - 1 where N is an integer greater than 0 i.e. 255, 511, etc
#' @param make_fake_vcf_with_sites_list Whether to output a list of sites as a minimal VCF, for example to use with GATK 3 to genotype given sites
#' @param output_sites_filename If make_fake_vcf_with_sites_list is TRUE, optional desired filename where to output sites VCF
#' @param expRate Expected recombination rate in cM/Mb
#' @param maxRate Maximum recomb rate cM/Mb
#' @param minRate Minimum recomb rate cM/Mb
#' @param print_extra_timing_information Print extra timing information, i.e. how long sub-processes take, to better understand why things take as long as they do
#' @param small_ref_panel_block_gibbs_iterations What iterations to perform block Gibbs sampling for the Gibbs sampler
#' @param small_ref_panel_gibbs_iterations How many iterations to run the Gibbs sampler for each time it is run (i.e. how many full passes to run the Gibbs sampler over all the reads)
#'
#' @param plot_per_sample_likelihoods Plot per sample likelihoods i.e. the likelihood as the method progresses through the Gibbs sampling iterations
#' @param use_small_eHapsCurrent_tc For testing purposes only
#' @param pbwtL How many neighouring haplotypes to select forward and backwards at each grid. Automatically detected.
#' @param pbwtS How many grids as one step
#' @param zilong Using zilong's solution
#' @param use_mspbwt Use msPBWT to select new haplotypes
#' @param mspbwt_nindices How many mspbwt indices to build
#' @param use_splitreadgl Use split real GL in hap selection and imputation
#' @param use_eMatDH_special_symbols Whether to use RAM efficient version (not for general use)
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT <- function(
    outputdir,
    chr,
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    bamlist = "",
    cramlist = "",
    sampleNames_file = "",
    reference = "",
    nCores = 1,
    nGibbsSamples = 7,
    n_seek_its = 3,
    n_burn_in_seek_its = NA,
    Ksubset = 400,
    Knew = 400,
    K_top_matches = 5,
    output_gt_phased_genotypes = TRUE,
    heuristic_match_thin = 0.1,
    output_filename = NULL,
    RData_objects_to_save = NULL,
    output_RData_filename = NULL,
    prepared_reference_filename = "",
    save_prepared_reference = FALSE,
    tempdir = NA,
    bqFilter = as.integer(17),
    panel_size = NA,
    posfile = "",
    genfile = "",
    phasefile = "",
    maxDifferenceBetweenReads = 1e10,
    make_plots = FALSE,
    make_plots_block_gibbs = FALSE,
    verbose = TRUE,
    shuffle_bin_radius = 5000,
    iSizeUpperLimit = 1e6,
    record_read_label_usage = FALSE,
    record_interim_dosages = FALSE,
    use_bx_tag = TRUE,
    bxTagUpperLimit = 50000,
    addOptimalHapsToVCF = FALSE,
    estimate_bq_using_truth_read_labels = FALSE,
    override_default_params_for_small_ref_panel = TRUE,
    gamma_physically_closest_to = NA,
    seed = NA,
    hla_run = FALSE,
    downsampleToCov = 30,
    minGLValue = 1e-10,
    minimum_number_of_sample_reads = 2,
    nGen = NA,
    reference_vcf_file = "",
    reference_haplotype_file = "",
    reference_legend_file = "",
    reference_sample_file = "",
    reference_populations = NA,
    reference_phred = 30,
    reference_exclude_samplelist_file = "",
    region_exclude_file = "",
    genetic_map_file = "",
    nMaxDH = NA,
    make_fake_vcf_with_sites_list = FALSE,
    output_sites_filename = NA,
    expRate = 1,
    maxRate = 100,
    minRate = 0.1,
    print_extra_timing_information = FALSE,
    small_ref_panel_block_gibbs_iterations = c(3, 6, 9),
    small_ref_panel_gibbs_iterations = 20,
    plot_per_sample_likelihoods = FALSE,
    use_small_eHapsCurrent_tc = FALSE,
    pbwtL = 0,
    pbwtS = 1,
    zilong = FALSE,
    use_mspbwt = FALSE,
    mspbwt_nindices = 4L,
    use_splitreadgl = FALSE,
    use_eMatDH_special_symbols = TRUE
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

    ##
    use_sample_is_diploid <- TRUE
    ## turn this off for now
    ## plot_p1 Plot first haplotype read sampling probabilities
    plot_p1 = FALSE

    ## turn this off as well
    small_ref_panel_skip_equally_likely_reads <- FALSE
    small_ref_panel_equally_likely_reads_update_iterations <- c(1,2,3,6,9,15)


    ## re-label these internally
    ## n_gibbs_burn_in_its <- small_ref_panel_gibbs_iterations
    ## block_gibbs_iterations <- small_ref_panel_block_gibbs_iterations

    ## #' @param make_plots_block_gibbs Whether to make some plots of per-sample imputation looking at how the block Gibbs is performing. This can be extremely slow so use for debugging or visualizing performance on one-off situations not for general runs

    options(digits.secs=6)
    options(scipen = 999)

    ## need for outputting
    options(scipen = 999) ## Dangit rounding
    regionName <- chr
    if (!is.na(regionStart) & !is.na(regionEnd)) {
        regionName <- paste0(chr, ".", regionStart,".", regionEnd)
    }

    if (make_plots | make_plots_block_gibbs | plot_per_sample_likelihoods) {
        dir.create(file.path(outputdir, "plots"), showWarnings = FALSE, recursive = TRUE)
    }

    ##
    ## validate parameters
    ##
    STITCH::validate_tempdir(tempdir)
    check_program_dependency("bgzip")
    check_program_dependency("tabix")


    ##
    ## local validate
    ##
    validate_panel_size(panel_size)
    validate_minimum_number_of_sample_reads(minimum_number_of_sample_reads)
    validate_niterations_and_small_ref_panel_block_gibbs(small_ref_panel_block_gibbs_iterations, small_ref_panel_gibbs_iterations)

    if (is.na(n_burn_in_seek_its)) {
        n_burn_in_seek_its <- n_seek_its - 1
        print_message(paste0("Auto-set n_burn_in_seek_its to ", n_burn_in_seek_its, " i.e. only sample one dosage per Gibbs sample"))
    }
    validate_n_seek_its_and_n_burn_in_seek_its(n_seek_its, n_burn_in_seek_its)

    ## if (make_plots && phasefile == "") {
    ##     stop("If you want to make plots using make_plots, you need to provide phase information using phasefile")
    ## }

    if (is.na(tempdir)) {
        ## tempdir <- tempdir()
        tempdir <- tempfile()
        dir.create(tempdir)
    }
    if (!is.null(output_filename)) {
        if (!dir.exists(dirname(output_filename))) {
            dir.create(dirname(output_filename), showWarnings = FALSE)
        }
    }
    output_filename <- STITCH::get_output_filename(
        output_filename = output_filename,
        outputdir = outputdir,
        regionName = regionName,
        output_format = "bgvcf",
        prefix = "quilt"
    )

    if (is.null(output_RData_filename)) {
        output_RData_filename <- file_quilt_output_RData(outputdir, regionName)
    }

    if (!is.na(as.integer(seed))) {
        print_message(paste0("Setting seed with seed:", seed))
        set.seed(seed)
    }

    ## getting some weird non-reproducible problems on the cluster
    ## hopefully this helps sort out what's going on
    if (!dir.exists(tempdir)) {
        dir.create(tempdir, showWarnings = FALSE)
        tempfile <- tempfile(tmpdir = tempdir)
        cat("test", file = tempfile)
        unlink(tempfile)
    }


    ##print(args)
    ##print(paste0(commandArgs(trailingOnly = TRUE), collapse = "', '"))
    ## chr <- args[1]
    ##regionStart <- as.numeric(args[2])
    ##regionEnd <- as.numeric(args[3])
    ##buffer <- as.numeric(args[4])
    ##scenario <- args[5]
    ## inbams <- args[6]
    ## inref <- args[7]
    ## result_file <- args[8]
    ## nCores <- as.integer(args[9])
    ##genfile <- args[10]
    ##depth <- args[11]
    ##down_type <- args[12]
    ##record_interim_dosages <- as.logical(args[13])
    ## bqFilter <- as.numeric(args[14]) ## e.g. 17 or 10
    ##have_truth_haplotypes <- as.logical(args[15])
    ## panel_size <- args[18]

    ## hacky but OK!
    ##if (nGibbsSamples == 10) {
    ##    record_interim_dosages <- TRUE
    ## }


    ##
    ## load quilt prepared filename
    ##
    ## if we are building this here, assume we don't want to save the prepared output
    ## this also removes the problem of a collision if two jobs use the same files with the same output directories from different locations
    ##
    if (prepared_reference_filename == "") {
        prepared_reference_filename <- file_quilt_prepared_reference(outputdir, regionName)
    }
    if (!file.exists(prepared_reference_filename)) {
        if (!reference_haplotype_file == "") {
            if (!save_prepared_reference) {
                prepared_reference_filename <- tempfile(fileext = ".RData")
            }
            QUILT_prepare_reference(
                outputdir = outputdir,
                chr = chr,
                nGen = nGen,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                reference_haplotype_file = reference_haplotype_file,
                reference_legend_file = reference_legend_file,
                reference_sample_file = reference_sample_file,
                reference_populations = reference_populations,
                reference_phred = reference_phred,
                reference_exclude_samplelist_file = reference_exclude_samplelist_file,
                region_exclude_file = region_exclude_file,
                genetic_map_file = genetic_map_file,
                nMaxDH = nMaxDH,
                tempdir = tempdir,
                make_fake_vcf_with_sites_list = make_fake_vcf_with_sites_list,
                output_sites_filename = output_sites_filename,
                expRate = expRate,
                maxRate = maxRate,
                minRate = minRate,
                use_mspbwt = use_mspbwt,
                use_pbwt_index = zilong,
                mspbwt_nindices =  mspbwt_nindices,
                reference_vcf_file = reference_vcf_file,
                output_file = prepared_reference_filename,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols
            )
        } else {
            stop(paste0("Cannot find prepared haplotype reference file, expecting:", prepared_reference_filename))
        }
    }




    ## always check regionStart, regionEnd and buffer, just in case
    ## require them to be the same to prevent problems
    new_regionStart <- regionStart
    new_regionEnd <- regionEnd
    new_buffer <- buffer

    load(prepared_reference_filename)

    if (use_mspbwt) {
        if (is.null(ms_indices)) {
            stop("To use mspbwt and QUILT, you must prepare the reference package using use_mspbt=TRUE")
        }
    }

    if (zilong && use_mspbwt) {
        stop("Please select only one of zilong or use_mspbwt")
    }

    ## now build PBWT using nSNPs, nrow(rhb_t) loaded from last step
    if(zilong && reference_vcf_file == "") stop("Zilong requires VCF file. Please feed vcf file!")

    if (zilong) {
        if (is.null(zilong_indices)) {
            stop("To use zilong pbwt and QUILT, you must prepare the reference package using use_pbwt_index=TRUE")
        }
    }


    if (hla_run) {

        ## print_message("SIMON HLA CODE - fix this eventually!")
        which_hapMatcher_0 <- which(hapMatcher == 0, arr.ind = TRUE) - 1
        eMatDH_special_grid_which <- integer(nGrids)
        special_grids <- unique(which_hapMatcher_0[, 2]) + 1 ## this-is-1-based
        eMatDH_special_grid_which[special_grids] <- as.integer(1:length(special_grids))
        if (nrow(which_hapMatcher_0) > 0) {
            ## now build list with them
            x <- which_hapMatcher_0[, 2]
            y <- which((x[-1] - x[-length(x)]) > 0) ## last entry that is OK
            starts <- c(1, y + 1)
            ends <- c(y, length(x))
            ##
            ## eMatDH_special_values
            ##   list of length the number of special grids
            ##   entries are which ones to re-do, and where they are in rhb_t
            ##   entries inside this are 0-based
            eMatDH_special_values_list <- lapply(1:length(starts), function(i) {
                return(as.integer(which_hapMatcher_0[starts[i]:ends[i], 1]))
            })
        } else {
            eMatDH_special_values_list <- list()
        }
        nrow_which_hapMatcher_0 <- nrow(which_hapMatcher_0) ## for testing
    }

    validate_quilt_use_of_region_variables(
        regionStart,
        new_regionStart,
        regionEnd,
        new_regionEnd,
        buffer,
        new_buffer
    )


    ##
    ## possibly reset values
    ##
    if (override_default_params_for_small_ref_panel) {
        K <- nrow(rhb_t)
        if (K < Ksubset) {
            print_message("Overriding default parameters for small reference panel")
            print_message(paste0("Observing number of reference haplotypes K=", K))
            print_message(paste0("Reset n_seek_its from ", n_seek_its, " to ", 1))
            n_seek_its <- 1
            print_message(paste0("Reset n_burn_in_seek_its from ", n_burn_in_seek_its, " to ", 0))
            n_burn_in_seek_its <- 0
            print_message(paste0("Set Ksubset to ", K, " (no longer necessary)"))
            Ksubset <- K
            print_message(paste0("Set Knew to ", K, " (no longer necessary)"))
            Knew <- K
        }
    }


    ##
    ## optionally load genotypes and phasevali
    ##
    if ((genfile != "") | (phasefile != "")) {
        if ((genfile != "") & (posfile == "")) {
            stop("You have given genfile, a set of high confidence genotypes, but not posfile, which confirms their physical positions. Please fix this, and see the help files for more information")
        }
        if ((phasefile != "") & (posfile == "")) {
            stop("You have given phasefile, with high confidence phased genotypes, but not posfile, which confirms their physical positions. Please fix this, and see the help files for more information")
        }
        out <- get_and_validate_pos_gen_and_phase(
            posfile = posfile,
            genfile = genfile,
            phasefile = phasefile,
            chr = chr,
            verbose = TRUE
        )
        posX <- out$pos
        genX <- out$gen
        phaseX <- out$phase
        ##
        ## shrink (if regionStart and regionEnd are NA)
        ##
        out <- shrink_region(
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            L = posX[, 2],
            pos = posX,
            gen = genX,
            phase = phaseX
        )
        posX <- out$pos
        gen <- out$gen
        phase <- out$phase
        key1 <- paste0(pos[, 1], pos[, 2], pos[, 3], pos[, 4])
        key2 <- paste0(posX[, 1], posX[, 2], posX[, 3], posX[, 4])
        err <- paste0("There was an error lining up the reference data from the reference legend file with the posfile. Please double check that posfile and genfile/phasefile are defined as exactly the same SNPs in the region to be imputed as the reference legend file")
        if (length(key1) != length(key2)) {
            stop(err)
        }
        if (sum(key1 != key2) > 0) {
            stop(err)
        }
    } else {
        gen <- NULL
        phase <- NULL
    }



    ##
    ## can request smaller panel size, only really useful for testing speed
    ##
    if (!is.na(panel_size)) {
        if (use_mspbwt | use_zilong) {
            stop("This functionality does not work when using mspwt or zilong. Sorry!")
        }
        rhb_t <- rhb_t[1:as.integer(panel_size), ] ## this is the number of HAPLOTYPES
        reference_samples <- reference_samples[1:as.integer(panel_size), ]
        ##nMaxDH <- 2 ** 8 - 1
        out <- make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]
        hapMatcher <- out[["hapMatcher"]]
        eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
        eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]

    }
    K <- nrow(rhb_t)
    ancAlleleFreqAll <- ref_alleleCount[, 3]


    ##
    ## get sample names
    ##
    out <- get_sample_names(
        bamlist = bamlist,
        cramlist = cramlist,
        nCores = nCores,
        outputdir = outputdir,
        regionName = regionName,
        originalRegionName = originalRegionName,
        sampleNames_file = sampleNames_file,
        save = FALSE,
        duplicate_name_behaviour = "warn"
    )
    N <- out$N
    sampleNames <- out$sampleNames
    bam_files <- out$bam_files
    cram_files <- out$cram_files

    ##
    ## check line up between
    ##
    out <- match_gen_and_phase_to_samples(
        sampleNames = sampleNames,
        gen = gen,
        phase = phase
    )
    highCovInLow <- out$highCovInLow
    samples_with_phase <- out$samples_with_phase
    if (length(samples_with_phase) > 0) {
        have_truth_haplotypes <- TRUE
    } else {
        have_truth_haplotypes <- FALSE
    }


    ##
    ## work on various recombination rates
    ##
    print_message(paste0("There are ", nrow(pos), " SNPs in this region"))
    small_transMatRate_tc_H <- get_transMatRate_m("pseudoHaploid", sigmaCurrent_m)
    X <- get_transMatRate_m("pseudoHaploid", sigmaCurrent_m)
    full_transMatRate_t_H <- array(NA, c(dim(X)[1], dim(X)[2]))
    ## ARGH R dropping dimensions
    full_transMatRate_t_H[, ] <- X[, , 1]
    rate2 <- -log(small_transMatRate_tc_H[1, , 1]) * 100
    smooth_cm <- rcpp_make_smoothed_rate(rate2, L_grid, shuffle_bin_radius, verbose = FALSE);
    smooth_cm <- smooth_cm / max(smooth_cm)


    ##
    ## work with truth haplotypes for Optimal output
    ##
    ## if (have_truth_haplotypes) {
    ##     ## I also here want to record everything basically
    ##     ## for gibbs and seek checks
    ##     record_read_label_usage <- TRUE
    ##     haps <- fread(cmd = paste0("gunzip -c genotypes/phase.", scenario, ".", chr, ".", regionStart, ".", regionEnd, ".phased.vcf.gz | grep -v '##'"), data.table = FALSE)
    ##     ## is basically a VCF
    ##     truth_haps_all <- haps[match(L, haps[, 2]), ]
    ##     if (scenario == "NA12878") {
    ##         truth_gen <- cbind(truth_gen, NA12878ONT = truth_gen[, "NA12878"])
    ##         truth_haps_all <- cbind(truth_haps_all, NA12878ONT = truth_haps_all[, "NA12878"])
    ##     }
    ## }

    ## sampleNames <- get_sample_names_from_bam_or_cram_files(
    ##     bam_files,
    ##     nCores = 1,
    ##     file_type = "BAM",
    ##     verbose = FALSE
    ## )

    ## if (scenario == "ont") {
    ##     ## fuck it, hack this for now
    ##     for(sampleName in sampleNames[-grep("10X", sampleNames)]) {
    ##         truth_gen <- cbind(truth_gen[, sampleName], truth_gen)
    ##         colnames(truth_gen)[1] <- paste0(sampleName, "10X")
    ##         truth_haps_all <- cbind(truth_haps_all[, sampleName], truth_haps_all)
    ##         colnames(truth_haps_all)[1] <- paste0(sampleName, "10X")
    ##     }
    ## }


    ## ##
    ## ## check all are OK
    ## ##
    ## if (have_truth_haplotypes) {
    ##     x <- match(sampleNames, colnames(truth_haps_all))
    ##     if (sum(is.na(x)) > 0) {
    ##         print(x)
    ##         print(sampleNames)
    ##         print(colnames(truth_haps_all))
    ##         stop("cannot match all sampleNames to truth phasing columns")
    ##     }
    ## }


    ##
    ## here initialize where to start and stop getting reads from
    ##
    out <- initialize_chrStart_and_chrEnd(
        chrStart = NA,
        chrEnd = NA,
        L = L,
        iSizeUpperLimit = iSizeUpperLimit
    )
    chrStart <- out$chrStart
    chrEnd <- out$chrEnd
    chrLength <- quilt_get_chromosome_length(iBam = 1, bam_files = bam_files, cram_files = cram_files, chr = chr)
    if (chrEnd > chrLength) {
        chrEnd <- chrLength
    }



    ##
    ## mclapply the runs!
    ##
    iCore <- 1
    sampleRanges <- getSampleRange(N, nCores)
    complete_set_of_results <- mclapply(1:length(sampleRanges), mc.cores = nCores, function(iCore) {

        sampleRange <- sampleRanges[[iCore]]

        ## for output
        hweCount <- array(0, c(nSNPs, 3))
        infoCount <- array(0, c(nSNPs, 2))
        afCount <- array(0, nSNPs)
        alleleCount <- array(0, c(nSNPs, 2))

        K <- nrow(rhb_t)
        full_alphaHat_t <- array(0, c(K, nGrids))
        ## full_betaHat_t <- array(0, c(K, nGrids))
        full_betaHat_t <- array(0, c(1, 1))
        if (make_plots) {
            full_gamma_t <- array(0, c(K, nGrids))
        } else {
            full_gamma_t <- array(0, c(1, 1))
        }
        iSample <- 1
        ww <- seq(1, nGrids, length.out = max(1, round(heuristic_match_thin * nGrids)))
        full_gammaSmall_cols_to_get <- array(-1, nGrids)
        full_gammaSmall_cols_to_get[ww] <- 0:(length(ww) - 1)
        full_gammaSmall_t <-  array(0, c(K, length(ww)))
        ## hmm, OK
        K <- Ksubset
        S <- 1
        alphaHat_t1 <- array(0, c(K, nGrids))
        betaHat_t1 <- array(0, c(K, nGrids))
        eMatGrid_t1 <- array(0, c(K, nGrids))
        alphaHat_t2 <- array(0, c(K, nGrids))
        betaHat_t2 <- array(0, c(K, nGrids))
        eMatGrid_t2 <- array(0, c(K, nGrids))
        alphaHat_t3 <- array(0, c(K, nGrids))
        betaHat_t3 <- array(0, c(K, nGrids))
        eMatGrid_t3 <- array(0, c(K, nGrids))
        gammaMT_t_local <- array(0, c(K, nGrids))
        gammaMU_t_local <- array(0, c(K, nGrids))
        gammaP_t_local <- array(0, c(K, nGrids))
        ##
        small_priorCurrent_m <- array(1 / K, c(K, S))
        small_alphaMatCurrent_tc <- array(1 / K, c(K, nGrids - 1, S))
        if (use_small_eHapsCurrent_tc) {
            small_eHapsCurrent_tc <- array(0, c(K, nSNPs, S))
        } else {
            small_eHapsCurrent_tc <- array(0, c(1, 1, 1))
        }




        results_across_samples <- as.list(sampleRange[2] - sampleRange[1] + 1)

        for(iSample in sampleRange[1]:sampleRange[2]) {

            print_message(paste0("Imputing sample: ", iSample))
            out <- get_and_impute_one_sample(
                rhb_t = rhb_t,
                outputdir = outputdir,
                nGibbsSamples = nGibbsSamples,
                n_seek_its = n_seek_its,
                n_burn_in_seek_its = n_burn_in_seek_its,
                full_alphaHat_t = full_alphaHat_t,
                full_betaHat_t = full_betaHat_t,
                full_gamma_t = full_gamma_t,
                full_gammaSmall_t = full_gammaSmall_t,
                full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                full_transMatRate_t_H = full_transMatRate_t_H,
                small_transMatRate_tc_H = small_transMatRate_tc_H,
                alphaHat_t1 = alphaHat_t1,
                betaHat_t1 = betaHat_t1,
                eMatGrid_t1 = eMatGrid_t1,
                alphaHat_t2 = alphaHat_t2,
                betaHat_t2 = betaHat_t2,
                eMatGrid_t2 = eMatGrid_t2,
                alphaHat_t3 = alphaHat_t3,
                betaHat_t3 = betaHat_t3,
                eMatGrid_t3 = eMatGrid_t3,
                gammaMT_t_local = gammaMT_t_local,
                gammaMU_t_local = gammaMU_t_local,
                gammaP_t_local = gammaP_t_local,
                small_alphaMatCurrent_tc = small_alphaMatCurrent_tc,
                small_priorCurrent_m = small_priorCurrent_m,
                small_eHapsCurrent_tc = small_eHapsCurrent_tc,
                bam_files = bam_files,
                L = L,
                pos = pos,
                chr = chr,
                tempdir = tempdir,
                regionName = regionName,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                verbose = verbose,
                gen = gen,
                phase = phase,
                iSample = iSample,
                grid = grid,
                ancAlleleFreqAll = ancAlleleFreqAll,
                L_grid = L_grid,
                shuffle_bin_radius = shuffle_bin_radius,
                Ksubset = Ksubset,
                Knew = Knew,
                K_top_matches = K_top_matches,
                heuristic_match_thin = heuristic_match_thin,
                record_interim_dosages = record_interim_dosages,
                have_truth_haplotypes = have_truth_haplotypes,
                bqFilter = bqFilter,
                record_read_label_usage = record_read_label_usage,
                sampleNames = sampleNames,
                smooth_cm = smooth_cm,
                iSizeUpperLimit = iSizeUpperLimit,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                make_plots = make_plots,
                ref_error = ref_error,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                hapMatcher = hapMatcher,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                eMatDH_special_grid_which = eMatDH_special_grid_which,
                eMatDH_special_values_list = eMatDH_special_values_list,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                af = af,
                use_bx_tag = use_bx_tag,
                bxTagUpperLimit = bxTagUpperLimit,
                addOptimalHapsToVCF = addOptimalHapsToVCF,
                make_plots_block_gibbs = make_plots_block_gibbs,
                estimate_bq_using_truth_read_labels = estimate_bq_using_truth_read_labels,
                chrStart = chrStart,
                chrEnd = chrEnd,
                gamma_physically_closest_to = gamma_physically_closest_to,
                hla_run = hla_run,
                downsampleToCov = downsampleToCov,
                minGLValue = minGLValue,
                minimum_number_of_sample_reads = minimum_number_of_sample_reads,
                print_extra_timing_information = print_extra_timing_information,
                small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
                small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                plot_per_sample_likelihoods = plot_per_sample_likelihoods,
                use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
                output_gt_phased_genotypes = output_gt_phased_genotypes,
                pbwtL = pbwtL,
                pbwtS = pbwtS,
                zilong = zilong,
                zilong_indices =  zilong_indices,
                use_mspbwt = use_mspbwt,
                ms_indices = ms_indices,
                use_splitreadgl = use_splitreadgl,
                use_sample_is_diploid = use_sample_is_diploid,
                plot_p1 = plot_p1,
                small_ref_panel_skip_equally_likely_reads = small_ref_panel_skip_equally_likely_reads,
                small_ref_panel_equally_likely_reads_update_iterations = small_ref_panel_equally_likely_reads_update_iterations
            )

            if (out[["sample_was_imputed"]]) {
                ## for summarization
                infoCount[, 1] <- infoCount[, 1, drop = FALSE] + out[["eij"]]
                infoCount[, 2] <- infoCount[, 2, drop = FALSE] + (out[["fij"]] - out[["eij"]]**2)
                afCount <- afCount + (out[["eij"]]) / 2
                hweCount[out[["max_gen"]]] <- hweCount[out[["max_gen"]]] + 1 ## hmmmmm not ideal
                alleleCount <- alleleCount + out[["per_sample_alleleCount"]]
                ## drop now - not useful anymore
                out[["eij"]] <- NULL
                out[["fij"]] <- NULL
                out[["per_sample_alleleCount"]] <- NULL
                out[["max_gen"]] <- NULL
            }

            results_across_samples[[iSample - sampleRange[1] + 1]] <- out

            rm(out)

            ## optionally, do some gc here, if longer running job
            if (as.numeric(K) * as.numeric(nSNPs) > (1e6)) {
                ## print("temporary")
                ## print(head(sort( sapply(ls(),function(x){object.size(get(x))}), decreasing = TRUE)))
                ## print(object.size(results_across_samples))
                ## print(gc(reset = TRUE))
                for(i in 1:5) {
                    gc(reset = TRUE)
                }
            }

        }

        return(
            list(
                results_across_samples = results_across_samples,
                infoCount = infoCount,
                afCount = afCount,
                hweCount = hweCount,
                alleleCount = alleleCount
            )
        )

    })

    check_mclapply_OK(complete_set_of_results)

    ##
    ## make and write output file
    ##
    make_and_write_output_file(
        output_filename = output_filename,
        sampleNames = sampleNames,
        nSNPs = nSNPs,
        N = N,
        pos = pos,
        nCores = nCores,
        complete_set_of_results = complete_set_of_results,
        inRegion2 = inRegion2,
        sampleRanges = sampleRanges,
        addOptimalHapsToVCF = addOptimalHapsToVCF,
        output_gt_phased_genotypes = output_gt_phased_genotypes
    )


    ##
    ## build a singular set of results
    ##
    if (!is.null(RData_objects_to_save)) {

        print_message("Begin saving extra RData objects to disk")

        final_set_of_results <- as.list(1:N)
        c <- 1
        for(i in 1:length(complete_set_of_results)) {
            x <- complete_set_of_results[[i]]
            for(j in 1:length(x[["results_across_samples"]])) {
                final_set_of_results[[c]] <- x[["results_across_samples"]][[j]]
                c <- c + 1
            }
        }

        ## these are properly in the VCF
        ## still could be exported
        ## imputed_dosages <- array(NA, c(nrow(pos), length(final_set_of_results)))
        ## for(i in 1:length(final_set_of_results)) {
        ##     imputed_dosages[, i] <- final_set_of_results[[i]]$dosage
        ## }

        for(object in RData_objects_to_save) {
            if (!exists(object)) {
                stop(paste0("You have asked to save object:", object, " as part of RData_objects_to_save but this is not a valid option"))
            }
        }
        save_text <- paste0(
            "save(",
            paste0(RData_objects_to_save, collapse = ", "),
            ", file = output_RData_filename)"
        )
        eval(parse(text = save_text))
        print_message("Done saving extra RData objects to disk")

    }

    print_message("Done QUILT")

    return(NULL)

}
