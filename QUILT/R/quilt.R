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
#' @param Ksubset How many haplotypes to use in the faster Gibbs sampling
#' @param Knew How many haplotypes to replace per-iteration after doing the full reference panel imputation
#' @param K_top_matches How many top haplotypes to store in each grid site when looking for good matches in the full haplotype reference panel. Large values potentially bring in more haplotype diversity, but risk losing haplotypes that are good matches over shorter distances
#' @param heuristic_match_thin What fraction of grid sites to use when looking for good matches in the full haplotype reference panel. Smaller values run faster but potentially miss haplotypes
#' @param output_filename Override the default bgzip-VCF / bgen output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple QUILT runs are processing on the same region then they may over-write each others inputs and outputs.
#' @param RData_objects_to_save Can be used to name interim and misc results from imputation to save an an RData file. Default NULL means do not save such output
#' @param output_RData_filename Override the default location for miscellaneous outputs saved in RData format
#' @param tempdir What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/
#' @param bqFilter Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used
#' @param panel_size Integer number of reference haplotypes to use, set to NA to use all of them
#'
#' @param posfile Optional, only needed when using genfile. File with positions of where to impute, lining up one-to-one with genfile. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>
#' @param genfile Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header
#' @param phasefile Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header. For NIPT imputation, there are 3 columns, representing maternal transmitted, maternal untransmitted, and paternal transmitted
#' @param maxDifferenceBetweenReads How much of a difference to allow the reads to make in the forward backward probability calculation. For example, if P(read | state 1)=1 and P(read | state 2)=1e-6, re-scale so that their ratio is this value. This helps prevent any individual read as having too much of an influence on state changes, helping prevent against influence by false positive SNPs
#' @param make_plots Whether to make some plots of per-sample imputation. Especially nice when truth data. This is pretty slow though so useful more for debugging and understanding and visualizing performance
#' @param verbose whether to be more verbose when running
#' @param shuffle_bin_radius Parameter that controls how to detect ancestral haplotypes that are shuffled during EM for possible re-setting. If set (not NULL), then recombination rate is calculated around pairs of SNPs in window of twice this value, and those that exceed what should be the maximum (defined by nGen and maxRate) are checked for whether they are shuffled
#' @param iSizeUpperLimit Do not use reads with an insert size of more than this value
#' @param record_read_label_usage Whether to store what read labels were used during the Gibbs samplings (i.e. whether reads were assigned to arbitrary labelled haplotype 1 or 2)
#' @param record_interim_dosages Whether to record interim dosages or not
#' @param use_bx_tag Whether to try and use BX tag in same to indicate that reads come from the same underlying molecule
#' @param bxTagUpperLimit When using BX tag, at what distance between reads to consider reads with the same BX tag to come from different molecules
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
    Ksubset = 400,
    Knew = 100,
    K_top_matches = 5,
    heuristic_match_thin = 0.1,
    output_filename = NULL,
    RData_objects_to_save = NULL,
    output_RData_filename = NULL,
    prepared_reference_filename = "",
    tempdir = NA,
    bqFilter = as.integer(17),
    panel_size = NA,
    posfile = "",
    genfile = "",
    phasefile = "",
    maxDifferenceBetweenReads = 1e10,
    make_plots = FALSE,
    verbose = TRUE,
    shuffle_bin_radius = 5000,
    iSizeUpperLimit = 1e6,
    record_read_label_usage = TRUE,
    record_interim_dosages = TRUE,
    use_bx_tag = TRUE,
    bxTagUpperLimit = 50000
) {

    ## init_method <- "simple"
    ## use_eMatDH <- TRUE
    ## maxDifferenceBetweenReads <- 1e10 ## OK hopefully?
    ## full_expRate <- NA ## not being used currently
    ##full_nGen <- NA ## not being used currently
    ## full_error <- 1e-3
    ## verbose <- TRUE
    ##tempdir <- tempdir()
    ## shuffle_bin_radius <- 5000
    ##Ksubset <- 400 ## how many good haps to start with
    ##Knew <- 100 ## when updating how many new haplotypes to get from each round
    ##K_top_matches <- 5 ## how many top haplotypes to priotize
    ##heuristic_match_thin <- 0.1 ## how much thinning for faster checking
    ##truth_haps_all <- NULL
    ##record_read_label_usage <- FALSE
    ##nGibbsSamples <- as.integer(args[16])
    ##n_seek_its <- as.integer(args[17])

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))
    
    options(digits.secs=6)
    options(scipen = 999)

    ## need for outputting
    regionName <- chr
    options(scipen = 999) ## Dangit rounding
    if(is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE)
        regionName <- paste0(chr, ".", regionStart,".", regionEnd)

    if (make_plots) {
        dir.create(file.path(outputdir, "plots"), showWarnings = FALSE)
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
    if (is.na(tempdir)) {
        tempdir <- tempdir()
    }
    if (!is.null(output_filename)) {
        if (!dir.exists(basename(output_filename))) {
            dir.create(basename(output_filename), showWarnings = FALSE)
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
    if (prepared_reference_filename == "") {
        prepared_reference_filename <- file_quilt_prepared_reference(outputdir, regionName) 
    }
    if (!file.exists(prepared_reference_filename)) {
        stop(paste0("Cannot find prepared haplotype reference file, expecting:", prepared_reference_filename))
    }
    load(prepared_reference_filename)


    ##
    ## optionally load genotypes and phasevali
    ##
    if (genfile != "") {
        if (posfile == "") {
            stop("You have given genfile, a set of high confidence genotypes, but not posfile, which confirms their physical positions. Please fix this, and see the help files for more information")
        }
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
    err <- paste0("There was an error lining up the reference data with the prepared reference haplotype file. Please double check that posfile and genfile are defined as exactly the same SNPs in the region to be imputed as the reference legend file")
    if (length(key1) != length(key2)) {
        stop(err)
    }
    if (sum(key1 != key2) > 0) {
        stop(err)
    }
    

    
    ##
    ## can request smaller panel size, only really useful for testing speed
    ##
    if (!is.na(panel_size)) {
        rhb_t <- rhb_t[1:as.integer(panel_size), ] ## this is the number of HAPLOTYPES
        reference_samples <- reference_samples[1:as.integer(panel_size), ]
        nMaxDH <- 2 ** 8 - 1
        out <- make_rhb_t_equality(
            rhb_t = rhb_t,
            nMaxDH = nMaxDH,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
        distinctHapsB <- out[["distinctHapsB"]]
        distinctHapsIE <- out[["distinctHapsIE"]]
        hapMatcher <- out[["hapMatcher"]]
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
        alphaHat_t2 <- array(0, c(K, nGrids))
        betaHat_t2 <- array(0, c(K, nGrids))
        alphaHat_t3 <- array(0, c(K, nGrids))
        betaHat_t3 <- array(0, c(K, nGrids))
        small_priorCurrent_m <- array(1 / K, c(K, S))
        small_alphaMatCurrent_tc <- array(1 / K, c(K, nGrids - 1, S))
        small_eHapsCurrent_tc <- array(0, c(K, nSNPs, S))

        results_across_samples <- as.list(sampleRange[2] - sampleRange[1] + 1)
        
        for(iSample in sampleRange[1]:sampleRange[2]) {

            print_message(paste0("Imputing sample: ", iSample))
            out <- get_and_impute_one_sample(
                rhb_t = rhb_t,
                outputdir = outputdir,
                nGibbsSamples = nGibbsSamples,
                n_seek_its = n_seek_its,
                full_alphaHat_t = full_alphaHat_t,
                full_betaHat_t = full_betaHat_t,
                full_gamma_t = full_gamma_t,
                full_gammaSmall_t = full_gammaSmall_t,
                full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                full_transMatRate_t_H = full_transMatRate_t_H,
                small_transMatRate_tc_H = small_transMatRate_tc_H,
                alphaHat_t1 = alphaHat_t1,
                betaHat_t1 = betaHat_t1,
                alphaHat_t2 = alphaHat_t2,
                betaHat_t2 = betaHat_t2,
                alphaHat_t3 = alphaHat_t3,
                betaHat_t3 = betaHat_t3,
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
                scenario = scenario,
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
                truth_gen = truth_gen,
                smooth_cm = smooth_cm,
                iSizeUpperLimit = iSizeUpperLimit,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                make_plots = make_plots,
                ref_error = ref_error,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                hapMatcher = hapMatcher,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                af = af,
                use_bx_tag = use_bx_tag,
                bxTagUpperLimit = bxTagUpperLimit
            )

            results_across_samples[[iSample - sampleRange[1] + 1]] <- out

            ## for summarization
            infoCount[, 1] <- infoCount[, 1, drop = FALSE] + out[["eij"]]
            infoCount[, 2] <- infoCount[, 2, drop = FALSE] + (out[["fij"]] - out[["eij"]]**2)
            afCount <- afCount + (out[["eij"]]) / 2
            hweCount[out[["max_gen"]]] <- hweCount[out[["max_gen"]]] + 1 ## hmmmmm not ideal
            alleleCount <- alleleCount + out[["per_sample_alleleCount"]]
            
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
        sampleRanges = sampleRanges
    )    

    
    ##
    ## build a singular set of results
    ##
    ## write out VCF now!

    
    if (!is.null(RData_objects_to_save)) {

        print_message("Begin saving extra RData objects to disk")
        
        final_set_of_results <- as.list(1:N)
        c <- 1
        for(i in 1:length(complete_set_of_results)) {
            x <- complete_set_of_results[[i]]
            for(j in 1:length(x)) {
                final_set_of_results[[c]] <- x[[j]]
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
        print(paste0("Saving using the following:", save_text))
        eval(parse(text = save_text))
        print_message("Done saving extra RData objects to disk")
        
    }
    
    print_message("Done QUILT")
    
    return(NULL)



}
