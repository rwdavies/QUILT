#' @title QUILT_prepare_reference
#' @param outputdir What output directory to use
#' @param chr What chromosome to run. Should match BAM headers
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param output_file Path to output RData file containing prepared haplotypes (has default value that works with QUILT)
#' @param reference_vcf_file Path to reference haplotype file in VCF format (values must be 0 or 1)
#' @param reference_haplotype_file Path to reference haplotype file in IMPUTE format (file with no header and no rownames, one row per SNP, one column per reference haplotype, space separated, values must be 0 or 1)
#' @param reference_legend_file Path to reference haplotype legend file in IMPUTE format (file with one row per SNP, and a header including position for the physical position in 1 based coordinates, a0 for the reference allele, and a1 for the alternate allele)
#' @param reference_sample_file Path to reference sample file (file with header, one must be POP, corresponding to populations that can be specified using reference_populations)
#' @param reference_populations Vector with character populations to include from reference_sample_file e.g. CHB, CHS
#' @param reference_phred Phred scaled likelihood or an error of reference haplotype. Higher means more confidence in reference haplotype genotypes, lower means less confidence
#' @param reference_exclude_samplelist_file File with one column of samples to exclude from reference samples e.g. in validation, the samples you are imputing
#' @param region_exclude_file File with regions to exclude from constructing the reference panel. Particularly useful for QUILT_HLA, where you want to exclude SNPs in the HLA genes themselves, so that reads contribute either to the read mapping or state inference. This file is space separated with a header of Name, Chr, Start and End, with Name being the HLA gene name (e.g. HLA-A), Chr being the chromosome (e.g. chr6), and Start and End are the 1-based starts and ends of the genes (i.e. where we don't want to consider SNPs for the Gibbs sampling state inference)
#' @param genetic_map_file Path to file with genetic map information, a file with 3 white-space delimited entries giving position (1-based), genetic rate map in cM/Mbp, and genetic map in cM. If no file included, rate is based on physical distance and expected rate (expRate)
#' @param nMaxDH Integer Maximum number of distinct haplotypes to store in reduced form. Recommended to keep as 2 ** N - 1 where N is an integer greater than 0 i.e. 255, 511, etc
#' @param tempdir What directory to use as temporary directory. If set to NA, use default R tempdir. If possible, use ramdisk, like /dev/shm/
#' @param make_fake_vcf_with_sites_list Whether to output a list of sites as a minimal VCF, for example to use with GATK 3 to genotype given sites
#' @param output_sites_filename If make_fake_vcf_with_sites_list is TRUE, optional desired filename where to output sites VCF
#' @param expRate Expected recombination rate in cM/Mb
#' @param maxRate Maximum recomb rate cM/Mb
#' @param minRate Minimum recomb rate cM/Mb
#' @param use_mspbwt Build mspbwt indices to be used in imputation
#' @param mspbwt_nindices How many mspbwt indices to build
#' @param override_use_eMatDH_special_symbols Not for general use. If NA will choose version appropriately depending on whether a PBWT flavour is used.
#' @param use_hapMatcherR Used for nMaxDH less than or equal to 255. Use R raw format to hold hapMatcherR. Lowers RAM use
#' @param impute_rare_common Whether to use common SNPs first for imputation, followed by a round of rare imputation
#' @param rare_af_threshold Allele frequency threshold under which SNPs are considered rare, otherwise they are considered common.
#' @param use_list_of_columns_of_A If when using mspbwt, use columns of A rather than the whole thing, to speed up this version
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
QUILT_prepare_reference <- function(
    outputdir,
    chr,
    nGen,
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    output_file = "",
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
    tempdir = NA,
    make_fake_vcf_with_sites_list = FALSE,
    output_sites_filename = NA,
    expRate = 1,
    maxRate = 100,
    minRate = 0.1,
    use_mspbwt = FALSE,
    mspbwt_nindices = 4L,
    override_use_eMatDH_special_symbols = NA,
    use_hapMatcherR = TRUE,
    impute_rare_common = FALSE,
    rare_af_threshold = 0.001,
    use_list_of_columns_of_A = TRUE
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_prepare_reference(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))


    ## mspbwtMAF = 0.0001,
    ##         @param mspbwtMAF Building mspbwt indices for common and rare variants seperately
    
    
    print_message("Program start")
    print_message("Begin converting reference haplotypes")

    ## need for outputting
    regionName <- chr
    options(scipen = 999) ## Dangit rounding
    if(is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE)
        regionName <- paste0(chr, ".", regionStart,".", regionEnd)


    ##
    ## validate parameters
    ##
    STITCH::validate_chr(chr)
    niterations <- NA ## not needed
    ## STITCH::validate_reference_files(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations)
    STITCH::validate_regionStart_regionEnd_and_buffer(regionStart, regionEnd, buffer)
    STITCH::validate_tempdir(tempdir)
    STITCH::validate_nGen(nGen)
    if (genetic_map_file != "") {
        if (!file.exists(genetic_map_file)) {
            stop(paste0("Cannot find file:", genetic_map_file))
        }
    }
    if (!(class(mspbwt_nindices) %in% c("integer", "numeric"))) {
        stop("mspbwt_nindices must be either numeric or integer")
    }
    if (mspbwt_nindices < 1) {
        stop("mspbwt_nindices must be at least 1")
    }
    if (round(mspbwt_nindices) != mspbwt_nindices) {
        stop("mspbwt_nindices must be an integer")
    }
    if (!is.na(nMaxDH)) {
        if (nMaxDH > 255 & use_hapMatcherR) {
            stop("Can only use hapMatcherR with nMaxDH <= 255")
        }
    }


    ##
    ## new validations
    ##
    if (!is.na(nMaxDH)) {
        validate_nMaxDH(nMaxDH)
    }
    ref_error <- 10 ** (-reference_phred / 10)
    dir.create(outputdir, showWarnings = FALSE)
    dir.create(file.path(outputdir, "RData"), showWarnings = FALSE)
    if (output_file == "") {
        output_file <- file_quilt_prepared_reference(outputdir, regionName)
    } else {
        if (!dir.exists(basename(output_file))) {
            dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
        }
        if (basename(output_file) == output_file) {
            output_file <- file.path(outputdir, output_file)
        }
    }
    if (is.na(tempdir)) {
        tempdir <- tempdir()
    }


    ##
    ## validate reference file options
    ##
    if (reference_vcf_file != "") {
        print_message(paste0("Using reference information from:", reference_vcf_file))
        use_reference_vcf <- TRUE
        if (region_exclude_file != "") {
            stop("You have given reference_vcf_file, however this does not yet work when using a reference VCF")
        }
        if (make_fake_vcf_with_sites_list) {
            stop("You have given make_fake_vcf_with_sites_list, however this does not yet work when using a reference VCF")
        }
        ##if (reference_sample_file !== "") {
        ##    stop("You have given reference_sample_file, however this does not yet work when using a reference VCF")
        ## }
        ##if (!is.na(reference_populations)) {
        ##    stop("You have given reference_populations, however this does not yet work when using a reference VCF")
        ##}
    } else if (reference_haplotype_file != "" && reference_legend_file != "") {
        if (use_mspbwt) {
            stop("Please supply reference_vcf_file when use_mspbwt")
        }
        print_message(paste0("Using reference information from:", reference_haplotype_file, " and ", reference_legend_file))
        use_reference_vcf <- FALSE
        if (reference_exclude_samplelist_file != "") {
            if (!file.exists(reference_exclude_samplelist_file)) {
                stop(paste0("Cannot find file:", reference_exclude_samplelist_file))
            }
        }
        if (reference_sample_file != "" & file.exists(reference_sample_file) == FALSE) 
            stop(paste0("Cannot find reference_sample_file:", reference_sample_file))
        if (reference_legend_file != "" & file.exists(reference_legend_file) == FALSE) 
            stop(paste0("Cannot find reference_legend_file:", reference_legend_file))
        if (reference_haplotype_file != "" & file.exists(reference_haplotype_file) == FALSE) 
            stop(paste0("Cannot find reference_haplotype_file:", reference_haplotype_file))
    } else {
        stop("Please include either reference_vcf_file, or reference_legend_file and reference_haplotype_file, as appropriate")
    }


    if (use_reference_vcf) {

        if (!impute_rare_common) {
            ## disable
            rare_af_threshold <- 0
        } else {
            print_message(paste0("Using strategy of first imputing using common SNPs and then using all SNPs, with allele frequency threshold:", rare_af_threshold))
        }
        
        ifelse(regionStart-buffer<1, samtoolslike <- paste0(chr, ":", 1, "-", regionEnd+buffer), samtoolslike <- paste0(chr, ":", regionStart-buffer, "-", regionEnd+buffer) )
        print_message("Begin get sites and haplotypes from reference vcf")


        if (reference_sample_file == "") {
            reference_samples <- NULL
            subsamples <- "-" ## all samples            
            if (!is.na(reference_populations[1])) {
                stop("You have selected reference populations to include, but have not included reference_sample_file with reference sample information")
            }
            if (reference_exclude_samplelist_file != "") {
              s <- read.table(reference_exclude_samplelist_file, h = FALSE)[,1]
              subsamples <- paste0("^", paste0(s, collapse = ",")) ## ^ means excluding in vcfpp.h
            }
        } else {
            reference_samples <- data.table::fread(reference_sample_file, header = TRUE, data.table = FALSE)
            if (reference_exclude_samplelist_file != "") {
                s <- read.table(reference_exclude_samplelist_file, h = FALSE)[,1]
                w <- (!reference_samples[, 1] %in% s)
                if (sum(w) == 0) {
                    stop("All of the reference samples would be excluded with this choice of reference_exclude_samplelist_file")
                }
                reference_samples <- reference_samples[w, ]
            }
            if (!is.na(reference_populations[1])) {
                w <- reference_samples[, 2] %in% reference_populations
                if (sum(w) == 0) {
                    stop("None of the retained reference samples have populations (column 2 in the reference_sample_file) that match the supplied reference_populations")
                }
                reference_samples <- reference_samples[w, ]
            }
            subsamples <- paste0(reference_samples[, 1], collapse = ",")
        }

      out <- STITCH::Rcpp_get_hap_info_from_vcf(
        vcffile = reference_vcf_file,
        af_cutoff = rare_af_threshold,
        region = samtoolslike,
        samples = subsamples
      )
        print_message("End get sites and haplotypes from reference vcf")

        ref_alleleCount <- out[["ref_alleleCount"]]
        rhb_t <- out[["rhb_t"]]
        pos <- cbind(chr, out[["pos"]]) ## make normal shape


        rare_per_hap_info <- out[["rare_per_hap_info"]]
        snp_is_common <- out[["snp_is_common"]]
        n_skipped <- out[["n_skipped"]]
        print_message(paste0("There were ", n_skipped, " skipped variants when processing the reference VCF (not bi-allelic, not a SNP, not unique position"))

        rm(out)
        gc(reset = TRUE); gc(reset = TRUE); gc(reset = TRUE);

        if (!impute_rare_common) {
            STITCH::validate_region_to_impute_when_using_regionStart(pos[, 2], regionStart, regionEnd, buffer)
        } else {
            rare_common_validate_region_to_impute_when_using_regionStart(pos[, 2], regionStart, regionEnd, buffer, snp_is_common)
        }
        
        ## re-define pos here, as being for common SNPs
        ## then keep around the all bit for later
        pos_all <- pos        
        pos <- pos[snp_is_common, ]
        ref_alleleCount_all <- ref_alleleCount
        ref_alleleCount <- ref_alleleCount[snp_is_common, ]
        

    } else {

        pos_to_use <- get_pos_to_use_from_reference_legend_file(reference_legend_file, chr, regionStart, regionEnd, buffer)

        ## potentially remove sites from consideration
        if (region_exclude_file != "") {
            pos_to_use <- remove_sites_from_pos_to_use(
                region_exclude_file = region_exclude_file,
                pos_to_use = pos_to_use,
                chr = chr
            )
        }

        ## (optional) make fake vcf with sites list
        if (make_fake_vcf_with_sites_list) {
            make_face_vcf_with_sites_list(outputdir, regionName, output_sites_filename, pos_to_use) 
        }

        ## extract haplotypes from reference
        out <- get_haplotypes_from_reference(
            reference_haplotype_file = reference_haplotype_file,
            reference_legend_file = reference_legend_file,
            reference_sample_file = reference_sample_file,
            reference_populations = reference_populations, ## this is done below, to synchronize
            pos = pos_to_use,
            tempdir = tempdir,
            regionName = regionName,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            chr = chr,
            niterations = 2,
            extraction_method = "hap_v3" ## to do - make both, then only the one I want
        )
        
        ref_alleleCount <- out[["ref_alleleCount3"]] ## defined at all SNPs
        rhb <- out[["rhb3"]]
        ## rh_in_L <- out[["rh_in_L"]]
        rhb_t <- t(rhb)
        pos <- pos_to_use
        rm(pos_to_use)
        rm(rhb)
        rm(out)
        gc(reset = TRUE); gc(reset = TRUE); gc(reset = TRUE);

        ## needed for rare common idea, only from VCF approach
        rare_per_hap_info <- NULL
        snp_is_common <- rep(TRUE, nrow(pos))
        n_skipped <- NULL
        pos_all <- pos
        ref_alleleCount_all <- ref_alleleCount


        ##
        ## possibly exclude samples, note, does not fix ref_alleleCount!
        ##
        if (reference_sample_file != "") {
            reference_samples <- data.table::fread(reference_sample_file, header = TRUE, data.table = FALSE)
            if (!use_reference_vcf) {
                if (!is.na(reference_populations[1])) {
                    keep_samples <- as.character(reference_samples[, "POP"]) %in% reference_populations
                    reference_samples <- reference_samples[keep_samples, ]
                }
                if (nrow(rhb_t) != (2 * nrow(reference_samples))) {
                    stop(paste0("The number of haplotypes from the reference haplotype file (N = ", nrow(rhb_t), ") does not match the inferred number of entries from the reference samples file (Nrows = ", nrow(reference_samples), ", N = ", 2 * nrow(reference_samples), ")"))
                }
            }
        } else {
            reference_samples <- NULL
        }
        if (reference_exclude_samplelist_file != "") {
            ## validated above
            if (is.null(reference_samples)) {
                stop(paste0("You have requested to exclude reference samples with reference_exclude_samplelist_file but you have not supplied reference_sample_file"))
            }
            exclude_samples <- read.table(reference_exclude_samplelist_file, stringsAsFactors = FALSE)[, 1]
            t1 <- which(rep(reference_samples[, 1], each = 2) %in% exclude_samples)
            rhb_t <- rhb_t[-t1, ]
            reference_samples <- reference_samples[-which(reference_samples[, 1] %in% exclude_samples), ]
        }
        
    }
        
    

    ##
    ## stuff common across number of samples in ref output
    ##
    L <- pos[, 2]
    nSNPs <- nrow(pos)
    af <- ref_alleleCount[, 3]
    ancAlleleFreqAll <- af
    ## returns NULL if no file
    cM <- load_validate_and_match_genetic_map(
        genetic_map_file = genetic_map_file,
        L = L,
        expRate = expRate
    )
    if (is.null(cM)) {
        cM <- c(0, cumsum(diff(L) * expRate)) / 1e6
    }


    ##
    ## grids
    ##
    if (is.na(regionStart)) {
        inRegion2 <- array(TRUE, length(L))
    } else {
        inRegion2 <- (regionStart <= L) & (L <= regionEnd)
    }
    out2 <- assign_positions_to_grid(
        L = pos[, 2],
        grid32 = TRUE,
        cM = cM
    )
    nGrids <- out2$nGrids
    grid <- out2$grid
    L_grid <- out2$L_grid
    dl <- diff(L_grid)


    ##
    ## might not be necessary
    ##
    S <- 1
    if (genetic_map_file != "") {
        genetic_map <- get_and_validate_genetic_map(genetic_map_file)
    } else {
        ## hack it
        genetic_map <- cbind(L, expRate, cM)
        colnames(genetic_map) <- c(
            "position", "COMBINED_rate.cM.Mb.", "Genetic_Map.cM."
        )
    }
    cM_grid <- match_genetic_map_to_L(
        genetic_map = genetic_map,
        L = L_grid,
        expRate = expRate
    )
    sigmaCurrent_m <- get_sigmaCurrent_m(nGen, cM_grid, L_grid, expRate, minRate, maxRate, nGrids, S)

    




    ##
    ## do compression here
    ## note now can work out recommended nMaxDH on the fly
    ##
    out <- STITCH::make_rhb_t_equality(
        rhb_t = rhb_t,
        nMaxDH = nMaxDH,
        nSNPs = nSNPs,
        ref_error = ref_error,
        use_hapMatcherR = use_hapMatcherR
    )
    distinctHapsB <- out[["distinctHapsB"]]
    distinctHapsIE <- out[["distinctHapsIE"]]
    hapMatcher <- out[["hapMatcher"]]
    hapMatcherR <- out[["hapMatcherR"]]
    eMatDH_special_grid_which <- out[["eMatDH_special_grid_which"]]
    eMatDH_special_values_list <- out[["eMatDH_special_values_list"]]

    if (!is.na(override_use_eMatDH_special_symbols)) {
        use_eMatDH_special_symbols <- override_use_eMatDH_special_symbols
    } else {
        if (use_mspbwt | impute_rare_common) {
            use_eMatDH_special_symbols <- TRUE
        } else {
            use_eMatDH_special_symbols <- FALSE
        }
    }

    if (use_eMatDH_special_symbols) {
        eMatDH_special_matrix_helper <- out[["eMatDH_special_matrix_helper"]]
        eMatDH_special_matrix <- out[["eMatDH_special_matrix"]]
        if (use_mspbwt) {
            rhb_t <- matrix(as.integer(1), 1, 1) ## nuke!
            gc(reset = TRUE); gc(reset = TRUE); 
        }
    } else {
        eMatDH_special_matrix_helper <- matrix(as.integer(1), 1, 1) ## nuke!
        eMatDH_special_matrix <- matrix(as.integer(1), 1, 1) ## nuke!
    }

    if (nGrids < mspbwt_nindices) {
        print_message(paste0("There are ", nGrids, " grids, so re-setting mspbwt_nindices to ", 1))
        mspbwt_nindices <- 1L
    }

    if (use_mspbwt) {
        print_message("Build mspbwt indices")
        all_symbols <- out$all_symbols
        ## if (temp_pre_save) {
        ##     print("temp pre save to /data/smew1/rdavies/ukbb_gel_2023_01_26/2023_06_22/quilt_mspbwt_FALSE_TRUE_488315_TRUE/RData/temp_presave.RData")
        ##     save(hapMatcher,hapMatcherR,mspbwt_nindices,use_hapMatcherR,all_symbols, file = "/data/smew1/rdavies/ukbb_gel_2023_01_26/2023_06_22/quilt_mspbwt_FALSE_TRUE_488315_TRUE/RData/temp_presave.RData")
        ##     print("end of temp pre save to /data/smew1/rdavies/ukbb_gel_2023_01_26/2023_06_22/quilt_mspbwt_FALSE_TRUE_488315_TRUE/RData/temp_presave.RData")
        ## }
        ms_indices <- build_mspbwt_indices(
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            mspbwt_nindices = mspbwt_nindices,
            use_hapMatcherR = use_hapMatcherR,
            all_symbols = all_symbols,
            use_list_of_columns_of_A = use_list_of_columns_of_A
        )
        print_message("Done building mspbwt indices")
    } else {
        ms_indices <- NULL
    }


    ##
    ## save here!
    ##
    print_message("Save converted reference haplotypes")
    ## rh_in_L,
    save(
        ref_error,
        hapMatcher,
        hapMatcherR,
        use_hapMatcherR,
        distinctHapsIE,
        distinctHapsB,
        eMatDH_special_grid_which,
        eMatDH_special_values_list,
        eMatDH_special_matrix,
        eMatDH_special_matrix_helper,
        inRegion2,
        rhb_t,
        reference_samples,
        af,
        ref_alleleCount,
        ref_alleleCount_all,
        pos_all,
        pos,
        L,
        nSNPs,
        nGrids,
        grid,
        L_grid,
        dl,
        cM,
        cM_grid,
        sigmaCurrent_m,
        regionStart,
        regionEnd,
        buffer,
        chr,
        ms_indices,
        nGen,
        use_eMatDH_special_symbols,
        rare_per_hap_info,
        snp_is_common,
        n_skipped,
        genetic_map,
        file = output_file,
        compress = FALSE
    )


    print_message("Done converting reference haplotypes")

}

