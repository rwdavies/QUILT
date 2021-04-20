#' @title QUILT_prepare_reference
#' @param outputdir What output directory to use
#' @param chr What chromosome to run. Should match BAM headers
#' @param nGen Number of generations since founding or mixing. Note that the algorithm is relatively robust to this. Use nGen = 4 * Ne / K if unsure
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param output_file Path to output RData file containing prepared haplotypes (has default value that works with QUILT)
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
    minRate = 0.1
) {

    x <- as.list(environment())
    command_line <- paste0(
        "QUILT_prepare_reference(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))

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
    STITCH::validate_reference_files(reference_haplotype_file, reference_legend_file, reference_sample_file, reference_populations, niterations)
    STITCH::validate_regionStart_regionEnd_and_buffer(regionStart, regionEnd, buffer)
    STITCH::validate_tempdir(tempdir)
    STITCH::validate_nGen(nGen)
    ##if (genetic_map_file == "") {
    ##    stop("Please include a genetic map file")
    ##}
    if (genetic_map_file != "") {
        if (!file.exists(genetic_map_file)) {
            stop(paste0("Cannot find file:", genetic_map_file))
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
    if (reference_exclude_samplelist_file != "") {
        if (!file.exists(reference_exclude_samplelist_file)) {
            stop(paste0("Cannot find file:", reference_exclude_samplelist_file))
        }
    }





    
    ##
    ## work on legend and identify SNPs to keep
    ##
    ref_legend <- fread(cmd = paste0("gunzip -c ", shQuote(reference_legend_file)), data.table = FALSE)
    colnames(ref_legend) <- c("CHR", "POS", "REF", "ALT")
    ref_legend[, 1] <- chr
    if (!is.na(regionStart)) {
        a <- unique(ref_legend[
        (ref_legend[, "POS"] >= (regionStart - buffer)) &
        (ref_legend[, "POS"] <= (regionEnd   + buffer)), "POS"])
        pos_to_use <- ref_legend[match(a, ref_legend[, "POS"]), ]
    } else {
        pos_to_use <- ref_legend
    }


    ##
    ## remove sites from consideration
    ## 
    if (region_exclude_file != "") {
        if (!file.exists(region_exclude_file)) {
            stop(paste0("Cannot find region_exclude_file:", region_exclude_file))
        }
        print_message("Loading list of regions to exclude")
        regionstoexclude <- read.table(region_exclude_file, header = TRUE, as.is = TRUE)
        regionstoexclude <- regionstoexclude[regionstoexclude[, 2] == chr,, drop = FALSE]
        if (nrow(regionstoexclude) == 0) {
            warning("No regions to exclude based on region_exclude_file. Perhaps the chr isn't the same?")
        } else {
            regionstoexclude <- as.matrix(regionstoexclude)
            excluderegions <- matrix(as.double(regionstoexclude[, 3:4]),ncol = 2)
            rownames(excluderegions) <- regionstoexclude[, 1]
            ## now do the removal
            keep <- array(TRUE, nrow(pos_to_use))
            oldsum <- 0
            for(i in 1:nrow(excluderegions)){
                to_remove <-
                    pos_to_use[,"POS"] >= excluderegions[i, 1] &
                    pos_to_use[,"POS"] <= excluderegions[i, 2]
                keep[to_remove] <- FALSE
                ##sentence=paste("Excluding ",sum(keep==0)-oldsum," SNPs from region ",rownames(excluderegions)[i])
                ##print(sentence)
                ##oldsum <- sum(keep == 0)
            }
            pos_to_use <- pos_to_use[keep, ]
            print_message(paste0("Excluded ", length(keep) - sum(keep), " out of ", length(keep), " SNPs from ", nrow(excluderegions), " regions"))
        }
    }

    
    ##
    ## (optional) make fake vcf with sites list
    ##
    if (make_fake_vcf_with_sites_list) {
        print_message("Make VCF with sites list")
        if (is.na(output_sites_filename)) {
            output_sites_filename <- file.path(outputdir, paste0("quilt.sites.", regionName, ".vcf.gz")) 
        } else {
            STITCH::validate_output_filename(output_sites_filename, output_format = "bgvcf")
        }
        cat(
            "##fileformat=VCFv4.2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            sep = "\n",
            file = gsub(".gz", "", output_sites_filename)
        )
        to_out <- cbind(
            pos_to_use[ , 1],
            pos_to_use[, 2],
            ".",
            pos_to_use[, 3],
            pos_to_use[, 4],
            1000,
            "PASS",
            'DP=1000',
            "GT",
            "0/0"
        )
        write.table(to_out, file = gsub(".gz", "", output_sites_filename), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
        check_program_dependency("bgzip")
        check_program_dependency("tabix")
        system(paste0("bgzip -f ", gsub(".gz", "", output_sites_filename)))
        system(paste0("tabix -f ", output_sites_filename))
        print_message("Done making VCF with sites list")        
    }


    ##
    ## extract haplotypes from reference
    ##
    out <- get_haplotypes_from_reference(
        reference_haplotype_file = reference_haplotype_file,
        reference_legend_file = reference_legend_file,
        reference_sample_file = reference_sample_file,
        reference_populations = reference_populations,
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
    
    ##outORI <- out ## in case
    ref_alleleCount <- out[["ref_alleleCount3"]] ## defined at all SNPs
    gc(reset =TRUE); gc(reset =TRUE); gc(reset =TRUE)



    
    ##
    ## stuff common across number of samples in ref output
    ##
    pos <- pos_to_use
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
        ## we need this here, build it!
        ## so expRate is cM / Mbp
        cM <- c(0, cumsum(diff(L) * expRate)) / 1e6
    }


    ##
    ##
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
    ## inRegion2 <- (L >= (regionStart - buffer) & L <= (regionEnd + buffer))
    nGrids <- out2$nGrids
    grid <- out2$grid
    L_grid <- out2$L_grid
    dl <- diff(L_grid)
    

    ##
    ## might not be necessary
    ##
    ## nGen <- 100 ## 4 * 20000 / nrow(rhb_t)
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


    ##
    ##
    ##
    rateBetweenGrids <- nGen * diff(cM_grid) / 100
    expRateBetweenGrids <- nGen * diff(L_grid) / 1e6 * (expRate / 100)
    minRateBetweenGrids <- nGen * diff(L_grid) / 1e6 * (minRate / 100)
    m <- cbind(expRateBetweenGrids, rateBetweenGrids, minRateBetweenGrids)
    m[, 1] <- m[, 1] * diff(L_grid)
    m[, 2] <- m[, 2] * diff(L_grid)
    m[, 3] <- m[, 3] * diff(L_grid)
    m <- cbind(diff(L_grid), m)
    ## 
    w <- rateBetweenGrids < minRateBetweenGrids
    ## print(paste0("There are ", sum(w), " regions below minimum rate, setting them to minimum rate"))
    rateBetweenGrids[w] <- minRateBetweenGrids[w]
    sigmaCurrent_m <- array(exp(-rateBetweenGrids), c(nGrids - 1, S))
    ## print(paste0("The probability of staying the same across the region is ", prod(sigmaCurrent_m)))
    
    ##
    ## stuff to do with the panel and possibly removing people
    ##
    ##     out <- outORI
    rhb <- out[["rhb3"]]
    rh_in_L <- out[["rh_in_L"]]
    ref_alleleCount <- out[["ref_alleleCount3"]] ## defined at all SNPs
    ##ref_samples <- reference_samples
    rhb_t <- t(rhb)

    ##
    ## try to grab NA12878 for convenience
    ##
    ## if ("sample" %in% colnames(reference_samples)) {
    ##     na_haps <- which(rep(reference_samples[, "sample"], each = 2) =="NA12878")
    ## } else {
    ##     na_haps <- c()
    ## }
    ## na12878_hap1 <- rcpp_int_expand(rhb_t[na_haps[1], ], nSNPs)
    ## na12878_hap2 <- rcpp_int_expand(rhb_t[na_haps[2], ], nSNPs)    
    ## rhb_t_no_NA12878 <- rhb_t[setdiff(1:nrow(rhb_t), na_haps), ]
    ## rhb_t <- rhb_t_no_NA12878
    ## ref_samples <- ref_samples[ref_samples[, "sample"] != "NA12878", ]


    ##
    ## possibly exclude samples
    ##
    if (reference_sample_file != "") {
        reference_samples <- fread(reference_sample_file, header = TRUE, data.table = FALSE)
        if (!is.na(reference_populations[1])) {
            keep_samples <- as.character(reference_samples[, "POP"]) %in% reference_populations
            reference_samples <- reference_samples[keep_samples, ]
        }
        if (nrow(rhb_t) != (2 * nrow(reference_samples))) {
            stop(paste0("The number of haplotypes from the reference haplotype file (N = ", nrow(rhb_t), ") does not match the inferred number of entries from the reference samples file (Nrows = ", nrow(reference_samples), ", N = ", 2 * nrow(reference_samples)))
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


    ##
    ## do compression here
    ## note now can work out recommended nMaxDH on the fly
    ## 
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


    ##
    ## save here!
    ##
    ##  na12878_hap1
    ##  na12878_hap2
    print_message("Save converted reference haplotypes")    
    save(
        ref_error,
        hapMatcher,
        distinctHapsIE,
        distinctHapsB,
        eMatDH_special_grid_which,
        eMatDH_special_values_list,
        inRegion2,
        rhb_t,
        reference_samples,
        rh_in_L,
        af,
        ref_alleleCount,    
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
        file = output_file,
        compress = FALSE
    )
    

    print_message("Done converting reference haplotypes")

}

