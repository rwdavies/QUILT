get_and_impute_one_sample <- function(
    rhb_t,
    outputdir,
    nGibbsSamples,
    n_seek_its,
    n_burn_in_seek_its,
    full_alphaHat_t,
    full_betaHat_t,
    full_gamma_t,
    full_gammaSmall_t,
    full_gammaSmall_cols_to_get,
    full_transMatRate_t_H,
    small_transMatRate_tc_H,
    alphaHat_t1,
    betaHat_t1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    eMatGrid_t3,
    gammaMT_t_local,
    gammaMU_t_local,
    gammaP_t_local,
    small_alphaMatCurrent_tc,
    small_priorCurrent_m,
    small_eHapsCurrent_tc,
    bam_files,
    cram_files,
    L,
    pos,
    chr,
    tempdir,
    regionName,
    regionStart,
    regionEnd,
    buffer,
    gen,
    phase,
    gen_all,
    phase_all,
    iSample,
    grid,
    ancAlleleFreqAll,
    L_grid,
    verbose,
    shuffle_bin_radius,
    Ksubset,
    Knew,
    K_top_matches,
    heuristic_match_thin,
    record_interim_dosages,
    have_truth_haplotypes,
    have_truth_genotypes, 
    bqFilter,
    record_read_label_usage,
    sampleNames,
    smooth_cm,
    iSizeUpperLimit,
    maxDifferenceBetweenReads,
    make_plots,
    ref_error,
    distinctHapsB,
    distinctHapsIE,
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    use_eMatDH_special_symbols,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
    eMatDH_special_grid_which,
    eMatDH_special_values_list,
    inRegion2,
    cM_grid,
    af,
    use_bx_tag,
    bxTagUpperLimit,
    addOptimalHapsToVCF,
    make_plots_block_gibbs,
    estimate_bq_using_truth_read_labels,
    chrStart,
    chrEnd,
    gamma_physically_closest_to,
    hla_run,
    downsampleToCov,
    minGLValue,
    minimum_number_of_sample_reads,
    print_extra_timing_information,
    small_ref_panel_gibbs_iterations,
    small_ref_panel_block_gibbs_iterations,
    plot_per_sample_likelihoods,
    use_small_eHapsCurrent_tc,
    output_gt_phased_genotypes,
    ff_values,
    method,
    mspbwtL,
    mspbwtM,
    msp,
    use_mspbwt,
    ms_indices,
    use_splitreadgl,
    use_sample_is_diploid,
    plot_p1,
    small_ref_panel_skip_equally_likely_reads,
    small_ref_panel_equally_likely_reads_update_iterations,
    shard_check_every_pair,
    use_eigen,
    pos_all,
    special_rare_common_objects,
    special_rare_common_objects_per_core,    
    impute_rare_common,
    make_heuristic_plot,
    heuristic_approach,
    calculate_gamma_on_the_fly
) { ## end of function def

    sample_name <- sampleNames[iSample]
    nSNPs <- nrow(pos)
    nGrids <- ncol(distinctHapsB)
    K <- nrow(hapMatcher)
    suppressOutput <- !print_extra_timing_information
    ff <- ff_values[iSample]
    
    if (impute_rare_common) {

        loadBamAndConvert(
            iBam = iSample,
            L = pos_all[, 2],
            pos = pos_all,
            nSNPs = nrow(pos_all),
            bam_files = bam_files,
            iSizeUpperLimit = iSizeUpperLimit,
            bqFilter = bqFilter,
            chr = chr,
            N = length(sampleNames),
            downsampleToCov = downsampleToCov,
            sampleNames = sampleNames,
            inputdir = tempdir,
            regionName = regionName,
            tempdir = tempdir,
            chrStart = chrStart,
            chrEnd = chrEnd,
            chrLength = NA,
            save_sampleReadsInfo = TRUE,
            use_bx_tag = use_bx_tag,
            bxTagUpperLimit = bxTagUpperLimit,
            default_sample_no_read_behaviour = "return_null"
        )

        load(file_sampleReads(tempdir, iSample, regionName))
        load(file_sampleReadsInfo(tempdir, iSample, regionName))
        removeTmpSamplesFile(tempdir, iSample, regionName, save_sampleReadsInfo = TRUE)
        allSNP_sampleReads <- sampleReads
        allSNP_sampleReadsInfo <- sampleReadsInfo
        rm(sampleReads, sampleReadsInfo)

        allSNP_sampleReads <- snap_sampleReads_to_grid(
            sampleReads = allSNP_sampleReads,
            grid = special_rare_common_objects[["grid"]]
        )

        nGrids_all <- length(special_rare_common_objects$L_grid)
        allSNP_wif0 <- as.integer(sapply(allSNP_sampleReads, function(x) x[[2]]))
        allSNP_grid_has_read <- rep(FALSE, nGrids_all)
        allSNP_grid_has_read[unique(allSNP_wif0) + 1] <- TRUE

        if (have_truth_haplotypes) {
            s <- sampleNames[iSample]
            if (!(s %in% dimnames(phase_all)[[2]])) {
                stop("Something went wrong with phase naming")
            }
            if (method == "diploid") {
                truth_haps_all <- cbind(phase_all[, s, 1], phase_all[, s, 2])
            } else {
                truth_haps_all <- cbind(phase_all[, s, 1], phase_all[, s, 2], phase_all[, s, 3])
            }
        } else {
            truth_haps_all <- NULL
        }

        if (have_truth_genotypes) {
            s <- sampleNames[iSample]
            if (!(s %in% colnames(gen))) {
                stop("Something went wrong with gen naming")
            }
            truth_gen_all <- gen_all[, s, drop = FALSE]
        } else {
            truth_gen_all <- NULL
        }

        if (!is.null(truth_haps_all)) {

            if (method == "diploid") {
                truth_label_set <- determine_a_set_of_truth_labels(
                    sampleReads = allSNP_sampleReads,
                    truth_hap1 = truth_haps_all[, 1],
                    truth_hap2 = truth_haps_all[, 2],
                    maxDifferenceBetweenReads = maxDifferenceBetweenReads
		)
            } else {
                truth_label_set <- determine_a_set_of_truth_labels_for_nipt(
                    sampleReads = allSNP_sampleReads,
                    truth_haps = truth_haps_all,
                    phase = phase,
                    ff = ff
                )
            }
            truth_labels_all <- truth_label_set[["truth_labels"]]
            uncertain_truth_labels_all <- truth_label_set[["uncertain_truth_labels"]]
            rm(truth_label_set)
        } else {
            truth_labels_all <- NULL
            uncertain_truth_labels_all <- NULL
        }

        nSNPs_all <- nrow(pos_all)

        if (method == "diploid") {
            dosage_all <- numeric(nSNPs_all)
            gp_t_all <- array(0, c(3, nSNPs_all))
        } else {
            mat_dosage_all <- numeric(nSNPs_all)
            fet_dosage_all <- numeric(nSNPs_all)            
            mat_gp_t_all <- array(0, c(3, nSNPs_all))
            fet_gp_t_all <- array(0, c(3, nSNPs_all))            
        }


    }

    ##
    ## sample read stuff - work off bam file!
    ##
    loadBamAndConvert(
        iBam = iSample,
        L = L,
        pos = pos,
        nSNPs = nSNPs,
        bam_files = bam_files,
        cram_files = cram_files,
        iSizeUpperLimit = iSizeUpperLimit,
        bqFilter = bqFilter,
        chr = chr,
        N = length(sampleNames),
        downsampleToCov = downsampleToCov,
        sampleNames = sampleNames,
        inputdir = tempdir,
        regionName = regionName,
        tempdir = tempdir,
        chrStart = chrStart,
        chrEnd = chrEnd,
        chrLength = NA,
        save_sampleReadsInfo = TRUE,
        use_bx_tag = use_bx_tag,
        bxTagUpperLimit = bxTagUpperLimit,
        default_sample_no_read_behaviour = "return_null"
    )

    load(file_sampleReads(tempdir, iSample, regionName))
    load(file_sampleReadsInfo(tempdir, iSample, regionName))
    removeTmpSamplesFile(tempdir, iSample, regionName, save_sampleReadsInfo = TRUE)

    if (length(sampleReads) < minimum_number_of_sample_reads) {
        print_message(paste0("Sample number ", iSample, " with sample name ", sampleNames[iSample], " has ", length(sampleReads), " reads which is fewer than the minimum ", minimum_number_of_sample_reads, ". This sample will therefore not be imputed and all results will be set to missing"))
        ## dummy up output here?
        ## note - useful, see writers.R,
        ## FORMAT <- "GT:GP:DS:HD" ## genotypes, posteriors, dosages, haploid-dosages
        per_sample_vcf_col <- "./.:.,.,.:.:.,."
        if (addOptimalHapsToVCF & have_truth_haplotypes) {
            FORMAT <- "GT:GP:DS:HD:OHD" ## (phased) genotypes, genotype posteriors, dosages, haploid-dosages, optimal haploid-dosages
            per_sample_vcf_col <- paste0(per_sample_vcf_col, ":.,.")
        }
        to_return <- list(
            sample_was_imputed = FALSE,
            per_sample_vcf_col = per_sample_vcf_col
        )
        return(to_return)
    }

    sample_alleleCount <- get_alleleCount(sampleReads, nrow(pos))
    print_message(paste0("The average depth of this sample is:", mean(sample_alleleCount[, 2])))
    print_message(paste0("There are ", length(sampleReads), " reads under consideration"))

    sampleReads <- snap_sampleReads_to_grid(
        sampleReads = sampleReads,
        grid = grid
    )

    i_gibbs_sample <- 1

    if (method == "nipt") {
        mat_dosage <- numeric(nSNPs)
        fet_dosage <- numeric(nSNPs)
        mat_gp_t <- array(0, c(3, nSNPs))
        fet_gp_t <- array(0, c(3, nSNPs))        
    } else {
        dosage <- numeric(nSNPs)        
        gp_t <- array(0, c(3, nSNPs))
    }
    nDosage <- 0
    nDosage_all <- 0

    wif0 <- as.integer(sapply(sampleReads, function(x) x[[2]]))
    grid_has_read <- rep(FALSE, nGrids)
    grid_has_read[wif0 + 1] <- TRUE

    nReads <- length(sampleReads)
    super_out_hap_dosages <- as.list(1:nGibbsSamples)
    super_out_read_labels <- as.list(1:nGibbsSamples)
    super_out_dosage_matrix <- as.list(1:nGibbsSamples)
    read_label_matrix_all <- array(NA, c(nReads, nGibbsSamples)) ## need this for phasing - store final read labels in it
    read_label_matrix_conf <- array(FALSE, c(nReads, nGibbsSamples)) ## need this for phasing - store whether these are confident reads

    return_gamma_t <- FALSE ## can turn on later, if required

    phasing_read_labels <- NULL
    phasing_haps <- NULL
    phasing_dosage <- NULL
    pse_mat <- NULL
    if (have_truth_haplotypes) {
        pse_mat <- array(NA, c(n_seek_its * (nGibbsSamples + 1), 6))
        colnames(pse_mat) <- c("ig", "ii", "pi", "r2", "pse", "disc")
        ## match
        s <- sampleNames[iSample]
        if (!(s %in% dimnames(phase)[[2]])) {
            stop("Something went wrong with phase naming")
        }
        if (method == "diploid") {
            truth_haps <- cbind(phase[, s, 1], phase[, s, 2])
        } else {
            truth_haps <- cbind(phase[, s, 1], phase[, s, 2], phase[, s, 3])
        }
    } else {
        truth_haps <- NULL
    }

    if (have_truth_genotypes) {
        s <- sampleNames[iSample]
        if (!(s %in% colnames(gen))) {
            stop("Something went wrong with gen naming")
        }
        truth_gen <- gen[, s, drop = FALSE]
    } else {
        truth_gen <- NULL
    }
    

    if (hla_run) {
        ## print_message("SPECIAL HLA CODE SIMON")
        gamma_total <- array(0, nrow(hapMatcher))
        list_of_gammas <- as.list(1:nGibbsSamples)
    }

    ## can just set this here
    if (is.null(phase) & !is.null(gen)) {
        have_truth_genotypes <- TRUE ## make UNIQUE to genotype having
    } else {
        have_truth_genotypes <- FALSE
    }

    ## leave NULL, add to it later if needed
    if (plot_per_sample_likelihoods) {
        for_likelihood_plotting <- as.list(1:(nGibbsSamples + 1))
    }

    if (plot_p1) {
        record_read_label_usage <- TRUE ## check this too
        p1_store <- lapply(1:(nGibbsSamples + 1), function(x) {
            as.list(1:n_seek_its)
        })
        read_store <- p1_store
        return_p1 <- TRUE
    } else {
        return_p1 <- FALSE
    }

    ## just in case
    hap3 <- NULL
    
    ## want this off for real usage (should be OK!)
    ## testing might want this off for consistency
    disable_read_category_usage <- FALSE

    ## don't need this for routine use - or do better matching!
    ## truth_g <- as.integer(truth_gen[, sampleNames[iSample]])
    for(i_gibbs_sample in 1:(nGibbsSamples + 1)) {

        if (K >100000) {
            gc(reset = TRUE);            gc(reset = TRUE);
        }

        if (i_gibbs_sample == (nGibbsSamples + 1)) {
            print_message("Phasing it")
            phasing_it <- TRUE
        } else {
            phasing_it <- FALSE
        }

        outplotprefix <- file.path(outputdir, "plots", paste0("haps.", sample_name, ".", regionName, "_igs.", i_gibbs_sample, "."))

        if (record_read_label_usage & !phasing_it) {
            ## init, truth, after each gibbs, after each seek
            cols <- c("minp", "maxp", "pos", "truth", "uncertain", "init", paste0("gibbs", 1:n_seek_its))
            read_label_matrix <- array(NA, c(length(sampleReads), length(cols)))
            colnames(read_label_matrix) <- cols
            rownames(read_label_matrix) <- sampleReadsInfo[, "qname"]
            read_label_matrix[, "pos"] <- L_grid[sapply(sampleReads, function(x) x[[2]]) + 1]
            x <- t(sapply(sampleReads, function(sampleRead) {L[range(sampleRead[[4]]) + 1]}))
            read_label_matrix[, "minp"] <- x[, 1]
            read_label_matrix[, "maxp"] <- x[, 2]
        }


        ##
        ## possibly do truth
        ## note: ignore this
        if (have_truth_haplotypes) {

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", truth"))
            }
            if (method == "diploid") {
                truth_label_set <- determine_a_set_of_truth_labels(
                    sampleReads = sampleReads,
                    truth_hap1 = truth_haps[, 1],
                    truth_hap2 = truth_haps[, 2],
                    maxDifferenceBetweenReads = maxDifferenceBetweenReads
		)
                truth_labels <- truth_label_set[["truth_labels"]]
                uncertain_truth_labels <- truth_label_set[["uncertain_truth_labels"]]
            } else {
                truth_label_set <- determine_a_set_of_truth_labels_for_nipt(
                    sampleReads = sampleReads,
                    truth_haps = truth_haps,
                    phase = phase,
                    ff = ff
                )
                truth_labels <- truth_label_set[["truth_labels"]]
                uncertain_truth_labels <- truth_label_set[["uncertain_truth_labels"]]
            }

            if ((i_gibbs_sample == 1) && estimate_bq_using_truth_read_labels) {
                bq_result <- estimate_bq(truth_labels = truth_labels, sampleReads = sampleReads, truth_haps = truth_haps)
                print(bq_result)
            }

            ## potentially for plotting
            previously_selected_haplotypes <- NULL

            truth_all <- impute_using_everything(
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                H = truth_labels,
                sampleReads = sampleReads,
                rhb_t = rhb_t,
                nSNPs = nSNPs,
                full_alphaHat_t = full_alphaHat_t,
                full_betaHat_t = full_betaHat_t,
                full_gamma_t = full_gamma_t,
                full_transMatRate_t_H  = full_transMatRate_t_H,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                hapMatcher = hapMatcher,
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                eMatDH_special_grid_which = eMatDH_special_grid_which,
                eMatDH_special_values_list = eMatDH_special_values_list,
                ref_error = ref_error,
                make_plots = make_plots,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                have_truth_genotypes = have_truth_genotypes,
                truth_haps = truth_haps,
                truth_gen = truth_gen,
                truth_labels = truth_labels,
                uncertain_truth_labels = uncertain_truth_labels,
                L_grid = L_grid,
                L = L,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                ancAlleleFreqAll = ancAlleleFreqAll,
                full_gammaSmall_t = full_gammaSmall_t,
                full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                plot_description = "0.truth",
                return_dosage = TRUE,
                sample_name = sample_name,
                smooth_cm = smooth_cm,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                minGLValue = minGLValue,
                suppressOutput = suppressOutput,
                method = method,
                Knew = Knew,
                previously_selected_haplotypes = previously_selected_haplotypes,
                use_eigen = use_eigen
            )

            if (verbose) {

                hap1 <- truth_all[["dosage1"]]
                hap2 <- truth_all[["dosage2"]]
                if (method == "nipt") {
                    hap3 <- truth_all[["dosage3"]]
                } else {
                    hap3 <- NULL
                }

                calculate_pse_and_r2_master(
                    method = method,
                    have_truth_haplotypes = have_truth_haplotypes,
                    have_truth_genotypes = have_truth_genotypes,                    
                    truth_haps = truth_haps,
                    truth_gen = truth_gen,
                    hap1 = hap1,
                    hap2 = hap2,
                    hap3 = hap3,
                    impute_rare_common = impute_rare_common,
                    checking_all_snps = FALSE,
                    verbose = verbose,
                    af = af,
                    inRegion2 = inRegion2
                )

                rm(hap1, hap2, hap3)
            }

        } else {
            truth_labels <- NULL
            uncertain_truth_labels <- NULL
        }


        ##
        ## prepare output matrix
        ##
        if (record_interim_dosages) {
            cols <- paste0(rep(c("gibbs", "all"), n_seek_its), rep(1:n_seek_its, each = 2))
            if (have_truth_haplotypes) {
                cols <- c("truth", cols)
            }
            dosage_matrix <- array(NA, c(nSNPs, length(cols)))
            colnames(dosage_matrix) <- cols
            if (have_truth_haplotypes) {
                dosage_matrix[, "truth"] <- truth_all[["dosage1"]] + truth_all[["dosage2"]]
            }
        }

        if (record_read_label_usage) {
            if (have_truth_haplotypes) {
                read_label_matrix[, "truth"] <- truth_labels
                read_label_matrix[, "uncertain"] <- as.numeric(uncertain_truth_labels)
            }
        }
        i_it <- 1

        ##
        ## now loop, first on subset, then on all
        ##
        for(i_it in 1:n_seek_its) {

            ## here it is 1
            n_gibbs_starts <- 1
            if (i_it == 1 & !phasing_it) {
                which_haps_to_use <- sort(sample(1:K, Ksubset))
                double_list_of_starting_read_labels <- list(
                    lapply(1:n_gibbs_starts, function(i) {
                        if (method == "diploid") {
                            H <- sample(c(1, 2), length(sampleReads), replace = TRUE)
                        } else {
                            H <- sample(c(1, 2, 3), prob = c(0.5, 0.5 - ff / 2, ff / 2),length(sampleReads), replace = TRUE)
                        }
                        return(H)
                    })
                )
                gibbs_initialize_iteratively <- TRUE
            } else {
                ## if phasing it will suck these up! same for "which_haps_to_use" (though is that bad?)
                ## which_haps_to_use <- which_haps_to_use
                double_list_of_starting_read_labels <- list(list(read_labels))
                gibbs_initialize_iteratively <- FALSE
            }

            if (i_it == 1 & record_read_label_usage) {
                read_label_matrix[, "init"] <- double_list_of_starting_read_labels[[1]][[1]]
            }

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " small gibbs"))
            }

            if (record_interim_dosages) {
                return_genProbs <- TRUE
                return_hapProbs <- TRUE
            } else {
                return_genProbs <- FALSE
                return_hapProbs <- FALSE
            }
            
            ## print_message(paste0("i_gibbs, which_haps_to_use = ", paste(which_haps_to_use, collapse = ",")))
            if (use_mspbwt) {
                return_hapProbs <- TRUE
            }

            gibbs_iterate <- impute_one_sample(
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                distinctHapsB = distinctHapsB,
                distinctHapsIE = distinctHapsIE,
                hapMatcher = hapMatcher,
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                rhb_t = rhb_t,
                ref_error = ref_error,
                nSNPs = nSNPs,
                sampleReads = sampleReads,
                ff = ff,
                small_eHapsCurrent_tc = small_eHapsCurrent_tc,
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
                smooth_cm = smooth_cm,
                which_haps_to_use = which_haps_to_use,
                n_gibbs_starts = n_gibbs_starts,
                small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
                n_gibbs_sample_its = 1,
                double_list_of_starting_read_labels = double_list_of_starting_read_labels,
                small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                perform_block_gibbs = TRUE,
                make_plots = make_plots,
                wif0 = wif0,
                grid_has_read = grid_has_read,
                plot_description = paste0("it", i_it, ".gibbs"),
                ancAlleleFreqAll = ancAlleleFreqAll,
                grid = grid,
                L_grid = L_grid,
                L = L,
                inRegion2 = inRegion2,
                cM_grid = cM_grid,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                truth_haps = truth_haps,
                have_truth_genotypes = have_truth_genotypes,
                truth_gen = truth_gen,
                truth_labels = truth_labels,
                uncertain_truth_labels = uncertain_truth_labels,
                verbose = FALSE,
                maxEmissionMatrixDifference = 1e100,
                return_p_store = FALSE,
                return_p1 = return_p1,
                return_extra = FALSE,
                return_genProbs = return_genProbs,
                return_hapProbs = return_hapProbs,
                return_gibbs_block_output = FALSE,
                gibbs_initialize_iteratively = gibbs_initialize_iteratively,
                gibbs_initialize_at_first_read = FALSE,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                rescale_eMatRead_t = TRUE,
                rescale_eMatGrid_t = FALSE,
                Jmax = 10000,
                suppressOutput = suppressOutput,
                shuffle_bin_radius = shuffle_bin_radius,
                make_plots_block_gibbs = make_plots_block_gibbs,
                sample_name = sample_name,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
                method = method,
                use_sample_is_diploid = use_sample_is_diploid,
                small_ref_panel_skip_equally_likely_reads = small_ref_panel_skip_equally_likely_reads,
                small_ref_panel_equally_likely_reads_update_iterations = small_ref_panel_equally_likely_reads_update_iterations,
                i_it = i_it,
                i_gibbs_sample = i_gibbs_sample,
                shard_check_every_pair = shard_check_every_pair,
                calculate_gamma_on_the_fly = calculate_gamma_on_the_fly,
                disable_read_category_usage = disable_read_category_usage
            )

            if (plot_p1) {
                p1_store[[i_gibbs_sample]][[i_it]] <- gibbs_iterate[["p1"]]
                read_store[[i_gibbs_sample]][[i_it]] <- gibbs_iterate[["pH"]]
            }

            if (hla_run) {
                ## final phasing it, save gamma
                if (
                (i_it == n_seek_its)
                ) {
                    ## print_message("HLA SPECIAL CODE SIMON")
                    return_gamma_t <- TRUE
                    full_gamma_t <- array(0, c(K, nGrids))
                } else {
                    return_gamma_t <- FALSE
                }
            }

            if (plot_per_sample_likelihoods) {
                if (i_it == 1) {
                    for_likelihood_plotting[[i_gibbs_sample]] <- as.list(1:n_seek_its)
                }
                for_likelihood_plotting[[i_gibbs_sample]][[i_it]] <- gibbs_iterate$per_it_likelihoods
            }

            ## continue to work on this. weird error below hmm...
            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " full"))
            }

            if (i_it < n_seek_its) {
                return_good_haps <- TRUE
            } else {
                ## keep on, keep selecting good haps, phasing sucks these up
                return_good_haps <- TRUE
            }

            read_labels <- gibbs_iterate$double_list_of_ending_read_labels[[1]][[1]]
            previously_selected_haplotypes <- sample(which_haps_to_use, Ksubset - Knew)
            # is equavaliant to hapProbs_t
            return_dosage <- (have_truth_haplotypes | record_interim_dosages | (i_it > n_burn_in_seek_its))

            igibbs <- (i_gibbs_sample - 1) * n_seek_its + i_it ## 1-based

            which_haps_to_use_zilong_A <- NULL
            which_haps_to_use_zilong_B <- NULL
            which_haps_to_use_zilong_mspbwt <- NULL
            which_haps_to_use_quilt1 <- NULL
            
            if ((!use_mspbwt && make_heuristic_plot)) {

                ##print(paste0("zilong = ", zilong))
                ##print(paste0("(!use_mspbwt && make_heuristic_plot) = ", (!use_mspbwt && make_heuristic_plot)))

                Kfull <- nrow(hapMatcher)
                hap1 <- gibbs_iterate$hapProbs_t[1, ]
                hap2 <- gibbs_iterate$hapProbs_t[2, ]
                if (method == "nipt") {
                    hap3 <- gibbs_iterate$hapProbs_t[3, ]
                }
                
                if (heuristic_approach == "A" | make_heuristic_plot) {

                    which_haps_to_use_zilong_A <- which_haps_to_use
                    
                }

                if (heuristic_approach == "B" | make_heuristic_plot) {


                    which_haps_to_use_zilong_B <- which_haps_to_use

                }
                
            }

            if (use_mspbwt | (make_heuristic_plot)) {

                ## for testing purposes
                if (use_splitreadgl) {
                    impute_all <- impute_using_split_reads_and_small_ref_panel(
                        H = read_labels,
                        which_haps_to_use = which_haps_to_use,
                        sampleReads = sampleReads,
                        rhb_t = rhb_t,
                        nSNPs = nSNPs,
                        full_alphaHat_t = full_alphaHat_t,
                        full_betaHat_t = full_betaHat_t,
                        full_gamma_t = full_gamma_t,
                        full_gammaSmall_t = full_gammaSmall_t,
                        full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                        full_transMatRate_t_H = full_transMatRate_t_H,
                        distinctHapsB = distinctHapsB,
                        distinctHapsIE = distinctHapsIE,
                        hapMatcher = hapMatcher,
                        eMatDH_special_grid_which = eMatDH_special_grid_which,
                        eMatDH_special_values_list = eMatDH_special_values_list,
                        ref_error = ref_error,
                        make_plots = FALSE, ## these plots are pretty useless? as plots?
                        outplotprefix = outplotprefix,
                        have_truth_haplotypes = have_truth_haplotypes,
                        truth_haps = truth_haps,
                        truth_labels = truth_labels,
                        uncertain_truth_labels = uncertain_truth_labels,
                        L_grid = L_grid,
                        L = L,
                        inRegion2 = inRegion2,
                        cM_grid = cM_grid,
                        ancAlleleFreqAll = ancAlleleFreqAll,
                        plot_description = paste0("it", i_it, ".full"),
                        return_good_haps = return_good_haps,
                        Knew = Knew,
                        return_dosage  = return_dosage,
                        previously_selected_haplotypes = previously_selected_haplotypes,
                        K_top_matches = K_top_matches,
                        heuristic_match_thin = heuristic_match_thin,
                        return_gamma_t = return_gamma_t,
                        sample_name = sample_name,
                        smooth_cm = smooth_cm ,
                        regionStart = regionStart,
                        regionEnd = regionEnd,
                        buffer = buffer,
                        minGLValue = minGLValue,
                        suppressOutput = suppressOutput
                    )
                    hap1 <- impute_all[["dosage1"]]
                    hap2 <- impute_all[["dosage2"]]
                    if (method == "nipt") {
                        hap3 <- impute_all[["dosage3"]]
                    }
                } else {
                    hap1 <- gibbs_iterate$hapProbs_t[1, ]
                    hap2 <- gibbs_iterate$hapProbs_t[2, ]
                    if (method == "nipt") {
                        hap3 <- gibbs_iterate$hapProbs_t[3, ]
                    }
                }

                ## this version here
                if (method == "diploid") {
                    hapProbs_t <- rbind(hap1, hap2)
                } else {
                    hapProbs_t <- rbind(hap1, hap2, hap3)
                }
                Kfull <- nrow(hapMatcher)

                if (heuristic_approach == "A" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_mspbwt_v3(
                        hapProbs_t = hapProbs_t,
                        hapMatcher = hapMatcher,
                        hapMatcherR = hapMatcherR,
                        use_hapMatcherR = use_hapMatcherR,
                        ms_indices = ms_indices,
                        Knew = Knew,
                        Kfull = Kfull,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        heuristic_approach = "A",
                        method = method
                    )
                    which_haps_to_use_mspbwt_A <- which_haps_to_use

                }

                if (heuristic_approach == "B" | make_heuristic_plot) {

                    which_haps_to_use <- select_new_haps_mspbwt_v3(
                        hapProbs_t = hapProbs_t,
                        hapMatcher = hapMatcher,
                        hapMatcherR = hapMatcherR,
                        use_hapMatcherR = use_hapMatcherR,
                        ms_indices = ms_indices,
                        Knew = Knew,
                        Kfull = Kfull,
                        mspbwtL = mspbwtL,
                        mspbwtM = mspbwtM,
                        heuristic_approach = "B",
                        method = method                        
                    )
                    which_haps_to_use_mspbwt_B <- which_haps_to_use

                }
                
                

            }

            if ((!use_mspbwt) | make_heuristic_plot) {
                
                impute_all <- impute_using_everything(
                    eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                    eMatDH_special_matrix = eMatDH_special_matrix,
                    use_eMatDH_special_symbols = use_eMatDH_special_symbols,
                    H = read_labels,
                    sampleReads = sampleReads,
                    rhb_t = rhb_t,
                    nSNPs = nSNPs,
                    full_alphaHat_t = full_alphaHat_t,
                    full_betaHat_t = full_betaHat_t,
                    full_gamma_t = full_gamma_t,
                    full_gammaSmall_t = full_gammaSmall_t,
                    full_gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
                    full_transMatRate_t_H = full_transMatRate_t_H,
                    distinctHapsB = distinctHapsB,
                    distinctHapsIE = distinctHapsIE,
                    hapMatcher = hapMatcher,
                    hapMatcherR = hapMatcherR,
                    use_hapMatcherR = use_hapMatcherR,
                    eMatDH_special_grid_which = eMatDH_special_grid_which,
                    eMatDH_special_values_list = eMatDH_special_values_list,
                    ref_error = ref_error,
                    make_plots = FALSE, ## these plots are pretty useless? as plots?
                    outplotprefix = outplotprefix,
                    have_truth_haplotypes = have_truth_haplotypes,
                    have_truth_genotypes = have_truth_genotypes,
                    truth_haps = truth_haps,
                    truth_gen = truth_gen,
                    truth_labels = truth_labels,
                    uncertain_truth_labels = uncertain_truth_labels,
                    L_grid = L_grid,
                    L = L,
                    inRegion2 = inRegion2,
                    cM_grid = cM_grid,
                    ancAlleleFreqAll = ancAlleleFreqAll,
                    plot_description = paste0("it", i_it, ".full"),
                    return_good_haps = return_good_haps,
                    Knew = Knew,
                    return_dosage  = return_dosage,
                    previously_selected_haplotypes = previously_selected_haplotypes,
                    K_top_matches = K_top_matches,
                    heuristic_match_thin = heuristic_match_thin,
                    return_gamma_t = return_gamma_t,
                    sample_name = sample_name,
                    smooth_cm = smooth_cm ,
                    regionStart = regionStart,
                    regionEnd = regionEnd,
                    buffer = buffer,
                    minGLValue = minGLValue,
                    suppressOutput = suppressOutput,
                    use_eigen = use_eigen,
                    method = method
                )
                which_haps_to_use <- c(previously_selected_haplotypes, impute_all$new_haps)
                hap1 <- impute_all[["dosage1"]]
                hap2 <- impute_all[["dosage2"]]
                if (method == "nipt") {
                    hap3 <- impute_all[["dosage3"]]
                }
                which_haps_to_use_quilt1 <- which_haps_to_use
            }

            if (make_heuristic_plot) {

                hapProbs_t <- gibbs_iterate$hapProbs_t 
                
                compare_heuristic_approaches(
                    hapProbs_t = hapProbs_t,
                    which_haps_to_use_quilt1 = which_haps_to_use_quilt1,
                    which_haps_to_use_zilong_A = which_haps_to_use_zilong_A,
                    which_haps_to_use_zilong_B = which_haps_to_use_zilong_B,
                    which_haps_to_use_mspbwt_A = which_haps_to_use_mspbwt_A,
                    which_haps_to_use_mspbwt_B = which_haps_to_use_mspbwt_B,                    
                    hapMatcherR = hapMatcherR,
                    hapMatcher = hapMatcher,
                    use_hapMatcherR = use_hapMatcherR,
                    distinctHapsB = distinctHapsB,
                    eMatDH_special_matrix = eMatDH_special_matrix,
                    eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                    outputdir = outputdir,
                    sample_name = sample_name,
                    regionName = regionName,
                    i_gibbs_sample = i_gibbs_sample,
                    i_it = i_it,
                    nGrids = nGrids
                )

            }

            
            if (record_interim_dosages) {
                dosage_matrix[, paste0("gibbs", i_it)] <- get_dosages_from_fbsoL(gibbs_iterate)
                dosage_matrix[, paste0("all", i_it)] <- hap1 + hap2
                stopifnot(method == "diploid") ## else fix above
            }

            if (record_read_label_usage) {
                read_label_matrix[, paste0("gibbs", i_it)] <- read_labels
            }

            ## if not a phasing iteration, and past burn in, store!
            if (!phasing_it && (i_it > n_burn_in_seek_its)) {

                if (method == "diploid") {

                    dosage <- dosage + hap1 + hap2
                    gp_t <- gp_t +
                        rbind((1 - hap1) * (1 - hap2), (1 - hap1) * hap2 + hap1 * (1 - hap2), hap1 * hap2)
                    nDosage <- nDosage + 1

                } else {
                    
                    mat_dosage <- mat_dosage + hap1 + hap2
                    fet_dosage <- fet_dosage + hap1 + hap3
                    mat_gp_t <- mat_gp_t + 
                        rbind((1 - hap1) * (1 - hap2), (1 - hap1) * hap2 + hap1 * (1 - hap2), hap1 * hap2)
                    fet_gp_t <- fet_gp_t + 
                        rbind((1 - hap1) * (1 - hap3), (1 - hap1) * hap3 + hap1 * (1 - hap3), hap1 * hap3)
                    nDosage <- nDosage + 1

                }

            }

            ## regardless, print out accuracy as going along
            calculate_pse_and_r2_master(
                method = method,
                have_truth_haplotypes = have_truth_haplotypes,
                have_truth_genotypes = have_truth_genotypes,                    
                truth_haps = truth_haps,
                truth_gen = truth_gen,
                hap1 = hap1,
                hap2 = hap2,
                hap3 = hap3,
                impute_rare_common = impute_rare_common,
                checking_all_snps = FALSE,
                verbose = verbose,
                af = af,
                inRegion2 = inRegion2
            )

        }


        if (impute_rare_common) {

            if (verbose) {
                print_message(paste0("i_gibbs=", i_gibbs_sample, ", i_it = ", i_it, " small gibbs (all SNPs)"))
            }

            out_rare_common <- impute_final_gibbs_with_rare_common(
                special_rare_common_objects = special_rare_common_objects,
                special_rare_common_objects_per_core = special_rare_common_objects_per_core,
                allSNP_sampleReads = allSNP_sampleReads,
                hap1 = hap1,
                hap2 = hap2,
                hap3 = hap3,
                pos_all = pos_all,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                hapMatcher = hapMatcher,                
                hapMatcherR = hapMatcherR,
                use_hapMatcherR = use_hapMatcherR,
                distinctHapsIE = distinctHapsIE,
                eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
                eMatDH_special_matrix = eMatDH_special_matrix,
                Ksubset = Ksubset,
                ref_error = ref_error,
                which_haps_to_use = which_haps_to_use,
                small_ref_panel_gibbs_iterations = small_ref_panel_gibbs_iterations,
                small_ref_panel_block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                allSNP_wif0 = allSNP_wif0,
                allSNP_grid_has_read = allSNP_grid_has_read,
                make_plots = make_plots,
                outplotprefix = outplotprefix,
                have_truth_haplotypes = have_truth_haplotypes,
                truth_haps_all = truth_haps_all,
                have_truth_genotypes = have_truth_genotypes,
                truth_gen_all = truth_gen_all,
                truth_labels_all = truth_labels_all,
                uncertain_truth_labels_all = uncertain_truth_labels_all,
                shuffle_bin_radius = shuffle_bin_radius,
                make_plots_block_gibbs = make_plots_block_gibbs,
                sample_name = sample_name,
                regionStart = regionStart,
                regionEnd = regionEnd,
                buffer = buffer,
                i_it = i_it,
                i_gibbs_sample = i_gibbs_sample,
                shard_check_every_pair = shard_check_every_pair,
                sampleNames = sampleNames,
                iSample = iSample,
                phase_all = phase_all,
                ff = ff,
                method = method
            )

            hap1_all <- out_rare_common[["hap1"]]
            hap2_all <- out_rare_common[["hap2"]]
            if (method == "nipt") {
                hap3_all <- out_rare_common[["hap3"]]
            }

            if (!phasing_it && (i_it > n_burn_in_seek_its)) {
                
                if (method == "diploid") {
                    
                    dosage_all <- dosage_all + hap1_all + hap2_all
                    gp_t_all <- gp_t_all + rbind(
                    (1 - hap1_all) * (1 - hap2_all),
                    (1 - hap1_all) * hap2_all + hap1_all * (1 - hap2_all),
                    hap1_all * hap2_all
                    )
                    nDosage_all <- nDosage_all + 1

                } else {

                    mat_dosage_all <- mat_dosage_all + hap1_all + hap2_all
                    fet_dosage_all <- fet_dosage_all + hap1_all + hap3_all ## was typo, second was hap2_all
                    mat_gp_t_all <- mat_gp_t_all + 
                        rbind((1 - hap1_all) * (1 - hap2_all), (1 - hap1_all) * hap2_all + hap1_all * (1 - hap2_all), hap1_all * hap2_all)
                    fet_gp_t_all <- fet_gp_t_all + 
                        rbind((1 - hap1_all) * (1 - hap3_all), (1 - hap1_all) * hap3_all + hap1_all * (1 - hap3_all), hap1_all * hap3_all)
                    nDosage_all <- nDosage_all + 1

                }

                calculate_pse_and_r2_master(
                    method = method,
                    have_truth_haplotypes = have_truth_haplotypes,
                    have_truth_genotypes = have_truth_genotypes,                    
                    truth_haps = truth_haps_all,
                    truth_gen = truth_gen_all,
                    hap1 = hap1_all,
                    hap2 = hap2_all,
                    hap3 = hap3_all,
                    impute_rare_common = TRUE, ## has to be true here!
                    checking_all_snps = TRUE,
                    verbose = verbose,
                    af = special_rare_common_objects[["ref_alleleCount_all"]][, 3],
                    inRegion2 = special_rare_common_objects[["inRegion2"]]
                )
                
            }
            
        }

        if (!phasing_it) {
            ## for phasing bit
            read_label_matrix_all[, i_gibbs_sample] <- read_labels
            read_label_matrix_conf[, i_gibbs_sample] <- assess_ability_of_reads_to_be_confident(
                hap1 = hap1,
                hap2 = hap2,
                hap3 = hap3,
                sampleReads = sampleReads,
                maxDifferenceBetweenReads = maxDifferenceBetweenReads,
                minrp = 0.95,
                minmp = 0.95,
                method = method
            )
            if (record_read_label_usage) {
                super_out_read_labels[[i_gibbs_sample]] <- read_label_matrix
            }
            if (record_interim_dosages) {
                super_out_dosage_matrix[[i_gibbs_sample]] <- dosage_matrix
            }
        }

        ## save last
        if (i_gibbs_sample == nGibbsSamples) {
            can_hap <- nGibbsSamples
            nReads <- length(sampleReads)
            if (method == "diploid") {
                out_best_labels <- determine_best_read_label_so_far(
                    read_label_matrix_all = read_label_matrix_all,
                    read_label_matrix_conf = read_label_matrix_conf,                
                    nReads = nReads,
                    nGibbsSamples = nGibbsSamples,
                    verbose = verbose,
                    can_hap = can_hap
                )
                read_labels <- out_best_labels$read_labels
                if (verbose) {
                    x <- out_best_labels$flip_matrix[, can_hap]
                    if (length(x) > 0) {
                        print_message(paste0("There are ", sum(x), " out of ", length(x), " regions that have been flipped by consensus"))
                    }
                }
            } else {
                ## this is not ideal
                ## ideally re-visit once I've written a bit more code
                out_best_labels <- determine_best_read_label_so_far_nipt(
                    read_label_matrix_all = read_label_matrix_all,
                    read_label_matrix_conf = read_label_matrix_conf,                
                    nReads = nReads,
                    nGibbsSamples = nGibbsSamples,
                    verbose = verbose,
                    can_hap = can_hap
                )
                read_labels <- out_best_labels$read_labels
            }
        }

        if (phasing_it) {
            if (!impute_rare_common) {
                if (method == "diploid") {                    
                    ## just save relevant stuff here
                    out <- recast_haps(hd1 = hap1, hd2 = hap2, gp = t(gp_t))
                    ## over-write here
                    hap1 <- out$hd1
                    hap2 <- out$hd2
                    ## 
                    phasing_haps <- cbind(hap1, hap2)
                    phasing_dosage <- hap1 + hap2
                } else {
                    ## phasing_haps <- cbind(hap1, hap2, hap3)
                    out <- recast_nipt_haps(
                        hap1 = hap1,
                        hap2 = hap2,
                        hap3 = hap3,
                        mat_gp_t = mat_gp_t,
                        fet_gp_t = fet_gp_t
                    )
                    phasing_haps <- cbind(
                        out[["hap1"]],
                        out[["hap2"]],
                        out[["hap3"]]
                    )
                }
            } else {
                if (method == "diploid") {                
                    ## just save relevant stuff here
                    out <- recast_haps(hd1 = hap1_all, hd2 = hap2_all, gp = t(gp_t_all))
                    ## over-write here
                    hap1_all <- out$hd1
                    hap2_all <- out$hd2
                    phasing_haps_all <- cbind(hap1_all, hap2_all)
                } else {
                    out <- recast_nipt_haps(
                        hap1 = hap1_all,
                        hap2 = hap2_all,
                        hap3 = hap3_all,
                        mat_gp_t = mat_gp_t_all,
                        fet_gp_t = fet_gp_t_all
                    )
                    ## phasing_haps <- cbind(hap1_all, hap2_all, hap3_all)
                    ##phasing_mat_dosage <- hap1_all + hap2_all
                    ##phasing_fet_dosage <- hap1_all + hap3_all
                    phasing_haps_all <- cbind(
                        out[["hap1"]],
                        out[["hap2"]],
                        out[["hap3"]]
                    )
                }
            }
        }

        if (hla_run) {
            gamma1 <- impute_all$fbsoL$gammaMT_t
            gamma2 <- impute_all$fbsoL$gammaMU_t
            if (is.na(gamma_physically_closest_to)) {
                iGrid <- round(nGrids / 2)
            } else {
                iGrid <- grid[which.min(abs(L - gamma_physically_closest_to))] + 1
            }
            gamma1 <- gamma1[, iGrid]
            gamma2 <- gamma2[, iGrid]
            if (!phasing_it) {
                gamma_total <- gamma_total + gamma1 + gamma2
            }
            if (i_gibbs_sample <= nGibbsSamples) {
                list_of_gammas[[i_gibbs_sample]] <- list(gamma1, gamma2)
            }
        }

    }

    if (plot_p1) {
        ##save(p1_store, read_store, super_out_read_labels, file = "~/temp.RData")
        p1 <- p1_store[[1]][[1]]
        print(paste0(sum(p1 == 0) / (prod(dim(p1)))))
        if (small_ref_panel_skip_equally_likely_reads) {
            save(p1_store, file = "~/temp.RData")
        }
    }
    

    if(have_truth_haplotypes) {
        imputed_truth_haplotypes <- cbind(truth_all[["dosage1"]], truth_all[["dosage2"]])
    } else {
        imputed_truth_haplotypes <- NULL
    }


    ##
    ## perform normalization and switchover here
    ##
    if (method == "diploid") {
        if (impute_rare_common) {
            dosage <- dosage_all / nDosage_all
            gp_t <- gp_t_all / nDosage_all
        } else {
            dosage <- dosage / nDosage
            gp_t <- gp_t / nDosage
        }
    } else {
        if (impute_rare_common) {
            mat_dosage <- mat_dosage_all / nDosage_all
            fet_dosage <- fet_dosage_all / nDosage_all
            mat_gp_t <- mat_gp_t_all / nDosage_all
            fet_gp_t <- fet_gp_t_all / nDosage_all
            
        } else {
            mat_dosage <- mat_dosage / nDosage
            fet_dosage <- fet_dosage / nDosage
            mat_gp_t <- mat_gp_t / nDosage
            fet_gp_t <- fet_gp_t / nDosage
        }
    }
    if (impute_rare_common) {
        phasing_haps <- phasing_haps_all    
        sampleReads <- allSNP_sampleReads
        nSNPs <- nrow(phasing_haps)
        af <- special_rare_common_objects[["ref_alleleCount_all"]][, 3]
        inRegion2 <- special_rare_common_objects[["inRegion2"]]
        truth_haps <- truth_haps_all
        truth_gen <- truth_gen_all
    }

    
    ##
    ## print final accuracies
    ##
    hap1 <- phasing_haps[, 1]
    hap2 <- phasing_haps[, 2]
    if (method == "nipt") {
        hap3 <- phasing_haps[, 3]
    }
    calculate_pse_and_r2_master(
        method = method,
        have_truth_haplotypes = have_truth_haplotypes,
        have_truth_genotypes = have_truth_genotypes,
        truth_haps = truth_haps,
        truth_gen = truth_gen,
        hap1 = hap1,
        hap2 = hap2,
        hap3 = hap3,
        impute_rare_common = impute_rare_common,
        checking_all_snps = TRUE,
        verbose = verbose,
        inRegion2 = inRegion2,
        af = af,
        dosage = dosage,
        mat_dosage = mat_dosage,
        fet_dosage = fet_dosage,
        prefix = paste0("Final imputation accuracy for sample ", sample_name, " ")
    )

    ## optionally plot here
    if (plot_per_sample_likelihoods) {
        plot_of_likelihoods_across_samplings_and_seek_its(
            outputdir = outputdir,
            for_likelihood_plotting = for_likelihood_plotting,
            nGibbsSamples = nGibbsSamples,
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            block_gibbs_iterations = block_gibbs_iterations,
            n_seek_its = n_seek_its,
            sample_name = sample_name,
            regionName = regionName,
            truth_likelihoods = NULL,
            ylab = "p_H_given_O_L_up_to_C"
        )
    }


    ##
    ## for allele count, from buildAlleleCount_subfunction in STITCH
    ##
    a <- unlist(sapply(sampleReads,function(x) x[[3]]))
    b <- unlist(sapply(sampleReads,function(x) x[[4]]))
    bqProbs <- STITCH::convertScaledBQtoProbs(matrix(a,ncol=1))
    ## y is numeric, z = integer counts (can have 0),
    c1 <- increment2N(
        y = as.numeric(bqProbs[,1]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    c2 <- increment2N(
        y = as.numeric(bqProbs[,2]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    per_sample_alleleCount <- cbind(c2, c1 + c2)


    ## make vcf character thing - NOTE - this is a HACK of previous one - but works OK here!
    ##
    ## 
    ##
    if (method == "diploid") {
        eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
        fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##
        max_gen <- get_max_gen_rapid(gp_t)
    } else {
        ## argh, just do mat for now
        eij <- round(mat_gp_t[2, ] + 2 * mat_gp_t[3, ], 3) ## prevent weird rounding issues
        fij <- round(mat_gp_t[2, ] + 4 * mat_gp_t[3, ], 3) ##
        max_gen <- get_max_gen_rapid(mat_gp_t)
    }
    if (method == "diploid") {
        if (addOptimalHapsToVCF & have_truth_haplotypes) {
            add_x_2_cols <- TRUE
            x_t <- t(imputed_truth_haplotypes)
        } else {
            add_x_2_cols <- FALSE
            x_t <- matrix()
        }
        per_sample_vcf_col <- STITCH::rcpp_make_column_of_vcf(
            gp_t = gp_t,
            use_read_proportions = FALSE,
            use_state_probabilities = TRUE,
            read_proportions = matrix(),
            q_t = t(phasing_haps),
            add_x_2_cols = add_x_2_cols,
            x_t = x_t
        )
        ## annoying but shouldn't be that slow I hope
        if (output_gt_phased_genotypes) {
            per_sample_vcf_col <-
                paste0(
                    round(phasing_haps[, 1]), "|", round(phasing_haps[, 2]),
                    substring(per_sample_vcf_col, first = 4, last = 100L)
                )
        }
    } else {

        ## do the slower in R version for now
        FORMAT <- "GT:MGP:MDS:FGP:FDS"
        ## so first want phased haplotypes
        ## next want first 
        per_sample_vcf_col <- paste0(
            round(phasing_haps[, 1]), "|",
            round(phasing_haps[, 2]), "|",
            round(phasing_haps[, 3]), ":",
            round(mat_gp_t[1, ], 3), ",",
            round(mat_gp_t[2, ], 3), ",",
            round(mat_gp_t[3, ], 3), ":",
            round(mat_dosage, 3), ":",
            round(fet_gp_t[1, ], 3), ",",
            round(fet_gp_t[2, ], 3), ",",
            round(fet_gp_t[3, ], 3), ":",
            round(fet_dosage, 3)
        )
        
    }

    ## what to return
    to_return <- list(
        sample_was_imputed = TRUE,
        eij = eij,
        fij = fij,
        max_gen = max_gen,
        per_sample_alleleCount = per_sample_alleleCount,
        per_sample_vcf_col = per_sample_vcf_col,
        super_out_hap_dosages = super_out_hap_dosages,
        super_out_read_labels = super_out_read_labels,
        super_out_dosage_matrix = super_out_dosage_matrix
    )

    ## phasing_read_labels = phasing_read_labels,
    ## imputed_truth_haplotypes = imputed_truth_haplotypes,
    ## pse_mat = pse_mat,
    ## dosage = dosage,
    ## super_out_dosage_matrix = super_out_dosage_matrix,
    ## phasing_haps = phasing_haps,
    ## phasing_dosage = phasing_dosage,
    ## gp_t = gp_t,
    ##

    if (hla_run) {
        ## print_message("More HLA SIMON code")
        to_return[["gamma1"]] <- gamma1
        to_return[["gamma2"]] <- gamma2
        to_return[["gamma_total"]] <- gamma_total
        to_return[["list_of_gammas"]] <- list_of_gammas
    }

    return(to_return)


}


#' @export
modified_calculate_pse <- function(
    test,
    truth,
    LL,
    seed = NULL
) {
    ## for testing
    if (is.null(seed) == FALSE)
        set.seed(seed)
    rownames(test) <- 1:nrow(test)
    which_sites <-
        rowSums(truth == 0 | truth == 1, na.rm = TRUE) == 2 &
        rowSums(truth, na.rm = TRUE) == 1 &
        rowSums(is.na(truth)) == 0
    truth <- truth[which_sites, ]
    test <- test[which_sites, ]
    if (nrow(test) == 0)
        return(NA)
    ## as these sites, discrepency as well
    disc <- sum(round(rowSums(test)) != 1)
    test <- round(test)
    ## round test data for now. choose hets at random
    ## for homs
    test[, 1] <- as.integer(round(test[, 1]))
    test[, 2] <- as.integer(round(test[, 2]))
    testO <- test
    truthO <- truth
    for(i_option in 1:2) {
        test <- testO
        truth <- truthO
        ## specifically remove from consideration double phase switch errors
        w <- rowSums(test) == 1
        w2 <- which(diff(abs(test[w,1] - truth[w,1])) != 0)
        to_remove <- NULL
        if (length(w2) > 0) {
            w3 <- which(diff(w2) == 1)
            if (length(w3) > 0) {
                for(a in w3) {
                    c <- w2[c(a, a + 1)]
                    to_remove <- c(to_remove, which(w)[c])
                }
            }
        }
        ## double pse are two consecutive
        if (length(to_remove) > 0) {
            test <- test[-to_remove, ]
            truth <- truth[-to_remove, ]
        }
        ##
        if (i_option == 1) {
            ## only consider non-discrepent sites
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            w <- rowSums(test) == 1
            if (sum(w) == 0) {
                print("Test has no hets! possibly an error or homo over region, possibly no record dosages turned on in impute_all")
                switches1 <- cbind(i1 = NA, i2 = NA, l1 = NA, l2 = NA)
                phase_errors_def1 <- 0
                phase_sites_def1 <- 0
            } else {
                y <- diff(abs(test[w,1] - truth[w,1])) != 0
                phase_errors_def1 <- sum(y)
                phase_sites_def1 <- sum(w) - 1
                s <- as.integer(rownames(test[w, , drop = FALSE][c(as.logical(y), FALSE), , drop = FALSE]))
                e <- as.integer(rownames(test[w, , drop = FALSE][c(FALSE, as.logical(y)), , drop = FALSE]))
                switches1 <- cbind(i1 = s, i2 = e, l1 = LL[s], l2 = LL[e])
            }
        }
        if (i_option == 2) {
            choose_at_random <- which(rowSums(test) != 1)
            if (length(choose_at_random) > 0) {
                test[choose_at_random, ] <- 0
                r <- sample(
                    c(1, 2),
                    length(choose_at_random),
                    replace = TRUE
                )
                test[cbind(choose_at_random, r)] <- 1
            }
            ## chose best start
            if (test[1, 1] != truth[1, 1])
                test <- test[, c(2, 1)]
            ## calculate number of differences
            phase_errors_def2 <- sum(diff(abs(test[,1] - truth[,1])) != 0)
            phase_sites_def2 <- nrow(test) - 1
        }
    }
    ##
    return(
        list(
            values = c(
                phase_errors_def1 = phase_errors_def1,
                phase_sites_def1 = phase_sites_def1,
                phase_errors_def2 = phase_errors_def2,
                phase_sites_def2 = phase_sites_def2,
                disc_errors = disc,
                dist_n = nrow(test)
            ),
            switches1 = switches1
        )
    )
}







assess_ability_of_reads_to_be_confident <- function(
    hap1,
    hap2,
    hap3,
    sampleReads,
    maxDifferenceBetweenReads,
    minrp = 0.95,
    minmp = 0.95,
    method = "diploid"
) {
    ## 
    p <- calculate_eMatRead_t_vs_haplotypes(
        sampleReads,
        hap1 = hap1,
        hap2 = hap2,
        hap3 = hap3,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads ,
        rescale_eMatRead_t = FALSE,
        method = method
    )
    if (method == "diploid") {
        p1 <- p[1, ]
        p2 <- p[2, ]
        ## maxmimum
        mp <- p1
        mp[p2 > p1] <- p2[p2 > p1]
        ## ratio
        rp <- p1 / (p1 + p2)
        rp[is.na(rp)] <- 0.5 ## if BOTH 0 i.e. suuuper unlikely 
        rp[rp < 0.5] <- 1 - rp[rp < 0.5]
        conf <- (rp > minrp)  ## & (mp >
    } else {
        d <- colSums(p)
        p1 <- p[1, ] / d
        p2 <- p[2, ] / d
        p3 <- p[3, ] / d
        ## maxmimum
        mp <- p1
        mp[p2 > p1] <- p2[p2 > p1]
        mp[p3 > mp] <- p3[p3 > mp]
        ## look for largest ratio
        mp[is.na(mp)] <- 1/3
        conf <- (mp > minrp)  ## & (mp >
    }
    return(conf)
}


## for(i_gibbs_sample in 1:nGibbsSamples) {
##     hap_matrix <- super_out_hap_dosages[[i_gibbs_sample]]
##     read_label_matrix <- super_out_read_labels[[i_gibbs_sample]]
##     ##
##     gp_t <- rbind(
##     (1 - h1) * (1 - h2),
##     (1 - h1) * (h2) + (h1) * (1 - h2),
##     (h1) * (h2)
##     )
##     eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
##     fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##
##     thetaHat <- eij / 2
##     info <- 1 - (fij - eij**2) / 2 * thetaHat * (1-thetaHat)
## }



determine_best_read_label_so_far <- function(
    read_label_matrix_all,
    read_label_matrix_conf,
    nReads,
    nGibbsSamples,
    verbose,
    can_hap = 1
 ) {
    ##
    ## for default
    read_labels <- as.integer(read_label_matrix_all[, can_hap])
    flip_matrix <- array(FALSE, c(0, nGibbsSamples)) ##
    default_out <- list(
        read_labels = read_labels,
        flip_matrix = flip_matrix,
        read_label_matrix_all = read_label_matrix_all
    )
    ## how to assess "switch error rate"?
    ## do the reads partition themselves nicely
    a <- read_label_matrix_all
    a[!read_label_matrix_conf] <- NA
    a <- a[rowSums(is.na(a)) == 0, , drop = FALSE]
    if (nrow(a) < 10) {
        ## does not really matter
        if (verbose) {
            print_message("Insufficiently many confident reads for aggregating across runs")
        }
        return(default_out)
    }
    ## consider differences vs can_hap
    ## can_hap <- 1
    can <- a[, can_hap]
    a <- a - can
    s <- which(diff(rowSums(abs(a))) != 0) ## where they change, e.g. if a 39 here, it switches between 39 and 40
    ## probably some minimum number here
    if (length(s) == 0) {
        ## doesn't really need a message, the runs are the same at confident reads
        ##if (verbose) {
        ##    print_message("No differences between runs, Insufficiently many confident reads for aggregating across runs")
        ## }
        return(default_out)
    }
    ## so, e.g, want majority vote in each?
    ## what if I do a majority vote option
    ##
    s1 <- c(1, s + 1) #[-(length(s) + 1)]
    e1 <- c(s1[-1] - 1, nrow(a))
    ##
    flip_matrix <- array(FALSE, c(length(s1), nGibbsSamples))
    a2 <- a
    a2[] <- NA
    ## fuck it, just flip em
    for(i in 2:length(s1)) {
        s2 <- s1[i]
        prev <- a[e1[i - 1], ]
        cur <- a[s1[i], ]
        changed <- which(cur != 0)
        ## changed <- which(prev != cur)
        w <- s1[i]:nrow(a)
        if (length(changed) > 0) {
            if (length(changed) <= (nGibbsSamples / 2)) {
                ## trust canonical
                for(c1 in changed) {
                    reverted <- a[w, c1] + can[w]
                    reverted <- 3 - reverted
                    a[s1[i]:nrow(a), c1] <- reverted - can[w]
                }
                stopifnot(!(can_hap %in% changed))
            } else {
                ## here we change canonical, then any others identified
                changed <- which(cur == 0)
                stopifnot((can_hap %in% changed))
                for(c1 in changed) {
                    reverted <- a[w, c1] + can[w]
                    reverted <- 3 - reverted
                    a[w, c1] <- reverted - can[w]
                }
                reverted <- a[w, ] + can[w]
                ## now change canonical
                can[w] <- 3 - can[w]
                ## now re-set all of a
                a[w, ] <- reverted - can[w]
            }
        }
        flip_matrix[i, changed] <- TRUE
    }
    ## now - try using one of these and re-run and check phase
    ##
    ## now perform this flipping
    for(i_col in 1:ncol(flip_matrix)) {
        to_flip <- which(flip_matrix[, i_col])
        if (length(to_flip) > 0) {
            w <- s1[i]:nReads
            read_label_matrix_all[w, i_col] <- 3 - read_label_matrix_all[w, i_col]
        }
    }
    read_labels <- as.integer(read_label_matrix_all[, can_hap])
    return(
        list(
            read_labels = read_labels,
            flip_matrix = flip_matrix,
            read_label_matrix_all = read_label_matrix_all
        )
    )
}



determine_best_read_label_so_far_nipt <- function(
    read_label_matrix_all,
    read_label_matrix_conf,
    nReads,
    nGibbsSamples,
    verbose,
    can_hap = 1
) {
    ## pretend
    i_it <- 1 ## ideally do three of these iterations?
    ## can_hap <- 3
    ## method <- "nipt"
    read_label_matrix_all_ori <- read_label_matrix_all
    read_label_matrix_conf_ori <- read_label_matrix_conf
    ## put 3 into 2 but call them not confident
    if (i_it == 1) {
        read_label_matrix_conf[read_label_matrix_all == 3] <- FALSE
        read_label_matrix_all[read_label_matrix_all == 3] <- 2
    }
    ##
    out_best_labels <- determine_best_read_label_so_far(
        read_label_matrix_all = read_label_matrix_all,
        read_label_matrix_conf = read_label_matrix_conf,                
        nReads = nReads,
        nGibbsSamples = nGibbsSamples,
        verbose = verbose,
        can_hap = can_hap
    )
    read_label_matrix_all <- out_best_labels[["read_label_matrix_all"]]
    flip_matrix <- out_best_labels[["flip_matrix"]]
    if (i_it == 1) {
        ## now re-set the 3 values
        read_label_matrix_all[read_label_matrix_all_ori == 3] <- 3
    }
    if (verbose) {
        x <- out_best_labels$flip_matrix[, can_hap]
        if (length(x) > 0) {
            print_message(paste0("There are ", sum(x), " out of ", length(x), " regions that have been flipped by consensus"))
        }
    }
    read_labels <- read_label_matrix_all[, can_hap]
    list(read_labels = read_labels)
    ##  w <- rowSums(read_label_matrix_conf_ori) == 3
    ## cbind(read_label_matrix_all[w, ], NA,     read_label_matrix_all_ori[w, ])
}





if (1 == 0) {

    head(out$dosage_matrix)

    breaks2 <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5),
        1 - seq(0.1, 0.5, length.out = 5),
        1 - c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1)
    )))
    d <- hap1 + hap2
    d <- dosage_matrix[, "truth"]
    r2_by_freq(breaks2, af, truth_g, d, flip = FALSE)


    dosage_matrix[, "gibbs0"]
    m <- round(cbind(
        r2_by_freq(breaks, af, truth_g, dosage_matrix[, "gibbs0"])[, c("n", "nA")],
        apply(dosage_matrix, 2, function(dosage) {
            r2_by_freq(breaks, af, truth_g, dosage, flip = TRUE)[, "simple"]
        })
    ), 3)
    m



}

determine_a_set_of_truth_labels <- function(
    sampleReads,
    truth_hap1,
    truth_hap2,
    truth_hap3 = NULL,
    maxDifferenceBetweenReads
) {
    ## do not count where both sites missing
    w <- is.na(truth_hap1) & is.na(truth_hap2)
    truth_hap1[w] <- 0.5
    truth_hap2[w] <- 0.5
    eMatRead_truth_t <- calculate_eMatRead_t_vs_haplotypes(
        sampleReads,
        hap1 = truth_hap1,
        hap2 = truth_hap2,
        hap3 = truth_hap3,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        rescale_eMatRead_t = FALSE
    )
    p1 <- eMatRead_truth_t[1, ]
    p2 <- eMatRead_truth_t[2, ]
    weird <- (colSums(is.na(eMatRead_truth_t)) > 0) | (colSums(eMatRead_truth_t == Inf) > 0) | (colSums(eMatRead_truth_t == -Inf) > 0)
    ## if both 0, set to 1 for both
    if (sum(p1 == 0 | p2 == 0) > 0) {
        ##
        w <- (p1 == 0) & (p2 == 0) & (!weird)
        p1[w] <- 1
        p2[w] <- 1
        w <- (p1 == 0) & (p2 !=0) & (!weird)
        p1[w] <- 1e-100
        p2[w][p2[w] < p1[w]] <- 1e-90
        w <- (p1 != 0) & (p2 ==0) & (!weird)
        p2[w] <- 1e-100
        p1[w][p1[w] < p2[w]] <- 1e-90
    }
    ##
    nReads <- length(sampleReads)
    H <- as.integer(runif(nReads) < (p2 / (p1 + p2))) + 1
    H[is.na(H)] <- sample(1:2, sum(is.na(H)), replace = TRUE)
    truth_labels <- H ## at least, one of them!
    fold_diff <- p1 / p2
    fold_diff[weird] <- 1 ## hmm
    fold_diff[fold_diff > 1] <- (p2 / p1)[fold_diff > 1]
    uncertain_truth_labels <- (fold_diff > 0.5) ## |         (p1 < 1e-5 & p2 < 1e-5)
    return(
        list(
            truth_labels = truth_labels,
            uncertain_truth_labels = uncertain_truth_labels
        )
    )
}



impute_using_everything <- function(
    eMatDH_special_matrix_helper,
    eMatDH_special_matrix,
    use_eMatDH_special_symbols,
    H,
    sampleReads,
    rhb_t,
    nSNPs,
    full_alphaHat_t,
    full_betaHat_t,
    full_gamma_t,
    full_gammaSmall_t,
    full_gammaSmall_cols_to_get,
    full_transMatRate_t_H,
    distinctHapsB,
    distinctHapsIE,
    hapMatcher,
    hapMatcherR,
    use_hapMatcherR,
    eMatDH_special_grid_which,
    eMatDH_special_values_list,
    ref_error,
    make_plots,
    outplotprefix,
    have_truth_haplotypes,
    have_truth_genotypes,
    truth_haps,
    truth_gen,
    truth_labels,
    uncertain_truth_labels,
    L_grid,
    L,
    inRegion2,
    cM_grid,
    ancAlleleFreqAll,
    return_good_haps = FALSE,
    plot_description,
    Knew,
    previously_selected_haplotypes,
    sample_name,
    smooth_cm,
    regionStart,
    regionEnd,
    buffer,
    minGLValue,
    return_dosage = FALSE,
    return_betaHat_t = FALSE,
    return_gamma_t = FALSE,
    K_top_matches = 5,
    heuristic_match_thin = 0.01,
    suppressOutput = 1,
    method = "diploid",
    use_eigen = FALSE    
) {

    ##
    dosage <- numeric(nSNPs)
    nGrids <- ncol(rhb_t)
    n <- c(diploid = 2, nipt = 3)[method]
    
    if (use_hapMatcherR) {
        K <- nrow(hapMatcherR)
    } else {
        K <- nrow(hapMatcher)
    }
    
    dosage <- numeric(nSNPs)
    nGrids <- ncol(distinctHapsB)
    if (return_good_haps) {
        return_gammaSmall_t <- FALSE
        get_best_haps_from_thinned_sites <- TRUE
        best_haps_stuff_list <- as.list(1:sum(full_gammaSmall_cols_to_get >= 0))
        new_haps <- as.list(1:n)
    } else {
        new_haps <- NULL
        ##gammaSmall_t <-  array(0, c(1, 1))
        ##gammaSmall_cols_to_get <- array(0, 1)
        return_gammaSmall_t <- FALSE
        get_best_haps_from_thinned_sites <- FALSE
        best_haps_stuff_list <- list()
    }
    ##
    if (make_plots | return_gamma_t) {
        fbsoL <- as.list(1:3)
        fbsoL$list_of_gammas <- as.list(1:n)
        fbsoL$hapProbs_t <- array(NA, c(n, nSNPs))
        return_gamma_t <- TRUE ## slow but needed!
    } else {
        fbsoL <- NULL
    }
    dosage3 <- NULL

    for(i_hap in 1:n) {
        ##
        u <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
        bq <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
        w <- bq != 0
        bq <- bq[w]
        u <- u[w]
        if (length(u) == 0) {
            print_message(paste0("Read label assignment includes no reads for haplotype ", i_hap))
        }
        gl <- make_gl_from_u_bq(u, bq, nSNPs, minGLValue = minGLValue)
        use_eMatDH <- TRUE
        c <-  array(1, c(nGrids))  ## more useful for debugging
        if (use_eigen) {
            arma_alphaHat_t <- array(0, c(1, 1))
            eigen_alphaHat_t <- full_alphaHat_t
        } else {
            arma_alphaHat_t <- full_alphaHat_t
            eigen_alphaHat_t <- array(0, c(1, 1))
        }
        Rcpp_haploid_dosage_versus_refs(
            gl = gl,
            arma_alphaHat_t = arma_alphaHat_t,
            eigen_alphaHat_t = eigen_alphaHat_t,
            betaHat_t = full_betaHat_t,
            c = c,
            gamma_t = full_gamma_t,
            dosage = dosage,
            transMatRate_t = full_transMatRate_t_H,
            rhb_t = rhb_t,
            ref_error = ref_error,
            use_eMatDH = use_eMatDH,
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            use_eMatDH_special_symbols = use_eMatDH_special_symbols,
            suppressOutput = suppressOutput,
            return_dosage = return_dosage,
            return_betaHat_t = return_betaHat_t,
            return_gamma_t = return_gamma_t,
            return_gammaSmall_t = return_gammaSmall_t,
            gammaSmall_t = full_gammaSmall_t,
            gammaSmall_cols_to_get = full_gammaSmall_cols_to_get,
            get_best_haps_from_thinned_sites = get_best_haps_from_thinned_sites,
            best_haps_stuff_list = best_haps_stuff_list,
            K_top_matches = K_top_matches,
            always_normalize = FALSE,
            is_version_2 = !use_eigen, 
            normalize_emissions = TRUE,
            use_eigen = use_eigen
        )
        ## do some checks here
        if ((min(dosage) < -1e-5) | (1 + 1e-5) < max(dosage)) {
            print(range(dosage))
            stop("Dosage observed outside of range of 0 to 1 on forward-backward full iteration. Something has gone wrong. Please report this")
        }
        dosageNew <- numeric(nSNPs)
        dosageNew <- dosage[]
        if (i_hap == 1) { dosage1 <- dosageNew}
        if (i_hap == 2) { dosage2 <- dosageNew}
        if (i_hap == 3) { dosage3 <- dosageNew}        
        ##
        ## if (return_good_haps) {
        ##     new_haps <- everything_per_hap_prepare_haps(
        ##         full_gammaSmall_t = full_gammaSmall_t,
        ##         K_top_matches = K_top_matches,
        ##         new_haps = new_haps,
        ##         i_hap = i_hap
        ##     )
        ## }
        if (return_good_haps) {
            new_haps[[i_hap]] <- everything_per_hap_rejig_haps(best_haps_stuff_list)
        }
        if ( return_gamma_t | make_plots) {
            ## requires gamma
            fbsoL$hapProbs_t[i_hap, ] <- dosageNew
            gamma_tNew <- array(0, c(K, nGrids))
            gamma_tNew[] <- full_gamma_t
            if (i_hap == 1) { fbsoL$gammaMT_t <- gamma_tNew}
            if (i_hap == 2) { fbsoL$gammaMU_t <- gamma_tNew}
            if (i_hap == 3) { fbsoL$gammaP_t <- gamma_tNew}
            rm(gamma_tNew)
        }
    }

    if (return_good_haps) {
        ## select some new haplotypes if possible!
        ## try all the depth-up-to-5 ones first, then go on
        new_haps <- everything_select_good_haps(
            Knew = Knew,
            K_top_matches = K_top_matches,
            new_haps = new_haps,
            previously_selected_haplotypes = previously_selected_haplotypes,
            K = K
        )
    }

    to_return <- list(
        dosage1 = dosage1,
        dosage2 = dosage2,
        dosage3 = dosage3,
        new_haps = new_haps,
        fbsoL = fbsoL
    )

    if (make_plots) {
        
        plot_single_gamma_dosage(
            sampleReads = sampleReads,
            fbsoL = fbsoL,
            L_grid = L_grid,
            L = L,
            cM_grid = cM_grid,
            inRegion2 = inRegion2,
            ancAlleleFreqAll = ancAlleleFreqAll,
            outname = paste0(outplotprefix, plot_description, ".png"),
            haps = truth_haps,
            truth_gen = truth_gen,
            truth_labels = truth_labels,
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            uncertain_truth_labels = uncertain_truth_labels,
            sample_name = sample_name,
            smooth_cm = smooth_cm,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            method = method,
            new_haps = c(new_haps, previously_selected_haplotypes)
        )

    }

    if (1 == 0) {
        check_accuracy_from_all(cheat_all)
    }
    return(to_return)
}



## do much less this time
everything_per_hap_rejig_haps <- function(
    best_haps_stuff_list
) {
    new_haps <- lapply(best_haps_stuff_list, function(x) {
        top_matches <- x[["top_matches"]] + 1 ## turn from 0 to 1 based for R stuff
        top_matches_values <- x[["top_matches_values"]]
        return(top_matches[order(-top_matches_values)])
    })
    return(new_haps)
}


##
## todo, make this work within the forward backward bit, then can post-process in R more quickly
##
everything_per_hap_prepare_haps <- function(
    full_gammaSmall_t,
    K_top_matches,
    new_haps,
    i_hap
) {
    ## OK, so this is slow, the quantile bit argh
    ## how can I make this fast?
    ## can I write c++ to do two passes over this
    ## otherwise
    n <- ncol(full_gammaSmall_t)
    to_out <- as.list(1:n)
    for(i in 1:n) {
        ## for these samples, return them with this
        x <- full_gammaSmall_t[, i]
        y <- R_get_top_K_or_more_matches(x, K_top_matches)
        to_out[[i]] <- y
    }
    ## now choose
    new_haps[[i_hap]] <- to_out
    return(new_haps)
}

R_get_top_K_or_more_matches <- function(x, K_top_matches) {
    val <- -rcpp_nth_partial_sort(-x, as.integer(K_top_matches))[K_top_matches]
    y <- which(x >= val) ## 1-based, surprisingly slow
    y <- y[order(-x[y])]
    return(y)
}

R_get_top_K_or_more_matches_while_building_gamma <- function(
    alphaHat_t,
    betaHat_t_col,
    gamma_t_col,
    iGrid,
    K,
    K_top_matches
) {
    ##
    top_K_values <- numeric(K_top_matches)
    beats_value <- 0;
    num_top_K <- 1
    for(k in 0:(K - 1)) {
        gamma_t_col[k + 1] <- alphaHat_t[k + 1, iGrid + 1] * betaHat_t_col[k + 1];
        if (gamma_t_col[k + 1] == top_K_values[0 + 1]) {
            num_top_K <- num_top_K + 1;
        } else if (gamma_t_col[k + 1] > top_K_values[0 + 1]) {
            num_top_K <- num_top_K + 1;
            beats_value <- 0;
            for(j in 0:(K_top_matches - 1)) {
                if (gamma_t_col[k + 1] > top_K_values[j + 1]) {
                    beats_value <- j
                }
            }
            if (beats_value > 0) {
                for(i in 0:(beats_value)) {
                    top_K_values[i + 1] = top_K_values[i + 1 + 1]
                }
            }
            top_K_values[beats_value + 1] <- gamma_t_col[k + 1]
        }
    }
    ##
    top_matches <- integer(num_top_K)
    top_matches_values <- numeric(num_top_K)
    count <- 0
    for(k in 0:(K - 1)) {
        if (gamma_t_col[k + 1] >= top_K_values[0 + 1]) {
            top_matches[count + 1] <- k + 1
            top_matches_values[count + 1] <- gamma_t_col[k + 1]
            count <- count + 1;
        }
    }
    ## note, can just leave for R in this way
    return(
        list(
            top_matches = top_matches[1:count] - 1, ## make 0-based
            top_matches_values = top_matches_values[1:count]
        )
    )
    ##top_matches[order(top_matches_values , decreasing = TRUE)])
}



everything_select_good_haps <- function(
    Knew,
    K_top_matches,
    new_haps,
    previously_selected_haplotypes,
    K
) {
    i <- 1
    to_keep <- NULL
    left <- Knew
    done <- FALSE
    while(!done) {
        if (i <= K_top_matches) {
            new <- unique(unlist(sapply(new_haps, function(x) lapply(x, function(y) y[i]))))
        } else {
            ## here you just take all because exhausted
            new <- unique(unlist(new_haps))
            done <- TRUE
        }
        new <- setdiff(new, previously_selected_haplotypes)
        new <- setdiff(new, to_keep)
        ## if (length(new) == 0) {
        ##     ## should not happen, have run out! then add some at random
        ##     print("AA")
        ##     toadd <- Knew - length(to_keep)
        ##     done <- TRUE
        ##     to_keep <- c(to_keep, sample(setdiff(setdiff(1:K, previously_selected_haplotypes), to_keep), toadd))
        if (length(new) < (Knew - length(to_keep))) {
            to_keep <- c(to_keep, new)
            i <- i + 1
        } else {
            toadd <- Knew - length(to_keep)
            to_keep <- c(to_keep, new[sample(1:length(new), toadd)])
            done <- TRUE
        }
    }
    ##     new_haps <- to_keep
    if (length(to_keep) < Knew) {
        ## if not enough, add some at random
        to_keep <- c(to_keep, sample(setdiff(1:K, c(to_keep, previously_selected_haplotypes)), Knew - length(to_keep)))
    }
    if (length(to_keep) != Knew) {
        print(paste0("length(to_keep) = ", length(to_keep)))
        print(paste0("Knew = ", Knew))
        save(Knew, K_top_matches, new_haps, previously_selected_haplotypes, K, file = "~/temp.RData")
        stop("Have returned too many haps, see debug file")
    }
    return(to_keep)
}


impute_one_sample <- function(
    eMatDH_special_matrix_helper = array(0, c(1, 1)),
    eMatDH_special_matrix = array(0,c(1, 1)),
    use_eMatDH_special_symbols = FALSE,
    distinctHapsB = array(0L, c(1, 1)),
    distinctHapsIE = array(0, c(1, 1)),
    hapMatcher = array(0L, c(1, 1)),
    hapMatcherR = array(as.raw(0), c(1, 1)),
    use_hapMatcherR = FALSE,
    rhb_t = array(0L, c(1, 1)),
    ref_error = 0.001,
    nSNPs,
    sampleReads,
    ff,
    small_eHapsCurrent_tc = array(0, c(1, 1, 1)),
    small_transMatRate_tc_H,
    alphaHat_t1,
    betaHat_t1,
    eMatGrid_t1,
    alphaHat_t2,
    betaHat_t2,
    eMatGrid_t2,
    alphaHat_t3,
    betaHat_t3,
    eMatGrid_t3,
    gammaMT_t_local,
    gammaMU_t_local,
    gammaP_t_local,
    small_alphaMatCurrent_tc,
    small_priorCurrent_m,
    smooth_cm,
    which_haps_to_use,
    n_gibbs_starts = 1,
    small_ref_panel_gibbs_iterations,
    n_gibbs_sample_its = 1,
    double_list_of_starting_read_labels,
    small_ref_panel_block_gibbs_iterations,
    perform_block_gibbs = TRUE,
    make_plots,
    maxDifferenceBetweenReads,
    wif0,
    grid_has_read,
    verbose = FALSE,
    shuffle_bin_radius,
    outplotprefix,
    plot_description,
    ancAlleleFreqAll,
    grid,
    L_grid,
    L,
    inRegion2,
    cM_grid = NULL,
    have_truth_haplotypes,
    truth_haps,
    have_truth_genotypes,
    truth_gen,
    truth_labels,
    sample_name,
    regionStart,
    regionEnd,
    buffer,
    uncertain_truth_labels,
    small_ref_panel_skip_equally_likely_reads = FALSE,
    small_ref_panel_equally_likely_reads_update_iterations = c(1,2,3,6,9,15),
    return_p_store = FALSE,
    return_p1 = FALSE,
    return_extra = FALSE,
    return_genProbs = TRUE,
    return_hapProbs = TRUE,
    return_gamma = FALSE,
    return_gibbs_block_output = FALSE,
    return_advanced_gibbs_block_output = FALSE,
    gibbs_initialize_iteratively = FALSE,
    gibbs_initialize_at_first_read = FALSE,
    maxEmissionMatrixDifference = 1e100,
    rescale_eMatRead_t = TRUE,
    rescale_eMatGrid_t = FALSE,
    Jmax = 10000,
    suppressOutput = 1,
    use_smooth_cm_in_block_gibbs = TRUE,
    block_gibbs_quantile_prob = 0.95,
    make_plots_block_gibbs = FALSE,
    use_small_eHapsCurrent_tc = TRUE,
    method = "diploid",
    use_provided_small_eHapsCurrent_tc = FALSE,
    use_sample_is_diploid = FALSE,
    i_it = NA,
    i_gibbs_sample = NA,
    shard_check_every_pair = FALSE,
    calculate_gamma_on_the_fly = FALSE,
    eMatRead_t = NULL,
    make_eMatRead_t_rare_common = FALSE,
    common_snp_index = integer(1),
    snp_is_common = logical(1),
    rare_per_hap_info = vector("list", 1),
    rare_per_snp_info = vector("list", 1),
    disable_read_category_usage = TRUE
) {


    ## file <- paste0("/well/davies/users/dcc832/werAwerBwerC.", i_it, ".", nSNPs, ".RData")
    ## print(paste0("saving to ", file))
    ## save(
    ## eMatDH_special_matrix_helper,
    ## eMatDH_special_matrix,
    ## use_eMatDH_special_symbols,
    ## distinctHapsB ,
    ## distinctHapsIE,
    ## hapMatcher ,
    ## hapMatcherR,
    ## use_hapMatcherR ,
    ## rhb_t,
    ## ref_error,
    ## nSNPs,
    ## sampleReads,
    ## ff,
    ## small_eHapsCurrent_tc,
    ## small_transMatRate_tc_H,
    ## alphaHat_t1,
    ## betaHat_t1,
    ## eMatGrid_t1,
    ## alphaHat_t2,
    ## betaHat_t2,
    ## eMatGrid_t2,
    ## alphaHat_t3,
    ## betaHat_t3,
    ## eMatGrid_t3,
    ## gammaMT_t_local,
    ## gammaMU_t_local,
    ## gammaP_t_local,
    ## small_alphaMatCurrent_tc,
    ## small_priorCurrent_m,
    ## smooth_cm,
    ## which_haps_to_use,
    ## n_gibbs_starts,
    ## small_ref_panel_gibbs_iterations,
    ## n_gibbs_sample_its,
    ## double_list_of_starting_read_labels,
    ## small_ref_panel_block_gibbs_iterations,
    ## perform_block_gibbs,
    ## make_plots,
    ## maxDifferenceBetweenReads,
    ## wif0,
    ## grid_has_read,
    ## verbose,
    ## shuffle_bin_radius,
    ## outplotprefix,
    ## plot_description,
    ## ancAlleleFreqAll,
    ## grid,
    ## L_grid,
    ## L,
    ## inRegion2,
    ## cM_grid,
    ## have_truth_haplotypes,
    ## truth_haps,
    ## have_truth_genotypes,
    ## truth_gen,
    ## truth_labels,
    ## sample_name,
    ## regionStart,
    ## regionEnd,
    ## buffer,
    ## uncertain_truth_labels,
    ## small_ref_panel_skip_equally_likely_reads ,
    ## small_ref_panel_equally_likely_reads_update_iterations,
    ## return_p_store,
    ## return_p1,
    ## return_extra,
    ## return_genProbs,
    ## return_hapProbs,
    ## return_gamma,
    ## return_gibbs_block_output,
    ## return_advanced_gibbs_block_output,
    ## gibbs_initialize_iteratively,
    ## gibbs_initialize_at_first_read,
    ## maxEmissionMatrixDifference,
    ## rescale_eMatRead_t,
    ## rescale_eMatGrid_t,
    ## Jmax,
    ## suppressOutput,
    ## use_smooth_cm_in_block_gibbs,
    ## block_gibbs_quantile_prob,
    ## make_plots_block_gibbs,
    ## use_small_eHapsCurrent_tc,
    ## method,
    ## use_provided_small_eHapsCurrent_tc,
    ## use_sample_is_diploid,
    ## i_it,
    ## i_gibbs_sample,
    ## shard_check_every_pair,
    ## calculate_gamma_on_the_fly,
    ## eMatRead_t,
    ## make_eMatRead_t_rare_common,
    ## common_snp_index,
    ## snp_is_common,
    ## rare_per_hap_info,
    ## rare_per_snp_info,
    ## disable_read_category_usage,
    ## compress = FALSE,
    ## file = file)
    ## ## stop("WER")
    ## ## print("SAVING")
    ## ## if (exit_after) {
    ## ## stop("WER")
    ## ## }

    ##
    K <- length(which_haps_to_use)
    S <- 1
    ## print(paste0("start = ", Sys.time()))
    if (use_small_eHapsCurrent_tc & !use_provided_small_eHapsCurrent_tc) {
        inflate_fhb_t_in_place(
            rhb_t,
            small_eHapsCurrent_tc,
            haps_to_get = which_haps_to_use - 1,
            nSNPs = nSNPs,
            ref_error = ref_error
        )
    }
    if (make_plots_block_gibbs) {
        return_gibbs_block_output <- TRUE
        return_advanced_gibbs_block_output <- TRUE
    }
    ##
    if (use_sample_is_diploid && ff == 0) {
        sample_is_diploid <- TRUE
    }else {
        sample_is_diploid <- FALSE
    }
    ## I dunno, could clearly be fixed
    if (is.null(eMatRead_t)) {
        pass_in_eMatRead_t <- FALSE
        eMatRead_t <- array(0, c(1, 1))
    } else {
        pass_in_eMatRead_t <- TRUE
    }

    if (ff == 0) {
        do_shard_block_gibbs <- TRUE
    } else {
        do_shard_block_gibbs <- FALSE
    }
    
    param_list <- list(
        return_alpha = FALSE,
        return_extra = return_extra,
        return_genProbs = return_genProbs,
        return_gamma = as.logical(return_gamma | make_plots),
        return_hapProbs = as.logical(return_hapProbs | make_plots),
        return_p_store = return_p_store,
        return_p1 = return_p1,
        return_gibbs_block_output = return_gibbs_block_output,
        return_advanced_gibbs_block_output = return_advanced_gibbs_block_output,
        use_starting_read_labels = TRUE,
        verbose = verbose,
        run_fb_subset = FALSE,
        haploid_gibbs_equal_weighting = TRUE,
        gibbs_initialize_iteratively = gibbs_initialize_iteratively,
        gibbs_initialize_at_first_read = gibbs_initialize_at_first_read,
        use_smooth_cm_in_block_gibbs = use_smooth_cm_in_block_gibbs,
        use_small_eHapsCurrent_tc = use_small_eHapsCurrent_tc,
        sample_is_diploid = sample_is_diploid,
        update_in_place = FALSE,
        do_shard_block_gibbs = do_shard_block_gibbs,
        shard_check_every_pair = shard_check_every_pair,
        force_reset_read_category_zero = FALSE,
        disable_read_category_usage = disable_read_category_usage,
        calculate_gamma_on_the_fly = calculate_gamma_on_the_fly,
        rescale_eMatRead_t = rescale_eMatRead_t,
        pass_in_eMatRead_t = pass_in_eMatRead_t,
        make_eMatRead_t_rare_common = make_eMatRead_t_rare_common,
        pass_in_alphaBeta = TRUE,
        update_hapSum = FALSE,
        record_read_set = TRUE,
        perform_block_gibbs = perform_block_gibbs,
        use_eMatDH_special_symbols = use_eMatDH_special_symbols
    )
    ## use_provided_small_eHapsCurrent_tc = use_provided_small_eHapsCurrent_tc
    ## this should catch hopefully rare underflow problems and re-run the samples
    done_imputing <- FALSE
    n_imputing <- 0
    n_gibbs_burn_in_its <- small_ref_panel_gibbs_iterations
    n_gibbs_full_its <- n_gibbs_burn_in_its + n_gibbs_sample_its
    skip_read_iteration <- rep(FALSE, n_gibbs_full_its)
    if (small_ref_panel_skip_equally_likely_reads) {
        skip_read_iteration <- rep(TRUE, n_gibbs_full_its)
        skip_read_iteration[small_ref_panel_equally_likely_reads_update_iterations] <- FALSE
    }
    while(!done_imputing) {

        out <- rcpp_forwardBackwardGibbsNIPT(
            sampleReads = sampleReads,
            eMatRead_t = eMatRead_t,
            priorCurrent_m = small_priorCurrent_m,
            alphaMatCurrent_tc = small_alphaMatCurrent_tc,
            eHapsCurrent_tc = small_eHapsCurrent_tc,
            transMatRate_tc_H = small_transMatRate_tc_H,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            use_hapMatcherR = use_hapMatcherR,
            distinctHapsB =distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            eMatDH_special_matrix = eMatDH_special_matrix,
            rhb_t = rhb_t,
            ref_error = ref_error,
            which_haps_to_use = which_haps_to_use,
            ff = ff,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax_local = Jmax,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            run_fb_grid_offset = FALSE,
            blocks_for_output = array(0, c(1, 1)),
            grid = grid,
            skip_read_iteration = skip_read_iteration,
            alphaHat_t1 = alphaHat_t1,
            alphaHat_t2 = alphaHat_t2,
            alphaHat_t3 = alphaHat_t3,
            betaHat_t1 = betaHat_t1,
            betaHat_t2 = betaHat_t2,
            betaHat_t3 = betaHat_t3,
            eMatGrid_t1 = eMatGrid_t1,
            eMatGrid_t2 = eMatGrid_t2,
            eMatGrid_t3 = eMatGrid_t3,
            gammaMT_t_local = gammaMT_t_local,
            gammaMU_t_local = gammaMU_t_local,
            gammaP_t_local = gammaP_t_local,
            hapSum_tc = array(0, c(1, 1, 1)),
            snp_start_1_based = -1,
            snp_end_1_based = -1,
            generate_fb_snp_offsets = FALSE,
            suppressOutput = suppressOutput, ## 
            n_gibbs_burn_in_its = n_gibbs_burn_in_its,
            n_gibbs_sample_its = n_gibbs_sample_its,
            n_gibbs_starts = n_gibbs_starts,
            double_list_of_starting_read_labels = double_list_of_starting_read_labels,
            prev_list_of_alphaBetaBlocks = as.list(c(1, 2)),
            i_snp_block_for_alpha_beta = -1,
            do_block_resampling = FALSE, ## turn off for now
            seed_vector = 0,
            class_sum_cutoff = 0.06, ## what is this
            wif0 = wif0,
            grid_has_read = grid_has_read,
            shuffle_bin_radius = shuffle_bin_radius,
            L_grid = L_grid,
            block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
            smooth_cm = smooth_cm,
            param_list = param_list,
            block_gibbs_quantile_prob = block_gibbs_quantile_prob,
            artificial_relabel = -1,
            common_snp_index = common_snp_index,
            snp_is_common = snp_is_common,
            rare_per_hap_info = rare_per_hap_info,
            rare_per_snp_info = rare_per_snp_info
        )

## ################################################################ werwerwer 

##         print(paste0("Saving - ", paste0("/data/smew1/rdavies/temp.", i_gibbs_sample, ".", i_it, ".RData")))
##         save(out, file = paste0("/data/smew1/rdavies/temp.", i_gibbs_sample, ".", i_it, ".RData"))
##         print(paste0("Done saving - ", i_it))        
##         print(names(out))
##         print(table(out$H))
##         print(table(out$H_class))
##         print(table(out$H_class, out$H))
##         print(table(            out$double_list_of_ending_read_labels[[1]][[1]], out$H))
        
##         print(out$per_it_likelihoods[, c("i_it", "p_O_given_H_L", "p_H_given_L", "p_H_class_given_L")])
        
##         x1 <- out$gibbs_block_output_list[[1]]$gibbs_block_output$block_results
##         x2 <- out$gibbs_block_output_list[[2]]$gibbs_block_output$block_results
##         x3 <- out$gibbs_block_output_list[[3]]$gibbs_block_output$block_results
##         print(round(rbind(
##             x1[c(1, nrow(x1) - 1), ],
##             x2[c(1, nrow(x2) - 1), ],
##             x3[c(1, nrow(x3) - 1), ]
##         ), 3))
## ################################################################ werwerwer 
        

        if (out[["underflow_problem"]]) {
            new_maxDifferenceBetweenReads <- max(1, maxDifferenceBetweenReads / 10)
            print_message(paste0("Underflow problem observed for sample ", sample_name, ". Re-setting maxDifferenceBetweenReads from ", maxDifferenceBetweenReads, " (log10 of ", log10(maxDifferenceBetweenReads), ") to ", new_maxDifferenceBetweenReads, " (log10 of ", log10(new_maxDifferenceBetweenReads), "). If this problem persists please reset maxDifferenceBetweenReads or downsampleToCov to lower coverage and/or reduce the impact of individual reads"))
            maxDifferenceBetweenReads <- new_maxDifferenceBetweenReads
            n_imputing <- n_imputing + 1
            if (n_imputing > 10) {
                stop(paste0("There were 5 consecutive underflow problems for sample ", sample_name, ". This is likely due to high coverage (at least locally) and/or reads that are very long. Please try setting downsampleToCov to a lower value and/or reduce maxDifferenceBetweenReads to decrease the probability of this happening again"))
            }
        } else {
            done_imputing <- TRUE
        }
    }
    ##
    genProbs_t <- out$genProbsM_t
    hapProbs_t <- out$happrobs_t
    dosage <- genProbs_t[2, ] + 2 * genProbs_t[3, ]
    out$dosage <- dosage
    ##
    if (make_plots) {
        plot_single_gamma_dosage(
            sampleReads = sampleReads,
            fbsoL = out,
            L_grid = L_grid,
            L = L,
            cM_grid = cM_grid,
            inRegion2 = inRegion2,
            ancAlleleFreqAll = ancAlleleFreqAll,
            outname = paste0(outplotprefix, plot_description, ".png"),
            method = method,
            haps = truth_haps,
            truth_labels = truth_labels,
            have_truth_haplotypes = have_truth_haplotypes,
            have_truth_genotypes = have_truth_genotypes,
            truth_gen = truth_gen,
            uncertain_truth_labels = uncertain_truth_labels,
            sample_name = sample_name,
            smooth_cm = smooth_cm,
            regionStart = regionStart,
            regionEnd = regionEnd,
            buffer = buffer,
            which_haps_to_use = which_haps_to_use
        )
        ## would like a plot here
    }
    if (make_plots_block_gibbs) {
        nGrids <- ncol(alphaHat_t1)
        for(n_block_it_to_plot in 1:length(small_ref_panel_block_gibbs_iterations)) {
            outname <- paste0(outplotprefix, "block", plot_description, ".n", n_block_it_to_plot, ".png")
            n_block_it_to_plot = n_block_it_to_plot
            plot_attempt_to_reblock_snps(
                out = out,
                nGrids = nGrids,
                block_gibbs_iterations = small_ref_panel_block_gibbs_iterations,
                outname = outname,
                L_grid = L_grid,
                L = L,
                uncertain_truth_labels = uncertain_truth_labels,
                truth_labels = truth_labels,
                have_truth_haplotypes = have_truth_haplotypes,
                sampleReads = sampleReads,
                n_block_it_to_plot = n_block_it_to_plot,
                wif0 = wif0,
                grid = grid,
                method = method,
                ff = ff
            )
        }
    }
    return(out)
}




get_alleleCount <- function(sampleReads , nSNPs) {
    alleleCount <- array(0,c(nSNPs, 3))
    a=unlist(sapply(sampleReads,function(x) x[[3]]))
    b=unlist(sapply(sampleReads,function(x) x[[4]]))
    bqProbs=convertScaledBQtoProbs(matrix(a,ncol=1))
    ## y is numeric, z = integer counts (can have 0),
    c1 <- increment2N(
        y = as.numeric(bqProbs[,1]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    c2 <- increment2N(
        y = as.numeric(bqProbs[,2]),
        z = as.numeric(b),
        yT = as.integer(nrow(bqProbs)),
        xT = as.integer(nSNPs - 1)
    )
    alleleCount[,1]=alleleCount[,1]+c2 # fixed nov 6 2015 - was backward
    alleleCount[,2]=alleleCount[,2]+c1+c2
    return(alleleCount)
}


#' @export
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE) {
    if (flip) {
        w <- af > 0.5
        af[w] <- 1 - af[w]
        truthG[w] <- 2 - truthG[w]
        testDS[w] <- 2 - testDS[w]
    }
    if (!is.null(which_snps)) {
        af <- af[which_snps]
        truthG <- truthG[which_snps]
        testDS <- testDS[which_snps]
    }
    x <- cut(af, breaks = breaks)
    cors_per_af <- tapply(1:length(x), x, function(w) {
        c(
            n = length(w),
            nA = sum(truthG[w], na.rm = TRUE),
            simple = cor(truthG[w], testDS[w], use = 'pairwise.complete') ** 2,
            norm = cor(truthG[w] - 2 * af[w], testDS[w] - 2 * af[w], use = 'pairwise.complete') ** 2
        )
    })
    cors_per_af <- t(sapply(cors_per_af[!sapply(cors_per_af, is.null)], I))
    return(cors_per_af)
}


rsync_to_rescomp <- function(ext = "png") {
    analysis_date <- basename(getwd()) ## probably
    system(paste0("rsync --progress -av plots/*", ext, " rescompNew:/gpfs3/well/davies/users/dcc832/single_imp_test/", analysis_date, "/plots/ "))
}



calc_metric3 <- function(
    eMatRead_t,
    rare_thresh = 0.99,
    hap_thresh = 0.10,
    read_length_thresh = 0.80
) {
    ## among longest 20% of reads
    ## count as rare those occuring with probability of at least that of 99th percentile
    ## but found in fewer than 5% of samples
    x2 <- apply(eMatRead_t, 2, function(x) {
        rare <- x >= quantile(x, probs = rare_thresh, na.rm = TRUE)
        return(sum(rare))
    })
    is_rare <- x2 <= hap_thresh * nrow(eMatRead_t)
    is_long <- nr >= quantile(nr, probs = read_length_thresh)
    which_reads_to_use <- which(is_rare & is_long)
    ## count the occurences of these things
    a <- unlist(lapply(which_reads_to_use, function(iRead) {
        x <- eMatRead_t[, iRead]
        rare <- x >= quantile(x, probs = rare_thresh)
        return(which(rare))
    }))
    x <- increment2N(
        as.integer(length(a)),
        as.integer(nrow(eMatRead_t)),
        as.integer(rep(1, length(a))),
        as.integer(a)
    )[-1] ## first entry is 0-match
    return(
        list(
            x = x,
            which_reads_to_use = which_reads_to_use
        )
    )
    ## alt, among these long, rare reads
    ## look for reads ful-filling a need
}


calculate_eMatRead_t_usig_rhb_t <- function(
    sampleReads,
    rhb_t,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- nrow(rhb_t)
    eMatRead_t_full  <- array(1, c(K, nReads))
    s <- 0
    ## slow but whatevs, will fix
    n <- 5000
    eMatRead_t_full  <- array(1, c(K, nReads))
    for(i in 1:ceiling(K / n)) {
        print(paste0(i, " / ", ceiling(K / n)))
        s <- n * (i - 1) + 1
        e <- min(n * i, K)
        KL <- e - s + 1
        ehc <- array(0, c(KL, nSNPs, 1))
        ## do it in parts
        ehc[, , 1] <- inflate_fhb_t(
            rhb_t,
            haps_to_get = s:e - 1,
            nSNPs
        )
        eMatRead_t_partial  <- array(1, c(KL, nReads))
        gc(reset = TRUE); gc(reset = TRUE);
        ##
        rcpp_make_eMatRead_t(
            eMatRead_t = eMatRead_t_partial,
            sampleReads = sampleReads,
            eHapsCurrent_tc = ehc,
            s = 0,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            Jmax = 10000,
            eMatHapOri_t = array(0, c(1, 1)), ## ugh
            pRgivenH1 = array(0),
            pRgivenH2 = array(0),
            run_pseudo_haploid = FALSE,
            prev = 0,
            suppressOutput = 1,
            prev_section = "text",
            next_section = "text",
            rescale_eMatRead_t = rescale_eMatRead_t
        )
        eMatRead_t_full[s:e, ] <- eMatRead_t_partial
        if ((i %% 10) == 0) {
            gc(reset = TRUE);    gc(reset = TRUE);
        }
    }
    gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t_full)
}



calculate_eMatRead_t_some_haplotypes <- function(
    sampleReads,
    rhb_t,
    which_haps_to_use,
    maxDifferenceBetweenReads,
    rescale_eMatRead_t = TRUE
) {
    nReads <- length(sampleReads)
    K <- length(which_haps_to_use) ## 1-based
    eMatRead_t  <- array(1, c(K, nReads))
    s <- 0
    ## slow but whatevs, will fix
    n <- 1000
    eMatRead_t  <- array(1, c(K, nReads))
    ehc <- array(0, c(K, nSNPs, 1))
    ## do it in parts
    ehc[, , 1] <- inflate_fhb_t(
        rhb_t,
        haps_to_get = which_haps_to_use - 1,
        nSNPs
    )
    gc(reset = TRUE); gc(reset = TRUE);
    ##
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = ehc,
        s = 0,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)), ## ugh
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        run_pseudo_haploid = FALSE,
        prev = 0,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text",
        rescale_eMatRead_t = rescale_eMatRead_t
    )
    gc(reset = TRUE);    gc(reset = TRUE);
    gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t)
}



calculate_eMatRead_t_vs_haplotypes <- function(
    sampleReads,
    hap1,
    hap2,
    hap3 = NULL,
    maxDifferenceBetweenReads,
    rescale_eMatRead_t = TRUE,
    method = "diploid"
) {
    nReads <- length(sampleReads)
    K <- c(diploid = 2, nipt = 3)[method]
    eMatRead_t <- array(1, c(K, nReads))
    s <- 0
    ## expand whole thing - yuck - fix this later!
    nSNPs <- length(hap1)
    ehc <- array(0, c(K, nSNPs, 1))    
    ehc[1, , 1] <- hap1
    ehc[2, , 1] <- hap2
    if (method == "nipt") {
        ehc[3, , 1] <- hap3        
    }
    ## shouldn't be necessary
    ## if (nSNPs > 10000) {
    ##     gc(reset = TRUE); gc(reset = TRUE);
    ## }
    ## 
    rcpp_make_eMatRead_t(
        eMatRead_t = eMatRead_t,
        sampleReads = sampleReads,
        eHapsCurrent_tc = ehc,
        s = s,
        maxDifferenceBetweenReads = maxDifferenceBetweenReads,
        Jmax = 1000,
        eMatHapOri_t = array(0, c(1, 1)), ## ugh
        pRgivenH1 = array(0),
        pRgivenH2 = array(0),
        run_pseudo_haploid = FALSE,
        prev = 0,
        suppressOutput = 1,
        prev_section = "text",
        next_section = "text",
        rescale_eMatRead_t = rescale_eMatRead_t
    )
    rm(ehc)
    ## gc(reset = TRUE);    gc(reset = TRUE);
    return(eMatRead_t)
}






##
## check accuracy quickly
##
check_accuracy_from_all <- function(all) {
    dosage1 <- all$dosage1
    dosage2 <- all$dosage2
    r2s <- c(
        cor(na12878_hap1, dosage1) ** 2,
        cor(na12878_hap2, dosage1) ** 2,
        cor(na12878_hap1, dosage2) ** 2,
        cor(na12878_hap2, dosage2) ** 2
    )
    testDS <- dosage1 + dosage2
    breaks <- sort(unique(c(
        seq(0, 1, length.out = 11)
    )))
    r2sb <- r2_by_freq(breaks, af, truth_g, testDS)
    r2sc <- r2_by_freq(c(0, 1), af, truth_g, testDS)
    afX <- af
    afX[afX > 0.5] <- (1 - af)[afX > 0.5]
    breaks <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5)
    )))
    r2sm <- r2_by_freq(breaks, afX, truth_g, testDS)
    r2sm2 <- r2_by_freq(breaks, af, truth_g, testDS, flip = TRUE)
    breaks2 <- sort(unique(c(
        c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
        seq(0.1, 0.5, length.out = 5),
        1 - seq(0.1, 0.5, length.out = 5),
        1 - c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
        1 - c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
        1 - c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1)
    )))
    afX <- af
    r2sn <- r2_by_freq(breaks2, afX, truth_g, testDS)
    to_return <- list(
        r2s = r2s,
        r2sb = r2sb,
        r2sc = r2sc,
        r2sm = r2sm,
        r2sm2 = r2sm2,
        r2sn = r2sn
    )
    to_return
}






if (1 == 0) {

        ## check them inidividually
    ## check them in aggregate
    f <- function(output) {
        d1 <-
            super_out[[1]][[output]][["out"]][[1]]$dosage +
            super_out[[1]][[output]][["out"]][[2]]$dosage
        d2 <-
            super_out[[2]][[output]][["out"]][[1]]$dosage +
            super_out[[2]][[output]][["out"]][[2]]$dosage
        d3 <-
            super_out[[3]][[output]][["out"]][[1]]$dosage +
            super_out[[3]][[output]][["out"]][[2]]$dosage
        d4 <-
            super_out[[4]][[output]][["out"]][[1]]$dosage +
            super_out[[4]][[output]][["out"]][[2]]$dosage
        breaks <- sort(unique(c(
            seq(0, 0.01, length.out = 6),
            seq(0, 0.05, length.out = 6),
            seq(0, 1, length.out = 11)
        )))
        ## just 1
        cbind(
            t(r2_by_freq(breaks, af, truth_g, d1))[, c("n","simple")],
            t(r2_by_freq(breaks, af, truth_g, (d1 + d2) / 2))[, "simple"],
            t(r2_by_freq(breaks, af, truth_g, (d1 + d2 + d3 + d4) / 4))[, "simple"]
        )
    }

    100 * f("cheat_all")
    100 * f("impute_all1")
    100 * f("impute_all2")
    100 * f("impute_all3")

    100 * f("cheat_all")
    100 * f("impute_all1")
    100 * f("impute_all2")
    100 * f("impute_all3")

    ## however - how well does this do for UNLINKED reads
    ## to some extent, how much do I care? (somewhat?)
    ## can try the above, but initialize differently somehow?
    ## 400 haplotypes at random?
    ## hopefully OK?



}


    get_dosages_from_fbsoL <- function(fbsoL) {
        return(1 * fbsoL$genProbsM_t[2, ] + 2 * fbsoL$genProbsM_t[3, ])
    }






estimate_bq <- function(truth_labels, sampleReads, truth_haps) {
    ## do my own bqsr versus truth
    H <- truth_labels
    i_hap <- 1
    u1 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq1 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    i_hap <- 2
    u2 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[4]])) + 1
    bq2 <- unlist(sapply(sampleReads[H == i_hap], function(x) x[[3]]))
    ##
    bqs <- c(bq1, bq2)
    vals <- c(truth_haps[u1, 1], truth_haps[u2, 2])
    ##
    cut_bqs <- cut(bqs, breaks = c(-Inf, -30, -20, -1, 0, 19, 29, Inf))
    results <- t(sapply(c("(-Inf,-30]", "(-30,-20]", "(-20,-1]", "(0,19]", "(19,29]", "(29, Inf]"), function(bqlevel) {
        w <- cut_bqs == bqlevel
        w0 <- cut_bqs == bqlevel & vals == 0
        w1 <- cut_bqs == bqlevel & vals == 1
        phred_stated <- round(exp(1 / sum(w) * sum(sum(log(abs(bqs[w]))))), 1)
        return(c(sum(w0), sum(w1), phred_stated))
    }))
    colnames(results)[1:3] <- c("N_ref", "N_alt", "phred_stated")
    rownames(results) <- c("(-Inf,-30]", "[-29,-20]", "[-19, 1]", "[1, 19]", "[20, 29]", "[30, Inf)")
    results <- cbind(results, accuracy = c(results[1:3, 1] / rowSums(results[1:3, 1:2]), results[4:6, 2] / rowSums(results[4:6, 1:2])))
    results <- cbind(results, phred_obs = -round(log10(1 - results[, "accuracy"]) * 10, 1))
    results <- results[, colnames(results) != "accuracy"]
    return(results)
}




## force phased haplotypes to agree with genotype posterior
recast_haps <- function(hd1, hd2, gp, force_round = FALSE) {
    if (force_round) {
        hd1 <- round(hd1)
        hd2 <- round(hd2)
    }
    gt1 <- round(hd1) + round(hd2)
    max_val <- gp[, 1]
    gt3 <- rep(0, nrow(gp))
    for(i in 2:3) {
        w <- gp[, i] > max_val
        gt3[w] <- i - 1
        max_val[w] <- gp[w, i]
    }
    ##
    to_change <- which(gt3 != gt1)
    ## um, round all of them? 
    ## easy one
    hd1[to_change][gt3[to_change] == 0] <- 0
    hd2[to_change][gt3[to_change] == 0] <- 0
    hd1[to_change][gt3[to_change] == 2] <- 1
    hd2[to_change][gt3[to_change] == 2] <- 1
    ## hard one
    a1 <- hd1[to_change][gt3[to_change] == 1]
    a2 <- hd2[to_change][gt3[to_change] == 1]
    hd1[to_change][gt3[to_change] == 1][a1 > a2] <- 1
    hd2[to_change][gt3[to_change] == 1][a1 > a2] <- 0
    hd1[to_change][gt3[to_change] == 1][a1 <= a2] <- 0
    hd2[to_change][gt3[to_change] == 1][a1 <= a2] <- 1
    return(list(hd1 = hd1, hd2 = hd2))
}


## this is all about the phased output
## so can re-do them all 
recast_nipt_haps <- function(
    hap1,
    hap2,
    hap3,
    mat_gp_t,
    fet_gp_t
) {
    ## base almost entirely on argmax genotypes
    ## obvious option: do both, keep track of those changed twice
    gt1A <- round(hap1) + round(hap2)
    gt1B <- round(hap1) + round(hap3)    
    gtMT <- rep(0, ncol(mat_gp_t)) ## genotype mother thresolded
    gtFT <- rep(0, ncol(mat_gp_t)) ## genotype fetus thresholded
    max_valA <- mat_gp_t[1, ]
    max_valB <- fet_gp_t[1, ]
    for(i in 2:3) {
        ##
        w <- mat_gp_t[i, ] > max_valA
        gtMT[w] <- i - 1
        max_valA[w] <- mat_gp_t[i, w]
        ##
        w <- fet_gp_t[i, ] > max_valB
        gtFT[w] <- i - 1
        max_valB[w] <- fet_gp_t[i, w]
    }
    ## OK let's rock and roll
    ## here first two cols are genotypes, last three are haps, in the normal order
    ## three exceptions:
    ##   two that shouldn't happen:
    ##     mat genotypes 0 and fet genotypes 2
    ##     mat genotypes 2 and fet genotypes 0
    ##   one that definitely could happen
    ##     mat genotype 1 and fet genotype 1
    ##   in all cases, favour the maternal one, as it has more data almost certainly
    conv_mat <- rbind(
        c(0, 0, 0, 0, 0),
        c(0, 1, 0, 0, 1),
        c(0, 2, 0, 0, 1),
        c(1, 0, 0, 1, 0),
        c(1, 2, 1, 0, 1),
        c(2, 0, 1, 1, 0),
        c(2, 1, 1, 1, 0),
        c(2, 2, 1, 1, 1)
    )
    for(i_row in 1:nrow(conv_mat)) {
        w <- (gtMT == conv_mat[i_row, 1]) & (gtFT == conv_mat[i_row, 2])
        hap1[w] <- conv_mat[i_row, 3]
        hap2[w] <- conv_mat[i_row, 4]
        hap3[w] <- conv_mat[i_row, 5]        
    }
    ## 
    ## for double 1, should either be 1 0 0 or 0 1 1
    ## so when those work, assign them
    w1 <- (gtMT == 1) & (gtFT == 1)
    ## 
    w2 <- round(hap1[w1]) == 1 & round(hap2[w1]) == 0 & round(hap3[w1]) == 0
    hap1[w1][w2] <- 1
    hap2[w1][w2] <- 0
    hap3[w1][w2] <- 0
    ##
    w3 <- round(hap1[w1]) == 0 & round(hap2[w1]) == 1 & round(hap3[w1]) == 1
    hap1[w1][w3] <- 0
    hap2[w1][w3] <- 1
    hap3[w1][w3] <- 1
    ## otherwise, favour the maternal one, assume error in fetal
    w4 <- !w2 & !w3
    hap1[w1][w4] <- round(hap1[w1][w4])
    hap2[w1][w4] <- round(hap2[w1][w4])
    hap3[w1][w4] <- 1 - hap1[w1][w4]
    ## now round them all just in case
    hap1 <- round(hap1)
    hap2 <- round(hap2)
    hap3 <- round(hap3)
    return(list(hap1 = hap1, hap2 = hap2, hap3 = hap3))
}


## chunker

#' @export
quilt_chunk_map <- function(chr, genetic_map_file, min.bp = 3e6, min.cm = 4, ex.cnt = 10) {

  qmap <- data.table::fread(genetic_map_file, data.table = F)

  ## min.cm <- 5
  ## min.bp <- 3000000
  ## ex.cnt <- 10 ## extra overlap number of sites for ligation
  min.buf <- 0  ## don't know it at the moment. assume 500 kb

  OUT <- matrix(ncol=2)
  i <- 0
  start <- 1
  end <- 1
  while(end < max(qmap[,1])){
    end <- start + min.bp + min.buf
    chunk <- subset(qmap, position >= start & position <= end)
    ## if no maps, extend first chunk further
    while(nrow(chunk)==0){ 
      end <- end + min.bp 
      chunk <- subset(qmap, position >= start & position <= end)
    }
    ## check if cM > min.cm
    while(diff(range(chunk[,3])) < min.cm) {
      end <- end + min.bp / 3  ## extended by 1/3
      chunk <- subset(qmap, position >= start & position <= end)
      if(chunk[nrow(chunk),1] == qmap[nrow(qmap),1]) break;
    }
    a <- matrix(range(chunk[,1]), ncol=2)
    OUT <- rbind(OUT, a)
    start <- chunk[nrow(chunk)-ex.cnt, 1]
    i <- i+1
  }

  OUT <- OUT[-1,] ## remove first row. NA!
  stopifnot(i==nrow(OUT))

  ## check if we should merge the last two rows
  if(OUT[i,2] - OUT[i-1,2] < min.bp / 3) {
    OUT[i-1, 2] <- OUT[i, 2]
    OUT <- OUT[-i, ]
  }

  ## fix first start 
  OUT[1, 1] <- 1
  ## pad more for the last row!
  OUT[nrow(OUT),2] <- OUT[nrow(OUT), 2] + 5000000

  rg <- apply(OUT, 1, function(o) paste0(chr,":",paste0(o, collapse = "-")))
  dat <- data.frame(chunk = 1:length(rg)-1, chr = rep(chr, length(rg)), region = rg)
  return(dat)
  ## write.table(dat, "chunk.txt", row.names=F, col.names = F, quote = F, sep = "\t")
}

